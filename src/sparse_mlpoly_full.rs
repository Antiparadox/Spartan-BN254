//! Full sparse multilinear polynomial implementation for SNARK mode
//! 
//! This module provides the complete machinery for R1CS evaluation proofs,
//! including memory checking via product circuits.
//!
//! Key components:
//! - SparseMatPolynomial: Sparse matrix as multilinear polynomial
//! - MultiSparseMatPolynomialAsDense: Dense representation for commitment
//! - AddrTimestamps: Memory access tracking for permutation arguments
//! - Derefs: Dereference values for memory operations
//! - Layers/ProductLayer: Product circuits for memory checking
//! - PolyEvalNetwork: Combines row and column layers
//! - HashLayerProof: Proves hash layer consistency
//! - ProductLayerProof: Proves product layer consistency
//! - SparseMatPolyEvalProof: Final sparse matrix evaluation proof

#![allow(clippy::too_many_arguments)]
#![allow(clippy::type_complexity)]
#![allow(non_snake_case)]

use crate::dense_mlpoly::{DensePolynomial, EqPolynomial, PolyCommitment, PolyCommitmentGens, PolyEvalProof};
use crate::errors::ProofVerifyError;
use crate::math::Math;
use crate::product_tree::{DotProductCircuit, ProductCircuit, ProductCircuitEvalProofBatched};
use crate::random::RandomTape;
use crate::scalar::Scalar;
use crate::transcript::{AppendToTranscript, ProofTranscript};
use merlin::Transcript;
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;

// ============================================================================
// Sparse Matrix Entry
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SparseMatEntry {
    pub row: usize,
    pub col: usize,
    pub val: Scalar,
}

impl SparseMatEntry {
    pub fn new(row: usize, col: usize, val: Scalar) -> Self {
        SparseMatEntry { row, col, val }
    }
}

// ============================================================================
// Sparse Matrix Polynomial
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SparseMatPolynomial {
    num_vars_x: usize,
    num_vars_y: usize,
    M: Vec<SparseMatEntry>,
}

impl SparseMatPolynomial {
    pub fn new(num_vars_x: usize, num_vars_y: usize, M: Vec<SparseMatEntry>) -> Self {
        Self { num_vars_x, num_vars_y, M }
    }

    pub fn get_num_nz_entries(&self) -> usize {
        self.M.len().next_power_of_two()
    }

    pub fn num_vars_x(&self) -> usize {
        self.num_vars_x
    }

    pub fn num_vars_y(&self) -> usize {
        self.num_vars_y
    }

    pub fn entries(&self) -> &[SparseMatEntry] {
        &self.M
    }

    fn sparse_to_dense_vecs(&self, N: usize) -> (Vec<usize>, Vec<usize>, Vec<Scalar>) {
        assert!(N >= self.get_num_nz_entries());
        let mut ops_row: Vec<usize> = vec![0; N];
        let mut ops_col: Vec<usize> = vec![0; N];
        let mut val: Vec<Scalar> = vec![Scalar::zero(); N];

        for i in 0..self.M.len() {
            ops_row[i] = self.M[i].row;
            ops_col[i] = self.M[i].col;
            val[i] = self.M[i].val;
        }
        (ops_row, ops_col, val)
    }

    pub fn evaluate_with_tables(&self, rx_evals: &[Scalar], ry_evals: &[Scalar]) -> Scalar {
        self.M
            .iter()
            .map(|entry| rx_evals[entry.row] * ry_evals[entry.col] * entry.val)
            .sum()
    }

    pub fn multi_evaluate(polys: &[&SparseMatPolynomial], rx: &[Scalar], ry: &[Scalar]) -> Vec<Scalar> {
        let eval_table_rx = EqPolynomial::new(rx.to_vec()).evals();
        let eval_table_ry = EqPolynomial::new(ry.to_vec()).evals();

        polys
            .iter()
            .map(|poly| poly.evaluate_with_tables(&eval_table_rx, &eval_table_ry))
            .collect()
    }

    fn multi_sparse_to_dense_rep(sparse_polys: &[&SparseMatPolynomial]) -> MultiSparseMatPolynomialAsDense {
        assert!(!sparse_polys.is_empty());
        for i in 1..sparse_polys.len() {
            assert_eq!(sparse_polys[i].num_vars_x, sparse_polys[0].num_vars_x);
            assert_eq!(sparse_polys[i].num_vars_y, sparse_polys[0].num_vars_y);
        }

        let N = sparse_polys
            .iter()
            .map(|poly| poly.get_num_nz_entries())
            .max()
            .unwrap();

        let mut ops_row_vec: Vec<Vec<usize>> = Vec::new();
        let mut ops_col_vec: Vec<Vec<usize>> = Vec::new();
        let mut val_vec: Vec<DensePolynomial> = Vec::new();
        
        for poly in sparse_polys {
            let (ops_row, ops_col, val) = poly.sparse_to_dense_vecs(N);
            ops_row_vec.push(ops_row);
            ops_col_vec.push(ops_col);
            val_vec.push(DensePolynomial::new(val));
        }

        let any_poly = &sparse_polys[0];
        let num_mem_cells = if any_poly.num_vars_x > any_poly.num_vars_y {
            any_poly.num_vars_x.pow2()
        } else {
            any_poly.num_vars_y.pow2()
        };

        let row = AddrTimestamps::new(num_mem_cells, N, ops_row_vec);
        let col = AddrTimestamps::new(num_mem_cells, N, ops_col_vec);

        // Combine polynomials into a single polynomial for commitment purposes
        let comb_ops = DensePolynomial::merge(
            row.ops_addr
                .iter()
                .chain(row.read_ts.iter())
                .chain(col.ops_addr.iter())
                .chain(col.read_ts.iter())
                .chain(val_vec.iter()),
        );
        let mut comb_mem = row.audit_ts.clone();
        comb_mem.extend(&col.audit_ts);

        MultiSparseMatPolynomialAsDense {
            batch_size: sparse_polys.len(),
            row,
            col,
            val: val_vec,
            comb_ops,
            comb_mem,
        }
    }

    pub fn multi_commit(
        sparse_polys: &[&SparseMatPolynomial],
        gens: &SparseMatPolyCommitmentGens,
    ) -> (SparseMatPolyCommitment, MultiSparseMatPolynomialAsDense) {
        let batch_size = sparse_polys.len();
        let dense = SparseMatPolynomial::multi_sparse_to_dense_rep(sparse_polys);

        let (comm_comb_ops, _blinds_comb_ops) = dense.comb_ops.commit(&gens.gens_ops, None);
        let (comm_comb_mem, _blinds_comb_mem) = dense.comb_mem.commit(&gens.gens_mem, None);

        (
            SparseMatPolyCommitment {
                batch_size,
                num_mem_cells: dense.row.audit_ts.len(),
                num_ops: dense.row.read_ts[0].len(),
                comm_comb_ops,
                comm_comb_mem,
            },
            dense,
        )
    }
}

// ============================================================================
// Address Timestamps for Memory Checking
// ============================================================================

#[derive(Clone, Serialize, Deserialize)]
pub struct AddrTimestamps {
    ops_addr_usize: Vec<Vec<usize>>,
    ops_addr: Vec<DensePolynomial>,
    read_ts: Vec<DensePolynomial>,
    audit_ts: DensePolynomial,
}

impl AddrTimestamps {
    pub fn new(num_cells: usize, num_ops: usize, ops_addr: Vec<Vec<usize>>) -> Self {
        for item in ops_addr.iter() {
            assert_eq!(item.len(), num_ops);
        }

        let mut audit_ts = vec![0usize; num_cells];
        let mut ops_addr_vec: Vec<DensePolynomial> = Vec::new();
        let mut read_ts_vec: Vec<DensePolynomial> = Vec::new();
        
        for ops_addr_inst in ops_addr.iter() {
            let mut read_ts = vec![0usize; num_ops];

            for i in 0..num_ops {
                let addr = ops_addr_inst[i];
                assert!(addr < num_cells);
                let r_ts = audit_ts[addr];
                read_ts[i] = r_ts;
                let w_ts = r_ts + 1;
                audit_ts[addr] = w_ts;
            }

            ops_addr_vec.push(DensePolynomial::from_usize(ops_addr_inst));
            read_ts_vec.push(DensePolynomial::from_usize(&read_ts));
        }

        AddrTimestamps {
            ops_addr: ops_addr_vec,
            ops_addr_usize: ops_addr,
            read_ts: read_ts_vec,
            audit_ts: DensePolynomial::from_usize(&audit_ts),
        }
    }

    fn deref_mem(addr: &[usize], mem_val: &[Scalar]) -> DensePolynomial {
        DensePolynomial::new(
            (0..addr.len())
                .map(|i| mem_val[addr[i]])
                .collect(),
        )
    }

    pub fn deref(&self, mem_val: &[Scalar]) -> Vec<DensePolynomial> {
        (0..self.ops_addr.len())
            .map(|i| AddrTimestamps::deref_mem(&self.ops_addr_usize[i], mem_val))
            .collect()
    }
}

// ============================================================================
// Dense Representation of Multiple Sparse Matrices
// ============================================================================

#[derive(Clone, Serialize, Deserialize)]
pub struct MultiSparseMatPolynomialAsDense {
    batch_size: usize,
    val: Vec<DensePolynomial>,
    row: AddrTimestamps,
    col: AddrTimestamps,
    comb_ops: DensePolynomial,
    comb_mem: DensePolynomial,
}

impl MultiSparseMatPolynomialAsDense {
    pub fn deref(&self, row_mem_val: &[Scalar], col_mem_val: &[Scalar]) -> Derefs {
        let row_ops_val = self.row.deref(row_mem_val);
        let col_ops_val = self.col.deref(col_mem_val);
        Derefs::new(row_ops_val, col_ops_val)
    }
}

// ============================================================================
// Derefs - Dereference values from memory
// ============================================================================

pub struct Derefs {
    row_ops_val: Vec<DensePolynomial>,
    col_ops_val: Vec<DensePolynomial>,
    comb: DensePolynomial,
}

impl Derefs {
    pub fn new(row_ops_val: Vec<DensePolynomial>, col_ops_val: Vec<DensePolynomial>) -> Self {
        assert_eq!(row_ops_val.len(), col_ops_val.len());
        let comb = DensePolynomial::merge(row_ops_val.iter().chain(col_ops_val.iter()));
        Derefs { row_ops_val, col_ops_val, comb }
    }

    pub fn commit(&self, gens: &PolyCommitmentGens) -> DerefsCommitment {
        let (comm_ops_val, _blinds) = self.comb.commit(gens, None);
        DerefsCommitment { comm_ops_val }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DerefsCommitment {
    comm_ops_val: PolyCommitment,
}

impl AppendToTranscript for DerefsCommitment {
    fn append_to_transcript(&self, label: &'static [u8], transcript: &mut Transcript) {
        transcript.append_message(b"derefs_commitment", b"begin_derefs_commitment");
        self.comm_ops_val.append_to_transcript(label, transcript);
        transcript.append_message(b"derefs_commitment", b"end_derefs_commitment");
    }
}

// ============================================================================
// DerefsEvalProof
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DerefsEvalProof {
    proof_derefs: PolyEvalProof,
}

impl DerefsEvalProof {
    fn protocol_name() -> &'static [u8] {
        b"Derefs evaluation proof"
    }

    fn prove_single(
        joint_poly: &DensePolynomial,
        r: &[Scalar],
        evals: Vec<Scalar>,
        gens: &PolyCommitmentGens,
        transcript: &mut Transcript,
        random_tape: &mut RandomTape,
    ) -> PolyEvalProof {
        assert_eq!(joint_poly.get_num_vars(), r.len() + evals.len().log_2());

        evals.append_to_transcript(b"evals_ops_val", transcript);

        // n-to-1 reduction
        let challenges = transcript.challenge_vector(b"challenge_combine_n_to_one", evals.len().log_2());
        let mut poly_evals = DensePolynomial::new(evals);
        for i in (0..challenges.len()).rev() {
            poly_evals.bound_poly_var_bot(&challenges[i]);
        }
        assert_eq!(poly_evals.len(), 1);
        let joint_claim_eval = poly_evals[0];
        let mut r_joint = challenges;
        r_joint.extend(r);

        joint_claim_eval.append_to_transcript(b"joint_claim_eval", transcript);
        let (proof_derefs, _comm_derefs_eval) = PolyEvalProof::prove(
            joint_poly,
            None,
            &r_joint,
            &joint_claim_eval,
            None,
            gens,
            transcript,
            random_tape,
        );

        proof_derefs
    }

    pub fn prove(
        derefs: &Derefs,
        eval_row_ops_val_vec: &[Scalar],
        eval_col_ops_val_vec: &[Scalar],
        r: &[Scalar],
        gens: &PolyCommitmentGens,
        transcript: &mut Transcript,
        random_tape: &mut RandomTape,
    ) -> Self {
        transcript.append_protocol_name(DerefsEvalProof::protocol_name());

        let evals = {
            let mut evals = eval_row_ops_val_vec.to_owned();
            evals.extend(eval_col_ops_val_vec);
            evals.resize(evals.len().next_power_of_two(), Scalar::zero());
            evals
        };
        let proof_derefs = DerefsEvalProof::prove_single(&derefs.comb, r, evals, gens, transcript, random_tape);

        DerefsEvalProof { proof_derefs }
    }

    fn verify_single(
        proof: &PolyEvalProof,
        comm: &PolyCommitment,
        r: &[Scalar],
        evals: Vec<Scalar>,
        gens: &PolyCommitmentGens,
        transcript: &mut Transcript,
    ) -> Result<(), ProofVerifyError> {
        evals.append_to_transcript(b"evals_ops_val", transcript);

        let challenges = transcript.challenge_vector(b"challenge_combine_n_to_one", evals.len().log_2());
        let mut poly_evals = DensePolynomial::new(evals);
        for i in (0..challenges.len()).rev() {
            poly_evals.bound_poly_var_bot(&challenges[i]);
        }
        assert_eq!(poly_evals.len(), 1);
        let joint_claim_eval = poly_evals[0];
        let mut r_joint = challenges;
        r_joint.extend(r);

        joint_claim_eval.append_to_transcript(b"joint_claim_eval", transcript);

        proof.verify_plain(gens, transcript, &r_joint, &joint_claim_eval, comm)
    }

    pub fn verify(
        &self,
        r: &[Scalar],
        eval_row_ops_val_vec: &[Scalar],
        eval_col_ops_val_vec: &[Scalar],
        gens: &PolyCommitmentGens,
        comm: &DerefsCommitment,
        transcript: &mut Transcript,
    ) -> Result<(), ProofVerifyError> {
        transcript.append_protocol_name(DerefsEvalProof::protocol_name());
        let mut evals = eval_row_ops_val_vec.to_owned();
        evals.extend(eval_col_ops_val_vec);
        evals.resize(evals.len().next_power_of_two(), Scalar::zero());

        DerefsEvalProof::verify_single(
            &self.proof_derefs,
            &comm.comm_ops_val,
            r,
            evals,
            gens,
            transcript,
        )
    }
}

// ============================================================================
// Commitment Generators for Sparse Matrix Polynomial
// ============================================================================

#[derive(Clone, Serialize, Deserialize)]
pub struct SparseMatPolyCommitmentGens {
    gens_ops: PolyCommitmentGens,
    gens_mem: PolyCommitmentGens,
    gens_derefs: PolyCommitmentGens,
}

impl SparseMatPolyCommitmentGens {
    pub fn new(
        label: &'static [u8],
        num_vars_x: usize,
        num_vars_y: usize,
        num_nz_entries: usize,
        batch_size: usize,
    ) -> Self {
        let num_vars_ops = num_nz_entries.next_power_of_two().log_2() 
            + (batch_size * 5).next_power_of_two().log_2();
        let num_vars_mem = std::cmp::max(num_vars_x, num_vars_y) + 1;
        let num_vars_derefs = num_nz_entries.next_power_of_two().log_2() 
            + (batch_size * 2).next_power_of_two().log_2();

        let gens_ops = PolyCommitmentGens::new(num_vars_ops, label);
        let gens_mem = PolyCommitmentGens::new(num_vars_mem, label);
        let gens_derefs = PolyCommitmentGens::new(num_vars_derefs, label);
        
        SparseMatPolyCommitmentGens { gens_ops, gens_mem, gens_derefs }
    }
}

// ============================================================================
// Sparse Matrix Polynomial Commitment
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SparseMatPolyCommitment {
    batch_size: usize,
    num_ops: usize,
    num_mem_cells: usize,
    comm_comb_ops: PolyCommitment,
    comm_comb_mem: PolyCommitment,
}

impl AppendToTranscript for SparseMatPolyCommitment {
    fn append_to_transcript(&self, _label: &'static [u8], transcript: &mut Transcript) {
        transcript.append_u64(b"batch_size", self.batch_size as u64);
        transcript.append_u64(b"num_ops", self.num_ops as u64);
        transcript.append_u64(b"num_mem_cells", self.num_mem_cells as u64);
        self.comm_comb_ops.append_to_transcript(b"comm_comb_ops", transcript);
        self.comm_comb_mem.append_to_transcript(b"comm_comb_mem", transcript);
    }
}

// ============================================================================
// Product Layer for Memory Checking
// ============================================================================

#[derive(Debug)]
struct ProductLayer {
    init: ProductCircuit,
    read_vec: Vec<ProductCircuit>,
    write_vec: Vec<ProductCircuit>,
    audit: ProductCircuit,
}

#[derive(Debug)]
struct Layers {
    prod_layer: ProductLayer,
}

impl Layers {
    fn build_hash_layer(
        eval_table: &[Scalar],
        addrs_vec: &[DensePolynomial],
        derefs_vec: &[DensePolynomial],
        read_ts_vec: &[DensePolynomial],
        audit_ts: &DensePolynomial,
        r_mem_check: &(Scalar, Scalar),
    ) -> (DensePolynomial, Vec<DensePolynomial>, Vec<DensePolynomial>, DensePolynomial) {
        let (r_hash, r_multiset_check) = r_mem_check;
        let r_hash_sqr = *r_hash * *r_hash;
        
        let hash_func = |addr: &Scalar, val: &Scalar, ts: &Scalar| -> Scalar {
            *ts * r_hash_sqr + *val * *r_hash + *addr
        };

        let num_mem_cells = eval_table.len();
        let poly_init_hashed = DensePolynomial::new(
            (0..num_mem_cells)
                .map(|i| hash_func(&Scalar::from_u64(i as u64), &eval_table[i], &Scalar::zero()) - *r_multiset_check)
                .collect(),
        );
        
        let poly_audit_hashed = DensePolynomial::new(
            (0..num_mem_cells)
                .map(|i| hash_func(&Scalar::from_u64(i as u64), &eval_table[i], &audit_ts[i]) - *r_multiset_check)
                .collect(),
        );

        let mut poly_read_hashed_vec: Vec<DensePolynomial> = Vec::new();
        let mut poly_write_hashed_vec: Vec<DensePolynomial> = Vec::new();
        
        for i in 0..addrs_vec.len() {
            let (addrs, derefs, read_ts) = (&addrs_vec[i], &derefs_vec[i], &read_ts_vec[i]);
            let num_ops = addrs.len();
            
            let poly_read_hashed = DensePolynomial::new(
                (0..num_ops)
                    .map(|j| hash_func(&addrs[j], &derefs[j], &read_ts[j]) - *r_multiset_check)
                    .collect(),
            );
            poly_read_hashed_vec.push(poly_read_hashed);

            let poly_write_hashed = DensePolynomial::new(
                (0..num_ops)
                    .map(|j| hash_func(&addrs[j], &derefs[j], &(read_ts[j] + Scalar::one())) - *r_multiset_check)
                    .collect(),
            );
            poly_write_hashed_vec.push(poly_write_hashed);
        }

        (poly_init_hashed, poly_read_hashed_vec, poly_write_hashed_vec, poly_audit_hashed)
    }

    pub fn new(
        eval_table: &[Scalar],
        addr_timestamps: &AddrTimestamps,
        poly_ops_val: &[DensePolynomial],
        r_mem_check: &(Scalar, Scalar),
    ) -> Self {
        let (poly_init_hashed, poly_read_hashed_vec, poly_write_hashed_vec, poly_audit_hashed) =
            Layers::build_hash_layer(
                eval_table,
                &addr_timestamps.ops_addr,
                poly_ops_val,
                &addr_timestamps.read_ts,
                &addr_timestamps.audit_ts,
                r_mem_check,
            );

        let prod_init = ProductCircuit::new(&poly_init_hashed);
        let prod_read_vec: Vec<ProductCircuit> = poly_read_hashed_vec
            .iter()
            .map(|p| ProductCircuit::new(p))
            .collect();
        let prod_write_vec: Vec<ProductCircuit> = poly_write_hashed_vec
            .iter()
            .map(|p| ProductCircuit::new(p))
            .collect();
        let prod_audit = ProductCircuit::new(&poly_audit_hashed);

        // Verify subset check: init × writes = reads × audit
        let hashed_writes: Scalar = prod_write_vec.iter().map(|p| p.evaluate()).product();
        let hashed_write_set = prod_init.evaluate() * hashed_writes;
        let hashed_reads: Scalar = prod_read_vec.iter().map(|p| p.evaluate()).product();
        let hashed_read_set = hashed_reads * prod_audit.evaluate();
        debug_assert_eq!(hashed_read_set, hashed_write_set);

        Layers {
            prod_layer: ProductLayer {
                init: prod_init,
                read_vec: prod_read_vec,
                write_vec: prod_write_vec,
                audit: prod_audit,
            },
        }
    }
}

// ============================================================================
// Polynomial Evaluation Network
// ============================================================================

#[derive(Debug)]
struct PolyEvalNetwork {
    row_layers: Layers,
    col_layers: Layers,
}

impl PolyEvalNetwork {
    pub fn new(
        dense: &MultiSparseMatPolynomialAsDense,
        derefs: &Derefs,
        mem_rx: &[Scalar],
        mem_ry: &[Scalar],
        r_mem_check: &(Scalar, Scalar),
    ) -> Self {
        let row_layers = Layers::new(mem_rx, &dense.row, &derefs.row_ops_val, r_mem_check);
        let col_layers = Layers::new(mem_ry, &dense.col, &derefs.col_ops_val, r_mem_check);

        PolyEvalNetwork { row_layers, col_layers }
    }
}

// ============================================================================
// Hash Layer Proof
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HashLayerProof {
    eval_row: (Vec<Scalar>, Vec<Scalar>, Scalar),
    eval_col: (Vec<Scalar>, Vec<Scalar>, Scalar),
    eval_val: Vec<Scalar>,
    eval_derefs: (Vec<Scalar>, Vec<Scalar>),
    proof_ops: PolyEvalProof,
    proof_mem: PolyEvalProof,
    proof_derefs: DerefsEvalProof,
}

impl HashLayerProof {
    fn protocol_name() -> &'static [u8] {
        b"Sparse polynomial hash layer proof"
    }

    fn prove_helper(
        rand: (&Vec<Scalar>, &Vec<Scalar>),
        addr_timestamps: &AddrTimestamps,
    ) -> (Vec<Scalar>, Vec<Scalar>, Scalar) {
        let (rand_mem, rand_ops) = rand;

        let eval_ops_addr_vec = addr_timestamps.ops_addr
            .iter()
            .map(|addr| addr.evaluate(rand_ops))
            .collect();

        let eval_read_ts_vec = addr_timestamps.read_ts
            .iter()
            .map(|ts| ts.evaluate(rand_ops))
            .collect();

        let eval_audit_ts = addr_timestamps.audit_ts.evaluate(rand_mem);

        (eval_ops_addr_vec, eval_read_ts_vec, eval_audit_ts)
    }

    pub fn prove(
        rand: (&Vec<Scalar>, &Vec<Scalar>),
        dense: &MultiSparseMatPolynomialAsDense,
        derefs: &Derefs,
        gens: &SparseMatPolyCommitmentGens,
        transcript: &mut Transcript,
        random_tape: &mut RandomTape,
    ) -> Self {
        transcript.append_protocol_name(HashLayerProof::protocol_name());

        let (rand_mem, rand_ops) = rand;

        // Prove derefs at rand_ops
        let eval_row_ops_val: Vec<Scalar> = derefs.row_ops_val
            .iter()
            .map(|row| row.evaluate(rand_ops))
            .collect();
        let eval_col_ops_val: Vec<Scalar> = derefs.col_ops_val
            .iter()
            .map(|col| col.evaluate(rand_ops))
            .collect();
        
        let proof_derefs = DerefsEvalProof::prove(
            derefs,
            &eval_row_ops_val,
            &eval_col_ops_val,
            rand_ops,
            &gens.gens_derefs,
            transcript,
            random_tape,
        );
        let eval_derefs = (eval_row_ops_val, eval_col_ops_val);

        // Evaluate components
        let (eval_row_addr_vec, eval_row_read_ts_vec, eval_row_audit_ts) =
            HashLayerProof::prove_helper((rand_mem, rand_ops), &dense.row);
        let (eval_col_addr_vec, eval_col_read_ts_vec, eval_col_audit_ts) =
            HashLayerProof::prove_helper((rand_mem, rand_ops), &dense.col);
        let eval_val_vec: Vec<Scalar> = dense.val
            .iter()
            .map(|v| v.evaluate(rand_ops))
            .collect();

        // Combine ops evaluations for single proof
        let mut evals_ops: Vec<Scalar> = Vec::new();
        evals_ops.extend(&eval_row_addr_vec);
        evals_ops.extend(&eval_row_read_ts_vec);
        evals_ops.extend(&eval_col_addr_vec);
        evals_ops.extend(&eval_col_read_ts_vec);
        evals_ops.extend(&eval_val_vec);
        evals_ops.resize(evals_ops.len().next_power_of_two(), Scalar::zero());
        evals_ops.append_to_transcript(b"claim_evals_ops", transcript);
        
        let challenges_ops = transcript.challenge_vector(b"challenge_combine_n_to_one", evals_ops.len().log_2());

        let mut poly_evals_ops = DensePolynomial::new(evals_ops);
        for i in (0..challenges_ops.len()).rev() {
            poly_evals_ops.bound_poly_var_bot(&challenges_ops[i]);
        }
        assert_eq!(poly_evals_ops.len(), 1);
        let joint_claim_eval_ops = poly_evals_ops[0];
        let mut r_joint_ops = challenges_ops;
        r_joint_ops.extend(rand_ops);
        
        joint_claim_eval_ops.append_to_transcript(b"joint_claim_eval_ops", transcript);
        let (proof_ops, _) = PolyEvalProof::prove(
            &dense.comb_ops,
            None,
            &r_joint_ops,
            &joint_claim_eval_ops,
            None,
            &gens.gens_ops,
            transcript,
            random_tape,
        );

        // Combine mem evaluations for single proof
        let evals_mem: Vec<Scalar> = vec![eval_row_audit_ts, eval_col_audit_ts];
        evals_mem.append_to_transcript(b"claim_evals_mem", transcript);
        let challenges_mem = transcript.challenge_vector(b"challenge_combine_two_to_one", evals_mem.len().log_2());

        let mut poly_evals_mem = DensePolynomial::new(evals_mem);
        for i in (0..challenges_mem.len()).rev() {
            poly_evals_mem.bound_poly_var_bot(&challenges_mem[i]);
        }
        assert_eq!(poly_evals_mem.len(), 1);
        let joint_claim_eval_mem = poly_evals_mem[0];
        let mut r_joint_mem = challenges_mem;
        r_joint_mem.extend(rand_mem);
        
        joint_claim_eval_mem.append_to_transcript(b"joint_claim_eval_mem", transcript);
        let (proof_mem, _) = PolyEvalProof::prove(
            &dense.comb_mem,
            None,
            &r_joint_mem,
            &joint_claim_eval_mem,
            None,
            &gens.gens_mem,
            transcript,
            random_tape,
        );

        HashLayerProof {
            eval_row: (eval_row_addr_vec, eval_row_read_ts_vec, eval_row_audit_ts),
            eval_col: (eval_col_addr_vec, eval_col_read_ts_vec, eval_col_audit_ts),
            eval_val: eval_val_vec,
            eval_derefs,
            proof_ops,
            proof_mem,
            proof_derefs,
        }
    }

    fn verify_helper(
        rand: &(&Vec<Scalar>, &Vec<Scalar>),
        claims: &(Scalar, Vec<Scalar>, Vec<Scalar>, Scalar),
        eval_ops_val: &[Scalar],
        eval_ops_addr: &[Scalar],
        eval_read_ts: &[Scalar],
        eval_audit_ts: &Scalar,
        r: &[Scalar],
        r_hash: &Scalar,
        r_multiset_check: &Scalar,
    ) -> Result<(), ProofVerifyError> {
        let r_hash_sqr = *r_hash * *r_hash;
        let hash_func = |addr: &Scalar, val: &Scalar, ts: &Scalar| -> Scalar {
            *ts * r_hash_sqr + *val * *r_hash + *addr
        };

        let (rand_mem, _rand_ops) = rand;
        let (claim_init, claim_read, claim_write, claim_audit) = claims;

        // Verify init claim
        let eval_init_addr = IdentityPolynomial::new(rand_mem.len()).evaluate(rand_mem);
        let eval_init_val = EqPolynomial::new(r.to_vec()).evaluate(rand_mem);
        let hash_init_at_rand_mem = hash_func(&eval_init_addr, &eval_init_val, &Scalar::zero()) - *r_multiset_check;
        
        if *claim_init != hash_init_at_rand_mem {
            eprintln!("DEBUG: verify_helper failed at init claim check");
            eprintln!("  claim_init = {:?}", claim_init);
            eprintln!("  hash_init_at_rand_mem = {:?}", hash_init_at_rand_mem);
            return Err(ProofVerifyError::InternalError);
        }

        // Verify audit claim
        let hash_audit_at_rand_mem = hash_func(&eval_init_addr, &eval_init_val, eval_audit_ts) - *r_multiset_check;
        if *claim_audit != hash_audit_at_rand_mem {
            eprintln!("DEBUG: verify_helper failed at audit claim check");
            eprintln!("  claim_audit = {:?}", claim_audit);
            eprintln!("  hash_audit_at_rand_mem = {:?}", hash_audit_at_rand_mem);
            return Err(ProofVerifyError::InternalError);
        }

        // Verify read and write claims
        for i in 0..eval_ops_val.len() {
            let hash_read_at_rand_ops = hash_func(&eval_ops_addr[i], &eval_ops_val[i], &eval_read_ts[i]) - *r_multiset_check;
            if claim_read[i] != hash_read_at_rand_ops {
                eprintln!("DEBUG: verify_helper failed at read claim {} check", i);
                eprintln!("  claim_read[{}] = {:?}", i, claim_read[i]);
                eprintln!("  hash_read_at_rand_ops = {:?}", hash_read_at_rand_ops);
                return Err(ProofVerifyError::InternalError);
            }

            let hash_write_at_rand_ops = hash_func(
                &eval_ops_addr[i],
                &eval_ops_val[i],
                &(eval_read_ts[i] + Scalar::one()),
            ) - *r_multiset_check;
            if claim_write[i] != hash_write_at_rand_ops {
                eprintln!("DEBUG: verify_helper failed at write claim {} check", i);
                eprintln!("  claim_write[{}] = {:?}", i, claim_write[i]);
                eprintln!("  hash_write_at_rand_ops = {:?}", hash_write_at_rand_ops);
                return Err(ProofVerifyError::InternalError);
            }
        }

        Ok(())
    }

    pub fn verify(
        &self,
        rand: (&Vec<Scalar>, &Vec<Scalar>),
        claims_row: &(Scalar, Vec<Scalar>, Vec<Scalar>, Scalar),
        claims_col: &(Scalar, Vec<Scalar>, Vec<Scalar>, Scalar),
        claims_dotp: &[Scalar],
        comm: &SparseMatPolyCommitment,
        comm_derefs: &DerefsCommitment,
        gens: &SparseMatPolyCommitmentGens,
        rx: &[Scalar],
        ry: &[Scalar],
        r_hash: &Scalar,
        r_multiset_check: &Scalar,
        transcript: &mut Transcript,
    ) -> Result<(), ProofVerifyError> {
        transcript.append_protocol_name(HashLayerProof::protocol_name());

        let (rand_mem, rand_ops) = rand;
        let (eval_row_ops_val, eval_col_ops_val) = &self.eval_derefs;

        // Verify derefs proof
        self.proof_derefs.verify(
            rand_ops,
            eval_row_ops_val,
            eval_col_ops_val,
            &gens.gens_derefs,
            comm_derefs,
            transcript,
        )?;

        let (eval_row_addr_vec, eval_row_read_ts_vec, eval_row_audit_ts) = &self.eval_row;
        let (eval_col_addr_vec, eval_col_read_ts_vec, eval_col_audit_ts) = &self.eval_col;

        // Verify row claims
        HashLayerProof::verify_helper(
            &(rand_mem, rand_ops),
            claims_row,
            eval_row_ops_val,
            eval_row_addr_vec,
            eval_row_read_ts_vec,
            eval_row_audit_ts,
            rx,
            r_hash,
            r_multiset_check,
        )?;

        // Verify col claims
        HashLayerProof::verify_helper(
            &(rand_mem, rand_ops),
            claims_col,
            eval_col_ops_val,
            eval_col_addr_vec,
            eval_col_read_ts_vec,
            eval_col_audit_ts,
            ry,
            r_hash,
            r_multiset_check,
        )?;

        // Verify dotp claims
        // claims_dotp contains triples: [left_0, right_0, weight_0, left_1, right_1, weight_1, ...]
        // We verify that these match the evaluations from the derefs proof
        let num_instances = eval_row_ops_val.len();
        assert_eq!(claims_dotp.len(), 3 * num_instances);
        
        for i in 0..num_instances {
            let claim_left = claims_dotp[3 * i];
            let claim_right = claims_dotp[3 * i + 1];
            let claim_weight = claims_dotp[3 * i + 2];
            
            if claim_left != eval_row_ops_val[i] {
                eprintln!("DEBUG: HashLayerProof::verify failed at dotp left claim {} check", i);
                eprintln!("  claim_left = {:?}", claim_left);
                eprintln!("  eval_row_ops_val[{}] = {:?}", i, eval_row_ops_val[i]);
                return Err(ProofVerifyError::InternalError);
            }
            if claim_right != eval_col_ops_val[i] {
                eprintln!("DEBUG: HashLayerProof::verify failed at dotp right claim {} check", i);
                eprintln!("  claim_right = {:?}", claim_right);
                eprintln!("  eval_col_ops_val[{}] = {:?}", i, eval_col_ops_val[i]);
                return Err(ProofVerifyError::InternalError);
            }
            if claim_weight != self.eval_val[i] {
                eprintln!("DEBUG: HashLayerProof::verify failed at dotp weight claim {} check", i);
                eprintln!("  claim_weight = {:?}", claim_weight);
                eprintln!("  self.eval_val[{}] = {:?}", i, self.eval_val[i]);
                return Err(ProofVerifyError::InternalError);
            }
        }

        // Verify ops proof
        let mut evals_ops: Vec<Scalar> = Vec::new();
        evals_ops.extend(eval_row_addr_vec);
        evals_ops.extend(eval_row_read_ts_vec);
        evals_ops.extend(eval_col_addr_vec);
        evals_ops.extend(eval_col_read_ts_vec);
        evals_ops.extend(&self.eval_val);
        evals_ops.resize(evals_ops.len().next_power_of_two(), Scalar::zero());
        evals_ops.append_to_transcript(b"claim_evals_ops", transcript);

        let challenges_ops = transcript.challenge_vector(b"challenge_combine_n_to_one", evals_ops.len().log_2());
        let mut poly_evals_ops = DensePolynomial::new(evals_ops);
        for i in (0..challenges_ops.len()).rev() {
            poly_evals_ops.bound_poly_var_bot(&challenges_ops[i]);
        }
        let joint_claim_eval_ops = poly_evals_ops[0];
        let mut r_joint_ops = challenges_ops;
        r_joint_ops.extend(rand_ops);

        joint_claim_eval_ops.append_to_transcript(b"joint_claim_eval_ops", transcript);
        self.proof_ops.verify_plain(
            &gens.gens_ops,
            transcript,
            &r_joint_ops,
            &joint_claim_eval_ops,
            &comm.comm_comb_ops,
        )?;

        // Verify mem proof
        let evals_mem = vec![*eval_row_audit_ts, *eval_col_audit_ts];
        evals_mem.append_to_transcript(b"claim_evals_mem", transcript);
        let challenges_mem = transcript.challenge_vector(b"challenge_combine_two_to_one", evals_mem.len().log_2());
        let mut poly_evals_mem = DensePolynomial::new(evals_mem);
        for i in (0..challenges_mem.len()).rev() {
            poly_evals_mem.bound_poly_var_bot(&challenges_mem[i]);
        }
        let joint_claim_eval_mem = poly_evals_mem[0];
        let mut r_joint_mem = challenges_mem;
        r_joint_mem.extend(rand_mem);

        joint_claim_eval_mem.append_to_transcript(b"joint_claim_eval_mem", transcript);
        self.proof_mem.verify_plain(
            &gens.gens_mem,
            transcript,
            &r_joint_mem,
            &joint_claim_eval_mem,
            &comm.comm_comb_mem,
        )?;

        Ok(())
    }
}

// Placeholder for IdentityPolynomial
struct IdentityPolynomial {
    num_vars: usize,
}

impl IdentityPolynomial {
    fn new(num_vars: usize) -> Self {
        IdentityPolynomial { num_vars }
    }

    fn evaluate(&self, r: &[Scalar]) -> Scalar {
        let mut result = Scalar::zero();
        let base = Scalar::from_u64(2);
        for i in 0..r.len() {
            result = result * base + r[i];
        }
        result
    }
}

// ============================================================================
// Product Layer Proof  
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProductLayerProof {
    eval_row: (Scalar, Vec<Scalar>, Vec<Scalar>, Scalar),
    eval_col: (Scalar, Vec<Scalar>, Vec<Scalar>, Scalar),
    eval_val: (Vec<Scalar>, Vec<Scalar>),
    proof_mem: ProductCircuitEvalProofBatched,
    proof_ops: ProductCircuitEvalProofBatched,
}

impl ProductLayerProof {
    fn protocol_name() -> &'static [u8] {
        b"Sparse polynomial product layer proof"
    }

    pub fn prove(
        row_prod_layer: &mut ProductLayer,
        col_prod_layer: &mut ProductLayer,
        dense: &MultiSparseMatPolynomialAsDense,
        derefs: &Derefs,
        eval: &[Scalar],
        transcript: &mut Transcript,
    ) -> (Self, Vec<Scalar>, Vec<Scalar>) {
        transcript.append_protocol_name(ProductLayerProof::protocol_name());

        let row_eval_init = row_prod_layer.init.evaluate();
        let row_eval_audit = row_prod_layer.audit.evaluate();
        let row_eval_read: Vec<Scalar> = row_prod_layer.read_vec.iter().map(|p| p.evaluate()).collect();
        let row_eval_write: Vec<Scalar> = row_prod_layer.write_vec.iter().map(|p| p.evaluate()).collect();

        // Subset check
        let ws: Scalar = row_eval_write.iter().cloned().product();
        let rs: Scalar = row_eval_read.iter().cloned().product();
        assert_eq!(row_eval_init * ws, rs * row_eval_audit);

        row_eval_init.append_to_transcript(b"claim_row_eval_init", transcript);
        row_eval_read.append_to_transcript(b"claim_row_eval_read", transcript);
        row_eval_write.append_to_transcript(b"claim_row_eval_write", transcript);
        row_eval_audit.append_to_transcript(b"claim_row_eval_audit", transcript);

        let col_eval_init = col_prod_layer.init.evaluate();
        let col_eval_audit = col_prod_layer.audit.evaluate();
        let col_eval_read: Vec<Scalar> = col_prod_layer.read_vec.iter().map(|p| p.evaluate()).collect();
        let col_eval_write: Vec<Scalar> = col_prod_layer.write_vec.iter().map(|p| p.evaluate()).collect();

        // Subset check
        let ws: Scalar = col_eval_write.iter().cloned().product();
        let rs: Scalar = col_eval_read.iter().cloned().product();
        assert_eq!(col_eval_init * ws, rs * col_eval_audit);

        col_eval_init.append_to_transcript(b"claim_col_eval_init", transcript);
        col_eval_read.append_to_transcript(b"claim_col_eval_read", transcript);
        col_eval_write.append_to_transcript(b"claim_col_eval_write", transcript);
        col_eval_audit.append_to_transcript(b"claim_col_eval_audit", transcript);

        // Prepare dot product circuits for batching
        // We need to create interleaved circuits for verifier compatibility:
        // (left_0, right_0, left_1, right_1, ...) ordering
        assert_eq!(eval.len(), derefs.row_ops_val.len());
        let mut dotp_circuits_owned: Vec<DotProductCircuit> = Vec::new();
        let mut eval_dotp_left_vec: Vec<Scalar> = Vec::new();
        let mut eval_dotp_right_vec: Vec<Scalar> = Vec::new();

        for i in 0..derefs.row_ops_val.len() {
            let left = derefs.row_ops_val[i].clone();
            let right = derefs.col_ops_val[i].clone();
            let weights = dense.val[i].clone();

            let mut dotp_circuit = DotProductCircuit::new(left, right, weights);
            let (dotp_circuit_left, dotp_circuit_right) = dotp_circuit.split();

            let (eval_dotp_left, eval_dotp_right) = (dotp_circuit_left.evaluate(), dotp_circuit_right.evaluate());

            eval_dotp_left.append_to_transcript(b"claim_eval_dotp_left", transcript);
            eval_dotp_right.append_to_transcript(b"claim_eval_dotp_right", transcript);
            assert_eq!(eval_dotp_left + eval_dotp_right, eval[i]);
            
            eval_dotp_left_vec.push(eval_dotp_left);
            eval_dotp_right_vec.push(eval_dotp_right);
            // Add in interleaved order: left_i, right_i
            dotp_circuits_owned.push(dotp_circuit_left);
            dotp_circuits_owned.push(dotp_circuit_right);
        }

        // For simplicity, we handle the general case
        // The original Spartan assumes batch_size = 3, but we'll be more flexible
        let num_instances = row_prod_layer.read_vec.len();
        
        // Collect all product circuits for ops proof
        let mut ops_circuits: Vec<&mut ProductCircuit> = Vec::new();
        for circ in row_prod_layer.read_vec.iter_mut() {
            ops_circuits.push(circ);
        }
        for circ in row_prod_layer.write_vec.iter_mut() {
            ops_circuits.push(circ);
        }
        for circ in col_prod_layer.read_vec.iter_mut() {
            ops_circuits.push(circ);
        }
        for circ in col_prod_layer.write_vec.iter_mut() {
            ops_circuits.push(circ);
        }

        // Collect dot product circuits
        let mut dotp_circuits: Vec<&mut DotProductCircuit> = dotp_circuits_owned.iter_mut().collect();

        let (proof_ops, rand_ops) = ProductCircuitEvalProofBatched::prove(
            &mut ops_circuits,
            &mut dotp_circuits,
            transcript,
        );

        // Prove mem circuits (init and audit for row and col)
        let mut mem_circuits: Vec<&mut ProductCircuit> = vec![
            &mut row_prod_layer.init,
            &mut row_prod_layer.audit,
            &mut col_prod_layer.init,
            &mut col_prod_layer.audit,
        ];

        let (proof_mem, rand_mem) = ProductCircuitEvalProofBatched::prove(
            &mut mem_circuits,
            &mut [],
            transcript,
        );

        (
            ProductLayerProof {
                eval_row: (row_eval_init, row_eval_read, row_eval_write, row_eval_audit),
                eval_col: (col_eval_init, col_eval_read, col_eval_write, col_eval_audit),
                eval_val: (eval_dotp_left_vec, eval_dotp_right_vec),
                proof_mem,
                proof_ops,
            },
            rand_mem,
            rand_ops,
        )
    }

    /// Returns:
    /// - claims_mem: reduced claims for mem circuits (row_init, row_audit, col_init, col_audit)
    /// - rand_mem: random challenges for mem
    /// - claims_ops: reduced claims for ops circuits (read/write for row and col)
    /// - claims_dotp: reduced claims for dotp circuits
    /// - rand_ops: random challenges for ops
    pub fn verify(
        &self,
        num_ops: usize,
        num_mem_cells: usize,
        evals: &[Scalar],
        transcript: &mut Transcript,
    ) -> Result<(Vec<Scalar>, Vec<Scalar>, Vec<Scalar>, Vec<Scalar>, Vec<Scalar>), ProofVerifyError> {
        transcript.append_protocol_name(ProductLayerProof::protocol_name());

        let num_instances = evals.len();
        let (row_eval_init, row_eval_read, row_eval_write, row_eval_audit) = &self.eval_row;
        let (col_eval_init, col_eval_read, col_eval_write, col_eval_audit) = &self.eval_col;
        let (eval_dotp_left_vec, eval_dotp_right_vec) = &self.eval_val;

        assert_eq!(row_eval_read.len(), num_instances);
        assert_eq!(row_eval_write.len(), num_instances);

        // Verify subset check for row
        let ws: Scalar = row_eval_write.iter().cloned().product();
        let rs: Scalar = row_eval_read.iter().cloned().product();
        if *row_eval_init * ws != rs * *row_eval_audit {
            eprintln!("DEBUG: ProductLayerProof::verify failed at row subset check");
            return Err(ProofVerifyError::InternalError);
        }

        row_eval_init.append_to_transcript(b"claim_row_eval_init", transcript);
        row_eval_read.append_to_transcript(b"claim_row_eval_read", transcript);
        row_eval_write.append_to_transcript(b"claim_row_eval_write", transcript);
        row_eval_audit.append_to_transcript(b"claim_row_eval_audit", transcript);

        // Verify subset check for col
        let ws: Scalar = col_eval_write.iter().cloned().product();
        let rs: Scalar = col_eval_read.iter().cloned().product();
        if *col_eval_init * ws != rs * *col_eval_audit {
            eprintln!("DEBUG: ProductLayerProof::verify failed at col subset check");
            return Err(ProofVerifyError::InternalError);
        }

        col_eval_init.append_to_transcript(b"claim_col_eval_init", transcript);
        col_eval_read.append_to_transcript(b"claim_col_eval_read", transcript);
        col_eval_write.append_to_transcript(b"claim_col_eval_write", transcript);
        col_eval_audit.append_to_transcript(b"claim_col_eval_audit", transcript);

        // Verify dotp claims
        let mut claims_dotp_circuit: Vec<Scalar> = Vec::new();
        for i in 0..num_instances {
            if eval_dotp_left_vec[i] + eval_dotp_right_vec[i] != evals[i] {
                eprintln!("DEBUG: ProductLayerProof::verify failed at dotp split check {}", i);
                return Err(ProofVerifyError::InternalError);
            }
            eval_dotp_left_vec[i].append_to_transcript(b"claim_eval_dotp_left", transcript);
            eval_dotp_right_vec[i].append_to_transcript(b"claim_eval_dotp_right", transcript);
            
            claims_dotp_circuit.push(eval_dotp_left_vec[i]);
            claims_dotp_circuit.push(eval_dotp_right_vec[i]);
        }

        // Build claims for ops proof
        let mut claims_prod_circuit: Vec<Scalar> = Vec::new();
        claims_prod_circuit.extend(row_eval_read);
        claims_prod_circuit.extend(row_eval_write);
        claims_prod_circuit.extend(col_eval_read);
        claims_prod_circuit.extend(col_eval_write);

        // Verify ops proof and get reduced claims
        let (claims_ops, claims_dotp, rand_ops) = self.proof_ops.verify(
            &claims_prod_circuit,
            &claims_dotp_circuit,
            num_ops,
            transcript,
        );

        // Build claims for mem proof
        let claims_prod_mem = vec![*row_eval_init, *row_eval_audit, *col_eval_init, *col_eval_audit];

        // Verify mem proof and get reduced claims
        let (claims_mem, _, rand_mem) = self.proof_mem.verify(
            &claims_prod_mem,
            &[],
            num_mem_cells,
            transcript,
        );

        Ok((claims_mem, rand_mem, claims_ops, claims_dotp, rand_ops))
    }
}

// ============================================================================
// Polynomial Evaluation Network Proof
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PolyEvalNetworkProof {
    proof_prod_layer: ProductLayerProof,
    proof_hash_layer: HashLayerProof,
}

impl PolyEvalNetworkProof {
    fn protocol_name() -> &'static [u8] {
        b"Sparse polynomial evaluation proof"
    }

    pub fn prove(
        network: &mut PolyEvalNetwork,
        dense: &MultiSparseMatPolynomialAsDense,
        derefs: &Derefs,
        evals: &[Scalar],
        gens: &SparseMatPolyCommitmentGens,
        transcript: &mut Transcript,
        random_tape: &mut RandomTape,
    ) -> Self {
        transcript.append_protocol_name(PolyEvalNetworkProof::protocol_name());

        let (proof_prod_layer, rand_mem, rand_ops) = ProductLayerProof::prove(
            &mut network.row_layers.prod_layer,
            &mut network.col_layers.prod_layer,
            dense,
            derefs,
            evals,
            transcript,
        );

        let proof_hash_layer = HashLayerProof::prove(
            (&rand_mem, &rand_ops),
            dense,
            derefs,
            gens,
            transcript,
            random_tape,
        );

        PolyEvalNetworkProof {
            proof_prod_layer,
            proof_hash_layer,
        }
    }

    pub fn verify(
        &self,
        comm: &SparseMatPolyCommitment,
        comm_derefs: &DerefsCommitment,
        evals: &[Scalar],
        gens: &SparseMatPolyCommitmentGens,
        rx: &[Scalar],
        ry: &[Scalar],
        r_mem_check: &(Scalar, Scalar),
        nz: usize,
        transcript: &mut Transcript,
    ) -> Result<(), ProofVerifyError> {
        transcript.append_protocol_name(PolyEvalNetworkProof::protocol_name());

        let num_instances = evals.len();
        let (r_hash, r_multiset_check) = r_mem_check;

        let num_ops = nz.next_power_of_two();
        let num_mem_cells = comm.num_mem_cells;

        // Verify product layer and get reduced claims
        let (claims_mem, rand_mem, claims_ops, claims_dotp, rand_ops) = self.proof_prod_layer.verify(
            num_ops,
            num_mem_cells,
            evals,
            transcript,
        )?;

        // claims_mem: [row_init, row_audit, col_init, col_audit] (4 elements)
        // claims_ops: [row_read_0..n, row_write_0..n, col_read_0..n, col_write_0..n] (4*n elements)
        // claims_dotp: [dotp_left_0, dotp_right_0, ...] (3*2*n = 6n elements for n instances)
        assert_eq!(claims_mem.len(), 4);
        assert_eq!(claims_ops.len(), 4 * num_instances);

        // Split claims_ops into row and col, read and write
        let (claims_ops_row, claims_ops_col) = claims_ops.split_at(2 * num_instances);
        let (claims_ops_row_read, claims_ops_row_write) = claims_ops_row.split_at(num_instances);
        let (claims_ops_col_read, claims_ops_col_write) = claims_ops_col.split_at(num_instances);

        // Build claims tuples for hash layer verification
        let claims_row = (
            claims_mem[0],
            claims_ops_row_read.to_vec(),
            claims_ops_row_write.to_vec(),
            claims_mem[1],
        );
        let claims_col = (
            claims_mem[2],
            claims_ops_col_read.to_vec(),
            claims_ops_col_write.to_vec(),
            claims_mem[3],
        );

        self.proof_hash_layer.verify(
            (&rand_mem, &rand_ops),
            &claims_row,
            &claims_col,
            &claims_dotp,
            comm,
            comm_derefs,
            gens,
            rx,
            ry,
            r_hash,
            r_multiset_check,
            transcript,
        )?;

        Ok(())
    }
}

// ============================================================================
// Sparse Matrix Polynomial Evaluation Proof (Final SNARK Component)
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SparseMatPolyEvalProof {
    comm_derefs: DerefsCommitment,
    poly_eval_network_proof: PolyEvalNetworkProof,
}

impl SparseMatPolyEvalProof {
    fn protocol_name() -> &'static [u8] {
        b"Sparse polynomial evaluation proof"
    }

    fn equalize(rx: &[Scalar], ry: &[Scalar]) -> (Vec<Scalar>, Vec<Scalar>) {
        match rx.len().cmp(&ry.len()) {
            Ordering::Less => {
                let diff = ry.len() - rx.len();
                let mut rx_ext = vec![Scalar::zero(); diff];
                rx_ext.extend(rx);
                (rx_ext, ry.to_vec())
            }
            Ordering::Greater => {
                let diff = rx.len() - ry.len();
                let mut ry_ext = vec![Scalar::zero(); diff];
                ry_ext.extend(ry);
                (rx.to_vec(), ry_ext)
            }
            Ordering::Equal => (rx.to_vec(), ry.to_vec()),
        }
    }

    pub fn prove(
        dense: &MultiSparseMatPolynomialAsDense,
        rx: &[Scalar],
        ry: &[Scalar],
        evals: &[Scalar],
        gens: &SparseMatPolyCommitmentGens,
        transcript: &mut Transcript,
        random_tape: &mut RandomTape,
    ) -> Self {
        transcript.append_protocol_name(SparseMatPolyEvalProof::protocol_name());

        assert_eq!(evals.len(), dense.batch_size);

        let (mem_rx, mem_ry) = {
            let (rx_ext, ry_ext) = SparseMatPolyEvalProof::equalize(rx, ry);
            let poly_rx = EqPolynomial::new(rx_ext).evals();
            let poly_ry = EqPolynomial::new(ry_ext).evals();
            (poly_rx, poly_ry)
        };

        let derefs = dense.deref(&mem_rx, &mem_ry);

        // Commit to derefs
        let comm_derefs = {
            let comm = derefs.commit(&gens.gens_derefs);
            comm.append_to_transcript(b"comm_poly_row_col_ops_val", transcript);
            comm
        };

        // Get random challenge for memory check
        let r_mem_check = transcript.challenge_vector(b"challenge_r_hash", 2);

        // Build evaluation network
        let mut net = PolyEvalNetwork::new(
            dense,
            &derefs,
            &mem_rx,
            &mem_ry,
            &(r_mem_check[0], r_mem_check[1]),
        );

        let poly_eval_network_proof = PolyEvalNetworkProof::prove(
            &mut net,
            dense,
            &derefs,
            evals,
            gens,
            transcript,
            random_tape,
        );

        SparseMatPolyEvalProof {
            comm_derefs,
            poly_eval_network_proof,
        }
    }

    pub fn verify(
        &self,
        comm: &SparseMatPolyCommitment,
        rx: &[Scalar],
        ry: &[Scalar],
        evals: &[Scalar],
        gens: &SparseMatPolyCommitmentGens,
        transcript: &mut Transcript,
    ) -> Result<(), ProofVerifyError> {
        transcript.append_protocol_name(SparseMatPolyEvalProof::protocol_name());

        let (rx_ext, ry_ext) = SparseMatPolyEvalProof::equalize(rx, ry);
        let (nz, num_mem_cells) = (comm.num_ops, comm.num_mem_cells);
        assert_eq!(rx_ext.len().pow2(), num_mem_cells);

        self.comm_derefs.append_to_transcript(b"comm_poly_row_col_ops_val", transcript);

        let r_mem_check = transcript.challenge_vector(b"challenge_r_hash", 2);

        self.poly_eval_network_proof.verify(
            comm,
            &self.comm_derefs,
            evals,
            gens,
            &rx_ext,
            &ry_ext,
            &(r_mem_check[0], r_mem_check[1]),
            nz,
            transcript,
        )
    }
}

// ============================================================================
// R1CS Evaluation Proof using Sparse Matrix Polynomial
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct R1CSEvalProofFull {
    proof: SparseMatPolyEvalProof,
}

impl R1CSEvalProofFull {
    pub fn prove(
        dense: &MultiSparseMatPolynomialAsDense,
        rx: &[Scalar],
        ry: &[Scalar],
        evals: &(Scalar, Scalar, Scalar),
        gens: &SparseMatPolyCommitmentGens,
        transcript: &mut Transcript,
        random_tape: &mut RandomTape,
    ) -> Self {
        let (eval_A, eval_B, eval_C) = evals;
        let evals_vec = vec![*eval_A, *eval_B, *eval_C];

        let proof = SparseMatPolyEvalProof::prove(
            dense,
            rx,
            ry,
            &evals_vec,
            gens,
            transcript,
            random_tape,
        );

        R1CSEvalProofFull { proof }
    }

    pub fn verify(
        &self,
        comm: &SparseMatPolyCommitment,
        rx: &[Scalar],
        ry: &[Scalar],
        evals: &(Scalar, Scalar, Scalar),
        gens: &SparseMatPolyCommitmentGens,
        transcript: &mut Transcript,
    ) -> Result<(), ProofVerifyError> {
        let (eval_A, eval_B, eval_C) = evals;
        let evals_vec = vec![*eval_A, *eval_B, *eval_C];

        self.proof.verify(comm, rx, ry, &evals_vec, gens, transcript)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sparse_mat_poly() {
        // Create a simple 4x4 sparse matrix
        let entries = vec![
            SparseMatEntry::new(0, 0, Scalar::one()),
            SparseMatEntry::new(1, 1, Scalar::from_u64(2)),
            SparseMatEntry::new(2, 2, Scalar::from_u64(3)),
            SparseMatEntry::new(3, 3, Scalar::from_u64(4)),
        ];

        let poly = SparseMatPolynomial::new(2, 2, entries);
        assert_eq!(poly.get_num_nz_entries(), 4);
    }
}

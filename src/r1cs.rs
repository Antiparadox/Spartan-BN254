//! R1CS constraint system for Spartan

use ark_serialize::{CanonicalSerialize, CanonicalDeserialize};
use crate::dense_mlpoly::DensePolynomial;
use crate::errors::{ProofVerifyError, R1CSError};
use crate::math::Math;
use crate::random::RandomTape;
use crate::scalar::Scalar;
use crate::sparse_mlpoly::{SparseMatEntry, SparseMatPolynomial};
use crate::sparse_mlpoly_full::{
    MultiSparseMatPolynomialAsDense, SparseMatPolyCommitment, SparseMatPolyCommitmentGens,
    SparseMatPolyEvalProof,
    SparseMatPolynomial as SparseMatPolynomialFull,
};
use crate::timer::Timer;
use crate::transcript::{AppendToTranscript, ProofTranscript};
use flate2::{write::ZlibEncoder, Compression};
use merlin::Transcript;
use serde::{Deserialize, Serialize};

/// R1CS constraint system shape (A, B, C matrices)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct R1CSShape {
    num_cons: usize,
    num_vars: usize,
    num_inputs: usize,
    pub A: SparseMatPolynomial,
    pub B: SparseMatPolynomial,
    pub C: SparseMatPolynomial,
}

impl R1CSShape {
    /// Create new R1CS shape from sparse matrix triplets
    pub fn new(
        num_cons: usize,
        num_vars: usize,
        num_inputs: usize,
        A: &[(usize, usize, Scalar)],
        B: &[(usize, usize, Scalar)],
        C: &[(usize, usize, Scalar)],
    ) -> Self {
        Timer::print(&format!("number_of_constraints {num_cons}"));
        Timer::print(&format!("number_of_variables {num_vars}"));
        Timer::print(&format!("number_of_inputs {num_inputs}"));
        Timer::print(&format!("number_non-zero_entries_A {}", A.len()));
        Timer::print(&format!("number_non-zero_entries_B {}", B.len()));
        Timer::print(&format!("number_non-zero_entries_C {}", C.len()));

        // Spartan requires power-of-2 dimensions
        assert_eq!(num_cons.next_power_of_two(), num_cons, "num_cons must be power of 2");
        assert_eq!(num_vars.next_power_of_two(), num_vars, "num_vars must be power of 2");
        assert!(num_inputs < num_vars, "num_inputs must be less than num_vars");

        let num_poly_vars_x = num_cons.log_2();
        let num_poly_vars_y = (2 * num_vars).log_2();

        let mat_A: Vec<SparseMatEntry> = A
            .iter()
            .map(|(row, col, val)| SparseMatEntry::new(*row, *col, *val))
            .collect();
        let mat_B: Vec<SparseMatEntry> = B
            .iter()
            .map(|(row, col, val)| SparseMatEntry::new(*row, *col, *val))
            .collect();
        let mat_C: Vec<SparseMatEntry> = C
            .iter()
            .map(|(row, col, val)| SparseMatEntry::new(*row, *col, *val))
            .collect();

        let poly_A = SparseMatPolynomial::new(num_poly_vars_x, num_poly_vars_y, mat_A);
        let poly_B = SparseMatPolynomial::new(num_poly_vars_x, num_poly_vars_y, mat_B);
        let poly_C = SparseMatPolynomial::new(num_poly_vars_x, num_poly_vars_y, mat_C);

        Self {
            num_cons,
            num_vars,
            num_inputs,
            A: poly_A,
            B: poly_B,
            C: poly_C,
        }
    }

    pub fn get_num_vars(&self) -> usize {
        self.num_vars
    }

    pub fn get_num_cons(&self) -> usize {
        self.num_cons
    }

    pub fn get_num_inputs(&self) -> usize {
        self.num_inputs
    }

    /// Compute digest of the R1CS shape
    pub fn get_digest(&self) -> Vec<u8> {
        let mut encoder = ZlibEncoder::new(Vec::new(), Compression::default());
        bincode::serialize_into(&mut encoder, &self).unwrap();
        encoder.finish().unwrap()
    }

    /// Check if the R1CS is satisfied by the given assignment
    /// z = (vars, 1, inputs) in Spartan format
    pub fn is_sat(&self, vars: &[Scalar], inputs: &[Scalar]) -> bool {
        assert_eq!(vars.len(), self.num_vars);
        assert_eq!(inputs.len(), self.num_inputs);

        // Build z = (vars, 1, inputs) in Spartan format
        let mut z = vars.to_vec();
        z.push(Scalar::one());
        z.extend(inputs);

        let num_cols = self.num_vars + self.num_inputs + 1;

        // Compute Az, Bz, Cz
        let Az = self.A.multiply_vec(self.num_cons, num_cols, &z);
        let Bz = self.B.multiply_vec(self.num_cons, num_cols, &z);
        let Cz = self.C.multiply_vec(self.num_cons, num_cols, &z);

        // Check Az ∘ Bz = Cz
        (0..self.num_cons).all(|i| Az[i] * Bz[i] == Cz[i])
    }

    /// Evaluate the R1CS polynomials at point (rx, ry)
    pub fn evaluate(&self, rx: &[Scalar], ry: &[Scalar]) -> (Scalar, Scalar, Scalar) {
        let evals = SparseMatPolynomial::multi_evaluate(&[&self.A, &self.B, &self.C], rx, ry);
        (evals[0], evals[1], evals[2])
    }

    /// Compute Az, Bz, Cz as dense polynomials
    pub fn multiply_vec(
        &self,
        num_rows: usize,
        num_cols: usize,
        z: &[Scalar],
    ) -> (DensePolynomial, DensePolynomial, DensePolynomial) {
        assert_eq!(num_rows, self.num_cons);
        assert_eq!(z.len(), num_cols);
        (
            DensePolynomial::new(self.A.multiply_vec(num_rows, num_cols, z)),
            DensePolynomial::new(self.B.multiply_vec(num_rows, num_cols, z)),
            DensePolynomial::new(self.C.multiply_vec(num_rows, num_cols, z)),
        )
    }

    /// Compute evaluation tables for sparse polynomials
    pub fn compute_eval_table_sparse(
        &self,
        num_rows: usize,
        num_cols: usize,
        evals: &[Scalar],
    ) -> (Vec<Scalar>, Vec<Scalar>, Vec<Scalar>) {
        assert_eq!(num_rows, self.num_cons);
        (
            self.A.compute_eval_table_sparse(evals, num_rows, num_cols),
            self.B.compute_eval_table_sparse(evals, num_rows, num_cols),
            self.C.compute_eval_table_sparse(evals, num_rows, num_cols),
        )
    }
}

/// R1CS Instance (wraps R1CSShape with digest)
#[derive(Clone)]
pub struct R1CSInstance {
    pub shape: R1CSShape,
    pub digest: Vec<u8>,
}

impl R1CSInstance {
    /// Create from raw triplets (for compatibility with original Spartan API)
    pub fn new(
        num_cons: usize,
        num_vars: usize,
        num_inputs: usize,
        A: &[(usize, usize, [u8; 32])],
        B: &[(usize, usize, [u8; 32])],
        C: &[(usize, usize, [u8; 32])],
    ) -> Result<Self, R1CSError> {
        // Pad to power of 2 if needed
        let num_vars_padded = std::cmp::max(num_vars, num_inputs + 1).next_power_of_two();
        let num_cons_padded = std::cmp::max(num_cons, 2).next_power_of_two();

        // Convert byte arrays to Scalars
        let convert = |triplets: &[(usize, usize, [u8; 32])], pad_offset: usize| -> Result<Vec<(usize, usize, Scalar)>, R1CSError> {
            triplets
                .iter()
                .map(|(row, col, val_bytes)| {
                    let scalar = Scalar::from_bytes(val_bytes)
                        .ok_or(R1CSError::InvalidScalar)?;
                    // Adjust column index for padding
                    let adjusted_col = if *col >= num_vars {
                        *col + pad_offset
                    } else {
                        *col
                    };
                    Ok((*row, adjusted_col, scalar))
                })
                .collect()
        };

        let pad_offset = num_vars_padded - num_vars;
        let A_scalar = convert(A, pad_offset)?;
        let B_scalar = convert(B, pad_offset)?;
        let C_scalar = convert(C, pad_offset)?;

        let shape = R1CSShape::new(
            num_cons_padded,
            num_vars_padded,
            num_inputs,
            &A_scalar,
            &B_scalar,
            &C_scalar,
        );

        let digest = shape.get_digest();

        Ok(R1CSInstance { shape, digest })
    }

    /// Create from Scalars directly (for use with R1CS reader)
    pub fn from_scalars(
        num_cons: usize,
        num_vars: usize,
        num_inputs: usize,
        A: Vec<(usize, usize, Scalar)>,
        B: Vec<(usize, usize, Scalar)>,
        C: Vec<(usize, usize, Scalar)>,
    ) -> Result<Self, R1CSError> {
        // Pad to power of 2 if needed
        let num_vars_padded = std::cmp::max(num_vars, num_inputs + 1).next_power_of_two();
        let num_cons_padded = std::cmp::max(num_cons, 2).next_power_of_two();

        let shape = R1CSShape::new(
            num_cons_padded,
            num_vars_padded,
            num_inputs,
            &A,
            &B,
            &C,
        );

        let digest = shape.get_digest();

        Ok(R1CSInstance { shape, digest })
    }

    pub fn is_sat(&self, vars: &[Scalar], inputs: &[Scalar]) -> Result<bool, R1CSError> {
        // Pad vars if needed
        let num_padded_vars = self.shape.get_num_vars();
        let padded_vars = if vars.len() < num_padded_vars {
            let mut v = vars.to_vec();
            v.resize(num_padded_vars, Scalar::zero());
            v
        } else {
            vars.to_vec()
        };

        Ok(self.shape.is_sat(&padded_vars, inputs))
    }
}

/// Generators for R1CS commitments in SNARK mode
/// Uses full sparse matrix polynomial commitment generators
pub struct R1CSCommitmentGens {
    pub gens: SparseMatPolyCommitmentGens,
}

impl R1CSCommitmentGens {
    /// Create new generators for R1CS commitment
    /// 
    /// IMPORTANT: num_nz_entries should be the max of nnz_a, nnz_b, nnz_c
    /// Using a smaller value will cause panics during commitment!
    pub fn new(label: &'static [u8], num_cons: usize, num_vars: usize, num_nz_entries: usize) -> Self {
        let num_poly_vars_x = num_cons.log_2();
        let num_poly_vars_y = (2 * num_vars).log_2();
        
        R1CSCommitmentGens {
            gens: SparseMatPolyCommitmentGens::new(
                label,
                num_poly_vars_x,
                num_poly_vars_y,
                num_nz_entries.next_power_of_two(),
                3, // 3 matrices: A, B, C
            ),
        }
    }
}

/// Commitment to R1CS matrices (A, B, C)
/// Uses full sparse matrix polynomial commitment
#[derive(Debug, Clone, Serialize, Deserialize, CanonicalSerialize, CanonicalDeserialize)]
pub struct R1CSCommitment {
    pub num_cons: usize,
    pub num_vars: usize,
    pub num_inputs: usize,
    pub comm: SparseMatPolyCommitment,
}

impl AppendToTranscript for R1CSCommitment {
    fn append_to_transcript(&self, _label: &'static [u8], transcript: &mut Transcript) {
        // Match arkworks-spartan's R1CSCommitment.append_to_transcript
        transcript.append_u64(b"num_cons", self.num_cons as u64);
        transcript.append_u64(b"num_vars", self.num_vars as u64);
        transcript.append_u64(b"num_inputs", self.num_inputs as u64);
        self.comm.append_to_transcript(b"comm", transcript);
    }
}

/// Decommitment for R1CS (used by prover)
/// Contains the dense representation needed for the full evaluation proof
pub struct R1CSDecommitment {
    /// Dense representation of sparse matrices for full evaluation proof
    pub dense: MultiSparseMatPolynomialAsDense,
}

impl R1CSShape {
    /// Commit to the R1CS shape for SNARK mode
    /// Uses full sparse matrix polynomial commitment for cryptographic soundness
    pub fn commit(&self, gens: &R1CSCommitmentGens) -> (R1CSCommitment, R1CSDecommitment) {
        let timer = Timer::new("R1CSShape::commit");
        
        // Convert to sparse_mlpoly_full format
        let mat_a = self.to_sparse_mat_poly_full(&self.A);
        let mat_b = self.to_sparse_mat_poly_full(&self.B);
        let mat_c = self.to_sparse_mat_poly_full(&self.C);
        
        // Create commitment using full sparse matrix polynomial commitment
        let (comm, dense) = SparseMatPolynomialFull::multi_commit(
            &[&mat_a, &mat_b, &mat_c],
            &gens.gens,
        );
        
        timer.stop();
        
        (
            R1CSCommitment {
                num_cons: self.num_cons,
                num_vars: self.num_vars,
                num_inputs: self.num_inputs,
                comm,
            },
            R1CSDecommitment { dense },
        )
    }
    
    /// Convert sparse matrix polynomial to the full format
    fn to_sparse_mat_poly_full(&self, mat: &SparseMatPolynomial) -> SparseMatPolynomialFull {
        use crate::sparse_mlpoly_full::SparseMatEntry as SparseMatEntryFull;
        
        let entries: Vec<SparseMatEntryFull> = mat
            .get_entries()
            .iter()
            .map(|e| SparseMatEntryFull::new(e.row(), e.col(), e.val()))
            .collect();
        
        SparseMatPolynomialFull::new(mat.num_vars_x(), mat.num_vars_y(), entries)
    }
}

/// Proof of evaluation of R1CS matrices at random point
/// Uses full sparse matrix polynomial evaluation proof (same as original Spartan)
#[derive(Debug, Clone, Serialize, Deserialize, CanonicalSerialize, CanonicalDeserialize)]
pub struct R1CSEvalProof {
    /// Full sparse matrix polynomial evaluation proof
    proof: SparseMatPolyEvalProof,
}

impl R1CSEvalProof {
    /// Prove evaluations of A, B, C at (rx, ry)
    /// Uses full sparse matrix polynomial evaluation proof with memory-checking
    pub fn prove(
        decomm: &R1CSDecommitment,
        rx: &[Scalar],
        ry: &[Scalar],
        evals: &(Scalar, Scalar, Scalar),
        gens: &R1CSCommitmentGens,
        transcript: &mut Transcript,
        random_tape: &mut RandomTape,
    ) -> Self {
        let timer = Timer::new("R1CSEvalProof::prove");
        
        let (eval_A, eval_B, eval_C) = evals;
        let evals_vec = vec![*eval_A, *eval_B, *eval_C];
        
        let proof = SparseMatPolyEvalProof::prove(
            &decomm.dense,
            rx,
            ry,
            &evals_vec,
            &gens.gens,
            transcript,
            random_tape,
        );
        
        timer.stop();
        
        R1CSEvalProof { proof }
    }
    
    /// Verify evaluations of A, B, C
    pub fn verify(
        &self,
        comm: &R1CSCommitment,
        rx: &[Scalar],
        ry: &[Scalar],
        evals: &(Scalar, Scalar, Scalar),
        gens: &R1CSCommitmentGens,
        transcript: &mut Transcript,
    ) -> Result<(), ProofVerifyError> {
        let timer = Timer::new("R1CSEvalProof::verify");
        
        let (eval_A, eval_B, eval_C) = evals;
        let evals_vec = vec![*eval_A, *eval_B, *eval_C];
        
        let result = self.proof.verify(
            &comm.comm,
            rx,
            ry,
            &evals_vec,
            &gens.gens,
            transcript,
        );
        
        timer.stop();
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_r1cs_shape_basic() {
        // Simple constraint: x * x = x (satisfied by x=0 or x=1)
        let A = vec![(0, 0, Scalar::one())];
        let B = vec![(0, 0, Scalar::one())];
        let C = vec![(0, 0, Scalar::one())];

        let shape = R1CSShape::new(1, 2, 0, &A, &B, &C);
        
        // x = 1 should satisfy
        let vars = vec![Scalar::one(), Scalar::zero()];
        let inputs: Vec<Scalar> = vec![];
        assert!(shape.is_sat(&vars, &inputs));

        // x = 0 should satisfy
        let vars = vec![Scalar::zero(), Scalar::zero()];
        assert!(shape.is_sat(&vars, &inputs));
    }
}

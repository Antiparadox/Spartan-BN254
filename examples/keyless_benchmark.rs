//! Unified Keyless Circuit Benchmark with Detailed Prover Breakdown
//!
//! Supports both Hyrax (default) and KZG modes.
//!
//! Usage:
//!   # Hyrax mode (no trusted setup)
//!   cargo run --release --example keyless_benchmark
//!
//!   # KZG mode (smaller proofs, requires SRS)
//!   cargo run --release --features kzg --example keyless_benchmark

use std::fs::File;
use std::io::Read;
use std::time::Instant;

use ark_serialize::CanonicalSerialize;
use spartan_bn254::{
    InputsAssignment, Instance, R1CSShape, Scalar, VarsAssignment, Math,
    r1cs_reader::R1CS,
    r1cs::R1CSCommitmentGens,
    R1CSGens,
    r1csproof::R1CSProof,
    hyrax::EqPolynomial,
    sparse_mlpoly_full::{SparseMatPolyEvalProof, PolyEvalNetwork, PolyEvalNetworkProof},
    random::RandomTape,
    transcript::{ProofTranscript, AppendToTranscript},
};

#[cfg(feature = "kzg")]
use spartan_bn254::kzg::KZGSrs;

const R1CS_PATH: &str = "/Users/yinuo/Desktop/Universal ZK/keyless-zk-proofs/main.r1cs";
const WTNS_PATH: &str = "/Users/yinuo/Desktop/Universal ZK/keyless-zk-proofs/witness_real.wtns";

#[cfg(feature = "kzg")]
const SRS_PATH: &str = "data/kzg_srs_keyless.bin";  // Local, ignored by git

fn parse_wtns(path: &str) -> Result<Vec<Scalar>, Box<dyn std::error::Error>> {
    let mut file = File::open(path)?;
    let mut data = Vec::new();
    file.read_to_end(&mut data)?;
    
    if data.len() < 4 || &data[0..4] != b"wtns" {
        return Err("Invalid wtns magic".into());
    }
    
    let num_sections = u32::from_le_bytes(data[8..12].try_into()?);
    let mut offset = 12;
    let mut witness_values = Vec::new();
    
    for _ in 0..num_sections {
        if offset + 12 > data.len() { break; }
        let section_id = u32::from_le_bytes(data[offset..offset+4].try_into()?);
        let section_size = u64::from_le_bytes(data[offset+4..offset+12].try_into()?) as usize;
        offset += 12;
        
        if section_id == 2 {
            let num_witnesses = section_size / 32;
            for i in 0..num_witnesses {
                let start = offset + i * 32;
                if start + 32 > data.len() { break; }
                let mut bytes = [0u8; 32];
                bytes.copy_from_slice(&data[start..start+32]);
                witness_values.push(Scalar::from_bytes(&bytes).unwrap_or_else(|| {
                    Scalar::from_u64(u64::from_le_bytes(bytes[0..8].try_into().unwrap()))
                }));
            }
        }
        offset += section_size;
    }
    Ok(witness_values)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    #[cfg(feature = "kzg")]
    let mode = "KZG";
    #[cfg(not(feature = "kzg"))]
    let mode = "HYRAX";

    println!("╔══════════════════════════════════════════════════════════════════════════╗");
    println!("║   SPARTAN-BN254 KEYLESS BENCHMARK ({} MODE)                       ║", mode);
    println!("╚══════════════════════════════════════════════════════════════════════════╝\n");

    // =========================================================================
    // PHASE 1: Load R1CS
    // =========================================================================
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    println!("  LOADING CIRCUIT");
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    
    let load_start = Instant::now();
    let r1cs = R1CS::from_file(R1CS_PATH)?;
    let load_time = load_start.elapsed();
    
    let stats = r1cs.stats();
    let num_cons = stats.num_constraints;
    let num_inputs = stats.num_pub_inputs;
    let num_prv_vars = r1cs.num_private_vars();
    let num_cons_padded = num_cons.next_power_of_two();
    let num_vars_padded = std::cmp::max(num_prv_vars, num_inputs + 1).next_power_of_two();
    let max_nnz = stats.nnz_a.max(stats.nnz_b).max(stats.nnz_c);
    let total_nnz = stats.nnz_a + stats.nnz_b + stats.nnz_c;

    println!("  Constraints:      {:>12} (padded: 2^{})", num_cons, num_cons_padded.log_2());
    println!("  Variables:        {:>12} (padded: 2^{})", stats.num_variables, num_vars_padded.log_2());
    println!("  Public inputs:    {:>12}", num_inputs);
    println!("  NNZ (A/B/C):      {}/{}/{}", stats.nnz_a, stats.nnz_b, stats.nnz_c);
    println!("  Total NNZ:        {:>12}", total_nnz);
    println!("  Load time:        {:>12.2?}", load_time);

    // =========================================================================
    // PHASE 2: Convert to Spartan Format
    // =========================================================================
    println!("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    println!("  PREPROCESSING");
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    
    let convert_start = Instant::now();
    let (a, b, c) = r1cs.to_sparse_matrices_padded(num_vars_padded);
    let shape = R1CSShape::new(num_cons_padded, num_vars_padded, num_inputs, &a, &b, &c);
    let inst = Instance::from_shape(shape);
    let convert_time = convert_start.elapsed();
    println!("  Convert time:     {:>12.2?}", convert_time);

    // Load witness
    let witness_start = Instant::now();
    let witness = parse_wtns(WTNS_PATH)?;
    let inputs_vec: Vec<Scalar> = witness[1..1+num_inputs].to_vec();
    let mut vars: Vec<Scalar> = witness[1+num_inputs..].to_vec();
    vars.resize(num_vars_padded, Scalar::zero());
    let witness_time = witness_start.elapsed();
    println!("  Witness load:     {:>12.2?}", witness_time);

    // Create generators
    let gens_start = Instant::now();
    let gens_r1cs_sat = R1CSGens::new(b"gens_r1cs_sat", num_cons_padded, num_vars_padded);
    
    #[cfg(not(feature = "kzg"))]
    let gens_r1cs_eval = R1CSCommitmentGens::new(b"gens_r1cs_eval", num_cons_padded, num_vars_padded, max_nnz);
    
    #[cfg(feature = "kzg")]
    let gens_r1cs_eval = {
        println!("  Loading KZG SRS...");
        let srs = KZGSrs::load_or_generate(SRS_PATH, max_nnz.next_power_of_two() * 8, 0xDEADBEEF)?;
        R1CSCommitmentGens::new_with_kzg(b"gens_r1cs_eval", num_cons_padded, num_vars_padded, max_nnz, srs)
    };
    
    let gens_time = gens_start.elapsed();
    println!("  Generators:       {:>12.2?}", gens_time);

    // Encode
    let encode_start = Instant::now();
    let (comm, decomm) = inst.inst.commit(&gens_r1cs_eval);
    let encode_time = encode_start.elapsed();
    println!("  Encode time:      {:>12.2?}", encode_time);

    // =========================================================================
    // PHASE 3: Proving with Detailed Breakdown
    // =========================================================================
    println!("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    println!("  PROVING (Detailed Breakdown)");
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");

    let total_prove_start = Instant::now();
    let mut random_tape = RandomTape::new(b"snark_proof");
    let mut transcript = merlin::Transcript::new(b"keyless_snark");
    
    transcript.append_protocol_name(b"Spartan SNARK proof");
    comm.append_to_transcript(b"comm", &mut transcript);

    // Part 1: R1CS Satisfiability Proof
    let r1cs_sat_start = Instant::now();
    let padded_vars = {
        let num_padded = inst.inst.get_num_vars();
        let mut v = vars.clone();
        v.resize(num_padded, Scalar::zero());
        v
    };
    let (r1cs_sat_proof, rx, ry) = R1CSProof::prove(
        &inst.inst, padded_vars, &inputs_vec, &gens_r1cs_sat, &mut transcript, &mut random_tape,
    );
    let r1cs_sat_time = r1cs_sat_start.elapsed();

    // Part 2: Instance Evaluations
    let inst_evals_start = Instant::now();
    let inst_evals = inst.inst.evaluate(&rx, &ry);
    let inst_evals_time = inst_evals_start.elapsed();

    // Part 3: R1CS Eval Proof (Lookup Argument) - Detailed Breakdown
    let r1cs_eval_start = Instant::now();
    transcript.append_protocol_name(b"Sparse polynomial evaluation proof");
    let evals_vec = vec![inst_evals.0, inst_evals.1, inst_evals.2];
    let dense = &decomm.dense;

    // 3a: EqPolynomial evaluation
    let eq_start = Instant::now();
    let (mem_rx, mem_ry) = {
        let (rx_ext, ry_ext) = SparseMatPolyEvalProof::equalize(&rx, &ry);
        (EqPolynomial::new(rx_ext).evals(), EqPolynomial::new(ry_ext).evals())
    };
    let eq_time = eq_start.elapsed();

    // 3b: Derefs computation
    let derefs_compute_start = Instant::now();
    let derefs = dense.deref(&mem_rx, &mem_ry);
    let derefs_compute_time = derefs_compute_start.elapsed();

    // 3c: Derefs commitment
    let derefs_commit_start = Instant::now();
    #[cfg(not(feature = "kzg"))]
    let comm_derefs = {
        let c = derefs.commit(&gens_r1cs_eval.gens.gens_derefs);
        c.append_to_transcript(b"comm_poly_row_col_ops_val", &mut transcript);
        c
    };
    #[cfg(feature = "kzg")]
    let comm_derefs = {
        let c = derefs.commit_kzg(&gens_r1cs_eval.gens.gens_derefs_kzg);
        c.append_to_transcript(b"comm_poly_row_col_ops_val", &mut transcript);
        c
    };
    let derefs_commit_time = derefs_commit_start.elapsed();

    // 3d: Network construction
    let network_start = Instant::now();
    let r_mem_check = transcript.challenge_vector(b"challenge_r_hash", 2);
    let mut net = PolyEvalNetwork::new(dense, &derefs, &mem_rx, &mem_ry, &(r_mem_check[0], r_mem_check[1]));
    let network_time = network_start.elapsed();

    // 3e: Network proof
    let network_proof_start = Instant::now();
    let poly_eval_network_proof = PolyEvalNetworkProof::prove(
        &mut net, dense, &derefs, &evals_vec, &gens_r1cs_eval.gens, &mut transcript, &mut random_tape,
    );
    let network_proof_time = network_proof_start.elapsed();

    let r1cs_eval_time = r1cs_eval_start.elapsed();
    let total_prove_time = total_prove_start.elapsed();

    // Construct proof for size measurement
    let proof = SparseMatPolyEvalProof::from_parts(comm_derefs.clone(), poly_eval_network_proof);
    let mut r1cs_eval_proof_bytes = Vec::new();
    proof.serialize_compressed(&mut r1cs_eval_proof_bytes)?;

    let mut r1cs_sat_proof_bytes = Vec::new();
    r1cs_sat_proof.serialize_compressed(&mut r1cs_sat_proof_bytes)?;

    // =========================================================================
    // PHASE 4: Generate Full SNARK Proof and Verify
    // =========================================================================
    println!("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    println!("  FULL SNARK (for verification test)");
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");

    use spartan_bn254::{SNARK, SNARKGens};
    let vars_assignment = VarsAssignment::from_scalars(vars.clone());
    let inputs_assignment = InputsAssignment::from_scalars(inputs_vec.clone());
    
    // Create generators (reusing the same SRS for KZG)
    #[cfg(not(feature = "kzg"))]
    let full_gens = SNARKGens::new(num_cons_padded, num_vars_padded, num_inputs, max_nnz);
    #[cfg(feature = "kzg")]
    let full_gens = {
        let srs = KZGSrs::load_or_generate(SRS_PATH, max_nnz.next_power_of_two() * 8, 0xDEADBEEF)?;
        SNARKGens::new_with_kzg(num_cons_padded, num_vars_padded, num_inputs, max_nnz, srs)
    };
    
    let (full_comm, _) = SNARK::encode(&inst, &full_gens);
    let mut prover_transcript = merlin::Transcript::new(b"keyless_snark");
    let full_proof = SNARK::prove(&inst, &full_comm, &decomm, vars_assignment, &inputs_assignment, &full_gens, &mut prover_transcript);
    
    let mut full_proof_bytes = Vec::new();
    full_proof.serialize_compressed(&mut full_proof_bytes)?;

    // Now time ONLY the verification
    println!("  Timing verification only...");
    let verify_start = Instant::now();
    let mut verifier_transcript = merlin::Transcript::new(b"keyless_snark");
    let verify_result = full_proof.verify(&full_comm, &inputs_assignment, &mut verifier_transcript, &full_gens);
    let verify_time = verify_start.elapsed();

    match verify_result {
        Ok(()) => println!("  ✅ Verification PASSED"),
        Err(e) => println!("  ❌ Verification FAILED: {:?}", e),
    }
    println!("  Verify time:      {:>12.2?}", verify_time);

    // =========================================================================
    // RESULTS
    // =========================================================================
    println!("\n");
    println!("╔══════════════════════════════════════════════════════════════════════════╗");
    println!("║                      PROVER TIME BREAKDOWN ({})                     ║", mode);
    println!("╠══════════════════════════════════════════════════════════════════════════╣");
    println!("║                                                                          ║");
    println!("║  R1CS SATISFIABILITY PROOF                                               ║");
    println!("║    Sumcheck + Witness Opening:          {:>10.2?}  ({:>5.1}%)      ║", 
             r1cs_sat_time, r1cs_sat_time.as_secs_f64() / total_prove_time.as_secs_f64() * 100.0);
    println!("║                                                                          ║");
    println!("║  INSTANCE EVALUATIONS                                                    ║");
    println!("║    Compute A(rx,ry), B(rx,ry), C(rx,ry): {:>10.2?}  ({:>5.1}%)      ║",
             inst_evals_time, inst_evals_time.as_secs_f64() / total_prove_time.as_secs_f64() * 100.0);
    println!("║                                                                          ║");
    println!("║  R1CS EVAL PROOF (Lookup Argument)                                       ║");
    println!("║    [a] EqPolynomial evaluation:         {:>10.2?}  ({:>5.1}%)      ║",
             eq_time, eq_time.as_secs_f64() / total_prove_time.as_secs_f64() * 100.0);
    println!("║    [b] Derefs computation:              {:>10.2?}  ({:>5.1}%)      ║",
             derefs_compute_time, derefs_compute_time.as_secs_f64() / total_prove_time.as_secs_f64() * 100.0);
    #[cfg(not(feature = "kzg"))]
    println!("║    [c] Derefs commitment (Hyrax MSM):   {:>10.2?}  ({:>5.1}%) ◀━━  ║",
             derefs_commit_time, derefs_commit_time.as_secs_f64() / total_prove_time.as_secs_f64() * 100.0);
    #[cfg(feature = "kzg")]
    println!("║    [c] Derefs commitment (KZG):         {:>10.2?}  ({:>5.1}%)      ║",
             derefs_commit_time, derefs_commit_time.as_secs_f64() / total_prove_time.as_secs_f64() * 100.0);
    println!("║    [d] Network construction:            {:>10.2?}  ({:>5.1}%)      ║",
             network_time, network_time.as_secs_f64() / total_prove_time.as_secs_f64() * 100.0);
    #[cfg(not(feature = "kzg"))]
    println!("║    [e] Network proof (Hyrax openings):  {:>10.2?}  ({:>5.1}%)      ║",
             network_proof_time, network_proof_time.as_secs_f64() / total_prove_time.as_secs_f64() * 100.0);
    #[cfg(feature = "kzg")]
    println!("║    [e] Network proof (KZG openings):    {:>10.2?}  ({:>5.1}%) ◀━━  ║",
             network_proof_time, network_proof_time.as_secs_f64() / total_prove_time.as_secs_f64() * 100.0);
    println!("║    ─────────────────────────────────────────────────────────             ║");
    println!("║    Subtotal R1CS Eval:                  {:>10.2?}  ({:>5.1}%)      ║",
             r1cs_eval_time, r1cs_eval_time.as_secs_f64() / total_prove_time.as_secs_f64() * 100.0);
    println!("║                                                                          ║");
    println!("╠══════════════════════════════════════════════════════════════════════════╣");
    println!("║  TOTAL PROVE TIME:                      {:>10.2?}  (100.0%)      ║", total_prove_time);
    println!("╚══════════════════════════════════════════════════════════════════════════╝");

    println!("\n");
    println!("╔══════════════════════════════════════════════════════════════════════════╗");
    println!("║                      PROOF SIZE BREAKDOWN ({})                      ║", mode);
    println!("╠══════════════════════════════════════════════════════════════════════════╣");
    println!("║  R1CS Sat Proof:                      {:>10} bytes ({:>6.1} KB)   ║", 
             r1cs_sat_proof_bytes.len(), r1cs_sat_proof_bytes.len() as f64 / 1024.0);
    println!("║  R1CS Eval Proof:                     {:>10} bytes ({:>6.1} KB)   ║", 
             r1cs_eval_proof_bytes.len(), r1cs_eval_proof_bytes.len() as f64 / 1024.0);
    println!("║  inst_evals (3 scalars):              {:>10} bytes               ║", 96);
    println!("║  ──────────────────────────────────────────────────────────────────      ║");
    println!("║  TOTAL PROOF:                         {:>10} bytes ({:>6.1} KB)   ║", 
             full_proof_bytes.len(), full_proof_bytes.len() as f64 / 1024.0);
    println!("╚══════════════════════════════════════════════════════════════════════════╝");

    println!("\n");
    println!("╔══════════════════════════════════════════════════════════════════════════╗");
    println!("║                           SUMMARY ({})                             ║", mode);
    println!("╠══════════════════════════════════════════════════════════════════════════╣");
    println!("║  Circuit:          Keyless (Aptos)                                       ║");
    println!("║  Constraints:      {:>12}                                         ║", num_cons);
    println!("║  Total NNZ:        {:>12}                                         ║", total_nnz);
    println!("║  ──────────────────────────────────────────────────────────────────      ║");
    println!("║  Encode Time:      {:>12.1}s                                         ║", encode_time.as_secs_f64());
    println!("║  Prove Time:       {:>12.1}s                                         ║", total_prove_time.as_secs_f64());
    println!("║  Verify Time:      {:>12.1}ms                                        ║", verify_time.as_secs_f64() * 1000.0);
    println!("║  Proof Size:       {:>12.1} KB                                       ║", full_proof_bytes.len() as f64 / 1024.0);
    #[cfg(not(feature = "kzg"))]
    println!("║  Trusted Setup:    No                                                    ║");
    #[cfg(feature = "kzg")]
    println!("║  Trusted Setup:    Yes (KZG SRS)                                         ║");
    println!("╚══════════════════════════════════════════════════════════════════════════╝");

    Ok(())
}

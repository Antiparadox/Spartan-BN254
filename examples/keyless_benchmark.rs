//! Spartan-BN254 SNARK Benchmark on Keyless Circuit
//!
//! Usage:
//!   cargo run --example keyless_benchmark --release
//!   cargo run --example keyless_benchmark --release -- --dummy  # Use dummy witness

use ark_serialize::CanonicalSerialize;
use spartan_bn254::{
    InputsAssignment, Instance, R1CSShape, Scalar, SNARK, SNARKGens, VarsAssignment,
    r1cs_reader::R1CS,
};
use std::fs::File;
use std::io::Read;
use std::time::Instant;

const R1CS_PATH: &str = "/Users/yinuo/Desktop/Universal ZK/keyless-zk-proofs/main.r1cs";
const WTNS_PATH: &str = "/Users/yinuo/Desktop/Universal ZK/keyless-zk-proofs/witness_real.wtns";

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
    let use_dummy = std::env::args().any(|a| a == "--dummy");
    
    println!("╔══════════════════════════════════════════════════════════════════╗");
    println!("║          Spartan-BN254 SNARK - Keyless Circuit Benchmark         ║");
    println!("╚══════════════════════════════════════════════════════════════════╝\n");

    // Load R1CS
    println!("📂 Loading R1CS...");
    let start = Instant::now();
    let r1cs = R1CS::from_file(R1CS_PATH)?;
    let r1cs_time = start.elapsed();
    
    let stats = r1cs.stats();
    println!("  Constraints:    {:>12}", stats.num_constraints);
    println!("  Variables:      {:>12}", stats.num_variables);
    println!("  Public inputs:  {:>12}", stats.num_pub_inputs);
    println!("  NNZ (A/B/C):    {}/{}/{}", stats.nnz_a, stats.nnz_b, stats.nnz_c);
    println!("  Load time:      {:?}", r1cs_time);

    // Dimensions
    let num_cons = stats.num_constraints;
    let num_inputs = stats.num_pub_inputs;
    let num_prv_vars = r1cs.num_private_vars();
    let num_cons_padded = num_cons.next_power_of_two();
    let num_vars_padded = std::cmp::max(num_prv_vars, num_inputs + 1).next_power_of_two();
    let max_nnz = stats.nnz_a.max(stats.nnz_b).max(stats.nnz_c);

    println!("\n  Padded: {} cons (2^{}), {} vars (2^{})", 
        num_cons_padded, num_cons_padded.trailing_zeros(),
        num_vars_padded, num_vars_padded.trailing_zeros());

    // Convert R1CS
    println!("\n⚙️  Converting to Spartan format...");
    let start = Instant::now();
    let (A, B, C) = r1cs.to_sparse_matrices_padded(num_vars_padded);
    let shape = R1CSShape::new(num_cons_padded, num_vars_padded, num_inputs, &A, &B, &C);
    let inst = Instance::from_shape(shape);
    let convert_time = start.elapsed();
    println!("  Done in {:?}", convert_time);

    // Witness
    let (vars, inputs_vec) = if use_dummy {
        println!("\n⚙️  Using DUMMY witness...");
        (vec![Scalar::zero(); num_vars_padded], vec![Scalar::zero(); num_inputs])
    } else {
        println!("\n📂 Loading real witness...");
        let witness = parse_wtns(WTNS_PATH)?;
        println!("  Loaded {} values", witness.len());
        
        let inputs_vec: Vec<Scalar> = witness[1..1+num_inputs].to_vec();
        let mut vars: Vec<Scalar> = witness[1+num_inputs..].to_vec();
        vars.resize(num_vars_padded, Scalar::zero());
        (vars, inputs_vec)
    };
    
    let vars_assignment = VarsAssignment::from_scalars(vars);
    let inputs_assignment = InputsAssignment::from_scalars(inputs_vec);

    // Check satisfaction
    println!("\n🔍 Checking R1CS satisfaction...");
    match inst.is_sat(&vars_assignment, &inputs_assignment) {
        Ok(true) => println!("  ✅ Witness SATISFIES R1CS"),
        Ok(false) => println!("  ⚠️  Witness does NOT satisfy R1CS"),
        Err(e) => println!("  ❌ Error: {:?}", e),
    }

    // Generators
    println!("\n⚙️  Creating generators (max_nnz={})...", max_nnz);
    let start = Instant::now();
    let gens = SNARKGens::new(num_cons_padded, num_vars_padded, num_inputs, max_nnz);
    let gens_time = start.elapsed();
    println!("  Done in {:?}", gens_time);

    // Encode
    println!("\n📝 Encoding (preprocessing)...");
    let start = Instant::now();
    let (comm, decomm) = SNARK::encode(&inst, &gens);
    let encode_time = start.elapsed();
    println!("  Encode time: {:?}", encode_time);

    // Prove
    println!("\n🔐 Proving...");
    let start = Instant::now();
    let mut prover_transcript = merlin::Transcript::new(b"keyless_snark");
    let proof = SNARK::prove(&inst, &comm, &decomm, vars_assignment, &inputs_assignment, &gens, &mut prover_transcript);
    let prove_time = start.elapsed();
    
    // Use CanonicalSerialize for consistent proof size with arkworks-spartan
    let mut proof_bytes = Vec::new();
    proof.serialize_compressed(&mut proof_bytes)?;
    println!("  Prove time:  {:?}", prove_time);
    println!("  Proof size:  {} bytes ({:.2} KB) [CanonicalSerialize]", proof_bytes.len(), proof_bytes.len() as f64 / 1024.0);

    // Verify
    println!("\n✔️  Verifying...");
    let start = Instant::now();
    let mut verifier_transcript = merlin::Transcript::new(b"keyless_snark");
    let result = proof.verify(&comm, &inputs_assignment, &mut verifier_transcript, &gens);
    let verify_time = start.elapsed();
    
    match result {
        Ok(()) => println!("  ✅ Verification PASSED ({:?})", verify_time),
        Err(e) => println!("  ❌ Verification FAILED: {:?} ({:?})", e, verify_time),
    }

    // Summary
    println!("\n╔══════════════════════════════════════════════════════════════════╗");
    println!("║                           SUMMARY                                ║");
    println!("╠══════════════════════════════════════════════════════════════════╣");
    println!("║  Constraints:     {:>10}                                     ║", num_cons);
    println!("║  Non-Zero Entries:{:>10}                                     ║", stats.nnz_a + stats.nnz_b + stats.nnz_c);
    println!("╠══════════════════════════════════════════════════════════════════╣");
    println!("║  Encode:          {:>10.1}s                                    ║", encode_time.as_secs_f64());
    println!("║  Prove:           {:>10.1}s                                    ║", prove_time.as_secs_f64());
    println!("║  Verify:          {:>10.1}ms                                   ║", verify_time.as_secs_f64() * 1000.0);
    println!("║  Proof Size:      {:>10.1} KB                                  ║", proof_bytes.len() as f64 / 1024.0);
    println!("╚══════════════════════════════════════════════════════════════════╝");

    Ok(())
}

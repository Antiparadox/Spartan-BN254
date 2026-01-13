//! Benchmark Spartan-BN254 SNARK on the Keyless circuit

use spartan_bn254::{
    InputsAssignment, Instance, R1CS, R1CSShape, Scalar, SNARK, SNARKGens, VarsAssignment,
};
use std::time::Instant;

fn main() {
    println!("╔════════════════════════════════════════════════════════════╗");
    println!("║       Spartan-BN254 Keyless Circuit SNARK Benchmark        ║");
    println!("╚════════════════════════════════════════════════════════════╝\n");

    let r1cs_path = "../keyless-zk-proofs/main.r1cs";

    println!("📋 Loading R1CS from: {}", r1cs_path);
    println!("────────────────────────────────────────────────────────────\n");

    let start = Instant::now();
    let r1cs = match R1CS::from_file(r1cs_path) {
        Ok(r1cs) => r1cs,
        Err(e) => {
            eprintln!("❌ Failed to load R1CS: {}", e);
            return;
        }
    };
    let load_time = start.elapsed();

    let stats = r1cs.stats();
    println!("{}", stats);
    println!("  Load time: {:?}", load_time);

    // Calculate padded dimensions for Spartan
    let num_cons = stats.num_constraints;
    let num_vars = stats.num_variables;
    let num_inputs = stats.num_pub_inputs;
    let num_cons_padded = num_cons.next_power_of_two();
    let num_vars_padded = std::cmp::max(num_vars, num_inputs + 1).next_power_of_two();

    println!("\n📊 Spartan Requirements:");
    println!("────────────────────────────────────────────────────────────");
    println!(
        "  Padded constraints: {} (2^{})",
        num_cons_padded,
        (num_cons_padded as f64).log2() as u32
    );
    println!(
        "  Padded variables:   {} (2^{})",
        num_vars_padded,
        (num_vars_padded as f64).log2() as u32
    );
    println!(
        "  Padding overhead:   {:.1}%",
        ((num_cons_padded as f64 / num_cons as f64) - 1.0) * 100.0
    );

    // Convert R1CS to Spartan format
    println!("\n⚙️  Converting R1CS to Spartan format...");
    let start = Instant::now();
    let (A, B, C) = r1cs.to_sparse_matrices();
    let convert_time = start.elapsed();
    println!("  Conversion time: {:?}", convert_time);

    // Create R1CS shape
    println!("\n⚙️  Creating R1CS shape...");
    let start = Instant::now();
    let shape = R1CSShape::new(num_cons_padded, num_vars_padded, num_inputs, &A, &B, &C);
    let shape_time = start.elapsed();
    println!("  Shape creation time: {:?}", shape_time);

    let inst = Instance::from_shape(shape);

    // Create dummy witness (all zeros for now - actual witness would come from circuit execution)
    println!("\n⚙️  Creating witness...");
    let mut vars = vec![Scalar::zero(); num_vars_padded];
    // Set first variable to 1 for a simple satisfying assignment (for testing)
    if num_vars_padded > 0 {
        vars[0] = Scalar::one();
    }
    let vars_assignment = VarsAssignment::from_scalars(vars);
    let inputs = InputsAssignment::from_scalars(vec![Scalar::zero(); num_inputs]);

    // Check satisfaction (expected to fail with dummy witness)
    println!("\n🔍 Checking R1CS satisfaction with dummy witness...");
    match inst.is_sat(&vars_assignment, &inputs) {
        Ok(true) => println!("  ✅ Dummy witness satisfies R1CS (unexpected for real circuit)"),
        Ok(false) => println!("  ⚠️  Dummy witness does NOT satisfy R1CS (expected)"),
        Err(e) => println!("  ❌ Error checking satisfaction: {}", e),
    }

    // Create generators
    println!("\n⚙️  Creating SNARK generators...");
    let start = Instant::now();
    let gens = SNARKGens::new(num_cons_padded, num_vars_padded, num_inputs);
    let gens_time = start.elapsed();
    println!("  Generator creation time: {:?}", gens_time);

    // Encode (preprocess)
    println!("\n⚙️  Encoding (preprocessing)...");
    let start = Instant::now();
    let (comm, decomm) = SNARK::encode(&inst, &gens);
    let encode_time = start.elapsed();
    println!("  Encode time: {:?}", encode_time);

    // Prove (will fail verification if witness doesn't satisfy, but timing is still useful)
    println!("\n🔐 Proving...");
    let start = Instant::now();
    let mut prover_transcript = merlin::Transcript::new(b"keyless_snark");
    let proof = SNARK::prove(
        &inst,
        &comm,
        &decomm,
        vars_assignment.clone(),
        &inputs,
        &gens,
        &mut prover_transcript,
    );
    let prove_time = start.elapsed();
    println!("  Prove time: {:?}", prove_time);

    // Serialize proof to get size
    let proof_bytes = bincode::serialize(&proof).unwrap();
    println!("  Proof size: {} bytes ({:.2} KB)", proof_bytes.len(), proof_bytes.len() as f64 / 1024.0);

    // Verify
    println!("\n✔️  Verifying...");
    let start = Instant::now();
    let mut verifier_transcript = merlin::Transcript::new(b"keyless_snark");
    let verify_result = proof.verify(&comm, &inputs, &mut verifier_transcript, &gens);
    let verify_time = start.elapsed();
    
    match verify_result {
        Ok(()) => println!("  ✅ Verification PASSED"),
        Err(e) => println!("  ❌ Verification FAILED: {:?} (expected with dummy witness)", e),
    }
    println!("  Verify time: {:?}", verify_time);

    // Summary
    println!("\n╔════════════════════════════════════════════════════════════╗");
    println!("║                        SUMMARY                             ║");
    println!("╚════════════════════════════════════════════════════════════╝");
    println!("  R1CS Load:       {:?}", load_time);
    println!("  Format Convert:  {:?}", convert_time);
    println!("  Shape Creation:  {:?}", shape_time);
    println!("  Generator Setup: {:?}", gens_time);
    println!("  Encode (Preproc):{:?}", encode_time);
    println!("  Prove Time:      {:?}", prove_time);
    println!("  Verify Time:     {:?}", verify_time);
    println!("  Proof Size:      {} bytes", proof_bytes.len());
    println!("");
    println!("  Total Prover Time: {:?}", load_time + convert_time + shape_time + gens_time + encode_time + prove_time);
    println!("  Constraints: {} (padded: {})", num_cons, num_cons_padded);
    println!("  Variables:   {} (padded: {})", num_vars, num_vars_padded);
}

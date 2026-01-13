//! Analyze proof size breakdown

use spartan_bn254::{
    InputsAssignment, Instance, R1CSShape, Scalar, SNARK, SNARKGens, VarsAssignment,
};

fn main() {
    // Create a simple test circuit
    let num_cons = 1024; // Small for quick analysis
    let num_vars = 1024;
    let num_inputs = 1;

    let one = Scalar::one();
    let A = vec![(0, 0, one)];
    let B = vec![(0, 0, one)];
    let C = vec![(0, 0, one)];

    let shape = R1CSShape::new(num_cons, num_vars, num_inputs, &A, &B, &C);
    let inst = Instance::from_shape(shape);

    let vars = vec![Scalar::one(); num_vars];
    let vars_assignment = VarsAssignment::from_scalars(vars);
    let inputs = InputsAssignment::from_scalars(vec![Scalar::one()]);

    let gens = SNARKGens::new(num_cons, num_vars, num_inputs);
    let (comm, decomm) = SNARK::encode(&inst, &gens);

    let mut transcript = merlin::Transcript::new(b"test");
    let proof = SNARK::prove(&inst, &comm, &decomm, vars_assignment, &inputs, &gens, &mut transcript);

    // Serialize and analyze
    let proof_bytes = bincode::serialize(&proof).unwrap();
    println!("Total proof size: {} bytes", proof_bytes.len());
    
    // The proof contains:
    // 1. r1cs_sat_proof (R1CSProof) - sumcheck proofs + polynomial commitments
    // 2. inst_evals - 3 scalars (eval of A, B, C)
    // 3. r1cs_eval_proof - evaluation proof for R1CS
    // 4. r (rx, ry) - random evaluation points
    
    println!("\nProof components (estimated):");
    println!("  - 3 evaluation scalars (A,B,C): ~96 bytes");
    println!("  - Random points rx, ry: ~{} bytes", (10 + 21) * 32); // log2(cons) + log2(2*vars)
    println!("  - Sumcheck phase 1 (20 rounds): ~{} bytes", 20 * 4 * 32); // degree-3 polynomials
    println!("  - Sumcheck phase 2 (21 rounds): ~{} bytes", 21 * 3 * 32); // degree-2 polynomials
    println!("  - Polynomial commitments: ~{} bytes", 10 * 32); // compressed group elements
    println!("  - ZK proofs (equality, product, etc.): ~{} bytes", 20 * 32);
}

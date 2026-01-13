//! Test the sparse matrix polynomial evaluation proof

use spartan_bn254::{
    dense_mlpoly::EqPolynomial,
    math::Math,
    random::RandomTape,
    scalar::Scalar,
    sparse_mlpoly_full::{
        SparseMatEntry, SparseMatPolynomial, SparseMatPolyCommitmentGens, SparseMatPolyEvalProof,
    },
    transcript::ProofTranscript,
};
use merlin::Transcript;
use std::time::Instant;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("╔══════════════════════════════════════════════════════════════════════╗");
    println!("║     Sparse Matrix Polynomial Evaluation Proof Test                    ║");
    println!("╚══════════════════════════════════════════════════════════════════════╝\n");

    // Create a simple sparse matrix
    // This represents a 4x4 identity matrix for testing
    let num_vars_x = 2; // log2(4) = 2
    let num_vars_y = 2;
    
    let sparse_a = vec![
        SparseMatEntry::new(0, 0, Scalar::one()),
        SparseMatEntry::new(1, 1, Scalar::from_u64(2)),
        SparseMatEntry::new(2, 2, Scalar::from_u64(3)),
        SparseMatEntry::new(3, 3, Scalar::from_u64(4)),
    ];
    let sparse_b = vec![
        SparseMatEntry::new(0, 1, Scalar::one()),
        SparseMatEntry::new(1, 2, Scalar::from_u64(2)),
        SparseMatEntry::new(2, 3, Scalar::from_u64(3)),
        SparseMatEntry::new(3, 0, Scalar::from_u64(4)),
    ];
    let sparse_c = vec![
        SparseMatEntry::new(0, 3, Scalar::from_u64(5)),
        SparseMatEntry::new(1, 0, Scalar::from_u64(6)),
        SparseMatEntry::new(2, 1, Scalar::from_u64(7)),
        SparseMatEntry::new(3, 2, Scalar::from_u64(8)),
    ];

    let mat_a = SparseMatPolynomial::new(num_vars_x, num_vars_y, sparse_a);
    let mat_b = SparseMatPolynomial::new(num_vars_x, num_vars_y, sparse_b);
    let mat_c = SparseMatPolynomial::new(num_vars_x, num_vars_y, sparse_c);

    println!("📊 Matrix sizes:");
    println!("  A: {} non-zero entries", mat_a.get_num_nz_entries());
    println!("  B: {} non-zero entries", mat_b.get_num_nz_entries());
    println!("  C: {} non-zero entries", mat_c.get_num_nz_entries());

    // Create commitment generators
    println!("\n🔧 Creating generators...");
    let start = Instant::now();
    let gens = SparseMatPolyCommitmentGens::new(
        b"test-sparse",
        num_vars_x,
        num_vars_y,
        4, // max non-zero entries
        3, // batch size (A, B, C)
    );
    println!("  Gens time: {:?}", start.elapsed());

    // Commit to matrices
    println!("\n📦 Committing to matrices...");
    let start = Instant::now();
    let (comm, dense) = SparseMatPolynomial::multi_commit(&[&mat_a, &mat_b, &mat_c], &gens);
    println!("  Commit time: {:?}", start.elapsed());

    // Random evaluation point
    let rx = vec![Scalar::from_u64(17), Scalar::from_u64(23)];
    let ry = vec![Scalar::from_u64(31), Scalar::from_u64(37)];

    // Compute expected evaluations
    let rx_evals = EqPolynomial::new(rx.clone()).evals();
    let ry_evals = EqPolynomial::new(ry.clone()).evals();
    
    let eval_a = mat_a.evaluate_with_tables(&rx_evals, &ry_evals);
    let eval_b = mat_b.evaluate_with_tables(&rx_evals, &ry_evals);
    let eval_c = mat_c.evaluate_with_tables(&rx_evals, &ry_evals);

    println!("\n📊 Evaluations:");
    println!("  A(rx, ry) = {:?}", eval_a);
    println!("  B(rx, ry) = {:?}", eval_b);
    println!("  C(rx, ry) = {:?}", eval_c);

    // Create proof
    println!("\n🚀 Generating evaluation proof...");
    let mut prover_transcript = Transcript::new(b"sparse-poly-test");
    let mut random_tape = RandomTape::new(b"sparse-poly-tape");
    
    let start = Instant::now();
    let proof = SparseMatPolyEvalProof::prove(
        &dense,
        &rx,
        &ry,
        &[eval_a, eval_b, eval_c],
        &gens,
        &mut prover_transcript,
        &mut random_tape,
    );
    let prove_time = start.elapsed();
    println!("  Prove time: {:?}", prove_time);

    // Serialize to measure size
    let proof_bytes = bincode::serialize(&proof)?;
    println!("\n📦 Proof size: {} bytes ({:.1} KB)", proof_bytes.len(), proof_bytes.len() as f64 / 1024.0);

    // Verify
    println!("\n🔍 Verifying...");
    let mut verifier_transcript = Transcript::new(b"sparse-poly-test");
    
    let start = Instant::now();
    proof.verify(
        &comm,
        &rx,
        &ry,
        &[eval_a, eval_b, eval_c],
        &gens,
        &mut verifier_transcript,
    )?;
    let verify_time = start.elapsed();
    println!("  ✅ Verification passed ({:?})", verify_time);

    println!("\n╔══════════════════════════════════════════════════════════════════════╗");
    println!("║                           SUCCESS                                     ║");
    println!("╚══════════════════════════════════════════════════════════════════════╝");

    Ok(())
}

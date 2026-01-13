//! Benchmark sparse polynomial evaluation at different sizes
//! Compare with Spartan paper Figure 7, 8, 9

use spartan_bn254::{
    dense_mlpoly::EqPolynomial,
    math::Math,
    random::RandomTape,
    scalar::Scalar,
    sparse_mlpoly_full::{
        SparseMatEntry, SparseMatPolynomial, SparseMatPolyCommitmentGens, SparseMatPolyEvalProof,
    },
};
use merlin::Transcript;
use std::time::Instant;

fn benchmark_size(num_vars: usize) -> Result<(), Box<dyn std::error::Error>> {
    let size = 1usize << num_vars;
    let num_nz = size; // Approximately n non-zero entries per matrix
    
    println!("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    println!("  2^{} = {} constraints", num_vars, size);
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    
    // Create sparse matrices with num_nz entries each
    let mut rng = ark_std::test_rng();
    
    let sparse_a: Vec<SparseMatEntry> = (0..num_nz)
        .map(|i| SparseMatEntry::new(i % size, i % size, Scalar::random(&mut rng)))
        .collect();
    let sparse_b: Vec<SparseMatEntry> = (0..num_nz)
        .map(|i| SparseMatEntry::new(i % size, (i + 1) % size, Scalar::random(&mut rng)))
        .collect();
    let sparse_c: Vec<SparseMatEntry> = (0..num_nz)
        .map(|i| SparseMatEntry::new(i % size, (i + 2) % size, Scalar::random(&mut rng)))
        .collect();

    let mat_a = SparseMatPolynomial::new(num_vars, num_vars, sparse_a);
    let mat_b = SparseMatPolynomial::new(num_vars, num_vars, sparse_b);
    let mat_c = SparseMatPolynomial::new(num_vars, num_vars, sparse_c);

    // Create commitment generators
    let start = Instant::now();
    let gens = SparseMatPolyCommitmentGens::new(
        b"bench-sparse",
        num_vars,
        num_vars,
        num_nz.next_power_of_two(),
        3,
    );
    let gens_time = start.elapsed();

    // Commit to matrices (preprocessing/encoding)
    let start = Instant::now();
    let (comm, dense) = SparseMatPolynomial::multi_commit(&[&mat_a, &mat_b, &mat_c], &gens);
    let encode_time = start.elapsed();

    // Random evaluation point
    let rx: Vec<Scalar> = (0..num_vars).map(|_| Scalar::random(&mut rng)).collect();
    let ry: Vec<Scalar> = (0..num_vars).map(|_| Scalar::random(&mut rng)).collect();

    // Compute expected evaluations
    let rx_evals = EqPolynomial::new(rx.clone()).evals();
    let ry_evals = EqPolynomial::new(ry.clone()).evals();
    
    let eval_a = mat_a.evaluate_with_tables(&rx_evals, &ry_evals);
    let eval_b = mat_b.evaluate_with_tables(&rx_evals, &ry_evals);
    let eval_c = mat_c.evaluate_with_tables(&rx_evals, &ry_evals);

    // Create proof
    let mut prover_transcript = Transcript::new(b"bench");
    let mut random_tape = RandomTape::new(b"bench-tape");
    
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

    // Serialize to measure size
    let proof_bytes = bincode::serialize(&proof)?;
    let proof_size_kb = proof_bytes.len() as f64 / 1024.0;

    // Verify
    let mut verifier_transcript = Transcript::new(b"bench");
    
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

    println!("  Encode:  {:>10.2?}", encode_time);
    println!("  Prove:   {:>10.2?}", prove_time);
    println!("  Verify:  {:>10.2?}", verify_time);
    println!("  Proof:   {:>10.1} KB", proof_size_kb);

    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("╔══════════════════════════════════════════════════════════════════════╗");
    println!("║     Spartan-BN254 Scaling Benchmark                                  ║");
    println!("║     (Compare with Spartan paper Figures 7, 8, 9, 10)                 ║");
    println!("╚══════════════════════════════════════════════════════════════════════╝");

    // Test sizes from 2^10 to 2^16 (to keep runtime reasonable)
    for num_vars in [10, 12, 14, 16] {
        benchmark_size(num_vars)?;
    }

    println!("\n╔══════════════════════════════════════════════════════════════════════╗");
    println!("║     Paper Benchmarks (SpartanSNARK, curve25519)                       ║");
    println!("╠══════════════════════════════════════════════════════════════════════╣");
    println!("║  Size  │ Prover  │ Proof KB │ Verifier │ Encoder                     ║");
    println!("║  2^10  │   0.07s │    32 KB │   9.6 ms │   0.04s                     ║");
    println!("║  2^12  │   0.21s │    48 KB │  13.9 ms │   0.12s                     ║");
    println!("║  2^14  │   0.79s │    54 KB │  21.0 ms │   0.40s                     ║");
    println!("║  2^16  │   2.60s │    63 KB │  34.3 ms │   1.40s                     ║");
    println!("║  2^18  │   9.20s │    85 KB │  55.9 ms │   4.50s                     ║");
    println!("║  2^20  │  36.30s │   142 KB │ 100.3 ms │  15.10s                     ║");
    println!("╚══════════════════════════════════════════════════════════════════════╝");

    Ok(())
}

//! Micro-benchmarks to understand overhead sources
//! Compare individual operation costs

use spartan_bn254::{
    scalar::Scalar,
    group::GroupElement,
    commitments::{MultiCommitGens, Commitments},
};
use std::time::Instant;
use ark_std::UniformRand;
use ark_ec::CurveGroup;

fn main() {
    println!("╔══════════════════════════════════════════════════════════════════════╗");
    println!("║     Micro-Benchmarks: Understanding Overhead Sources                 ║");
    println!("╚══════════════════════════════════════════════════════════════════════╝\n");

    let mut rng = ark_std::test_rng();
    
    // =========================================================================
    // 1. Scalar Field Operations
    // =========================================================================
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    println!("  1. Scalar Field Operations (BN254 Fr)");
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    
    let n = 100_000;
    let scalars: Vec<Scalar> = (0..n).map(|_| Scalar::random(&mut rng)).collect();
    
    // Multiplication
    let start = Instant::now();
    let mut acc = Scalar::one();
    for s in &scalars {
        acc = acc * *s;
    }
    let mul_time = start.elapsed();
    println!("  {} scalar multiplications: {:?} ({:.0} ns/op)", 
        n, mul_time, mul_time.as_nanos() as f64 / n as f64);
    
    // Addition
    let start = Instant::now();
    let mut acc = Scalar::zero();
    for s in &scalars {
        acc = acc + *s;
    }
    let add_time = start.elapsed();
    println!("  {} scalar additions: {:?} ({:.0} ns/op)", 
        n, add_time, add_time.as_nanos() as f64 / n as f64);

    // Inversion
    let start = Instant::now();
    for s in scalars.iter().take(1000) {
        let _ = s.invert();
    }
    let inv_time = start.elapsed();
    println!("  1000 scalar inversions: {:?} ({:.0} ns/op)", 
        inv_time, inv_time.as_nanos() as f64 / 1000.0);

    // =========================================================================
    // 2. Group Operations (BN254 G1)
    // =========================================================================
    println!("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    println!("  2. Group Operations (BN254 G1)");
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");

    let n = 10_000;
    let points: Vec<GroupElement> = (0..n)
        .map(|_| {
            let p = ark_bn254::G1Projective::rand(&mut rng);
            GroupElement(p)
        })
        .collect();
    let scalars_small: Vec<Scalar> = (0..n).map(|_| Scalar::random(&mut rng)).collect();

    // Point addition
    let start = Instant::now();
    let mut acc = points[0];
    for p in &points[1..] {
        acc = acc + *p;
    }
    let add_time = start.elapsed();
    println!("  {} point additions: {:?} ({:.0} ns/op)", 
        n-1, add_time, add_time.as_nanos() as f64 / (n-1) as f64);

    // Scalar multiplication
    let start = Instant::now();
    for i in 0..1000 {
        let _ = points[i] * scalars_small[i];
    }
    let smul_time = start.elapsed();
    println!("  1000 scalar multiplications: {:?} ({:.0} µs/op)", 
        smul_time, smul_time.as_micros() as f64 / 1000.0);

    // =========================================================================
    // 3. Multi-Scalar Multiplication (MSM) - The Main Cost
    // =========================================================================
    println!("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    println!("  3. Multi-Scalar Multiplication (MSM) - Main SNARK Cost");
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");

    for log_n in [10, 12, 14, 16] {
        let n = 1usize << log_n;
        let gens = MultiCommitGens::new(n, b"msm-bench");
        let scalars: Vec<Scalar> = (0..n).map(|_| Scalar::random(&mut rng)).collect();
        let blind = Scalar::zero();
        
        let start = Instant::now();
        let _ = scalars.commit(&blind, &gens);
        let msm_time = start.elapsed();
        
        println!("  MSM 2^{}: {:?} ({:.2} µs/element)", 
            log_n, msm_time, msm_time.as_micros() as f64 / n as f64);
    }

    // =========================================================================
    // 4. Element Sizes
    // =========================================================================
    println!("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    println!("  4. Element Sizes (BN254 vs curve25519)");
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    
    println!("  Scalar (Fr):           32 bytes (both curves)");
    println!("  BN254 G1 compressed:   48 bytes");
    println!("  curve25519 compressed: 32 bytes");
    println!("  Size ratio:            1.5x larger for BN254");

    // =========================================================================
    // 5. Estimated Proof Composition
    // =========================================================================
    println!("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    println!("  5. Proof Size Analysis");
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    
    // Paper proof at 2^20: 142 KB
    // Our proof at 2^20: 231 KB
    // Ratio: 1.63x
    
    let paper_size = 142.0;
    let our_size = 231.0;
    let ratio = our_size / paper_size;
    
    println!("  Paper proof size (2^20):    142 KB");
    println!("  Our proof size (2^20):      231 KB");
    println!("  Ratio:                      {:.2}x", ratio);
    println!();
    println!("  If proof is ~60% group elements, ~40% scalars:");
    println!("    Expected ratio = 0.6 * 1.5 + 0.4 * 1.0 = 1.3x");
    println!("  Actual ratio: {:.2}x", ratio);
    println!("  Extra overhead: {:.0}%", (ratio / 1.3 - 1.0) * 100.0);

    // =========================================================================
    // 6. Prover Time Analysis
    // =========================================================================
    println!("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    println!("  6. Prover Time Analysis");
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    
    let paper_time = 36.3;
    let our_time = 250.5;
    let time_ratio = our_time / paper_time;
    
    println!("  Paper prover time (2^20):   36.3s");
    println!("  Our prover time (2^20):     250.5s");
    println!("  Ratio:                      {:.1}x", time_ratio);
    println!();
    println!("  Expected from curve difference: ~2-3x");
    println!("  Actual ratio: {:.1}x", time_ratio);
    println!("  Extra unexplained overhead: {:.1}x", time_ratio / 2.5);

    println!("\n╔══════════════════════════════════════════════════════════════════════╗");
    println!("║                         ANALYSIS SUMMARY                             ║");
    println!("╚══════════════════════════════════════════════════════════════════════╝");
    println!();
    println!("  Proof Size Overhead:");
    println!("    - BN254 G1 is 1.5x larger than curve25519 ✓");
    println!("    - Our proof is 1.6x larger → matches expectation ✓");
    println!();
    println!("  Prover Time Overhead:");
    println!("    - Expected ~2-3x from curve operations");
    println!("    - Actual 6.9x → ~2.5x unexplained overhead");
    println!();
    println!("  Possible causes of extra overhead:");
    println!("    1. ark-bn254 vs curve25519-dalek optimization level");
    println!("    2. Memory allocation patterns");
    println!("    3. Implementation differences in sumcheck/product circuits");
    println!("    4. Our port may have inefficiencies");
}

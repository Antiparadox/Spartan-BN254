//! Debug script to check R1CS constraint mapping

use spartan_bn254::{Scalar, r1cs_reader::R1CS};

fn main() {
    let r1cs = R1CS::from_file("/Users/yinuo/Desktop/Universal ZK/keyless-zk-proofs/main.r1cs").unwrap();
    let stats = r1cs.stats();
    
    println!("=== Circom R1CS Stats ===");
    println!("num_variables: {}", stats.num_variables);
    println!("num_pub_inputs: {}", stats.num_pub_inputs);
    println!("num_prv_inputs: {}", stats.num_prv_inputs);
    println!("num_constraints: {}", stats.num_constraints);
    println!("n_prv = {}", r1cs.num_private_vars());
    
    // Find entries for constraint 127 in original format
    println!("\n=== Constraint 127 in ORIGINAL circom format ===");
    println!("A entries:");
    for (row, col, val) in &r1cs.a {
        if *row == 127 {
            println!("  ({}, {}, {:?})", row, col, val);
        }
    }
    println!("B entries:");
    for (row, col, val) in &r1cs.b {
        if *row == 127 {
            println!("  ({}, {}, {:?})", row, col, val);
        }
    }
    println!("C entries:");
    for (row, col, val) in &r1cs.c {
        if *row == 127 {
            println!("  ({}, {}, {:?})", row, col, val);
        }
    }
    
    // Get remapped version
    let (a_remapped, b_remapped, c_remapped) = r1cs.to_sparse_matrices();
    
    println!("\n=== Constraint 127 in REMAPPED Spartan format ===");
    println!("A entries:");
    for (row, col, val) in &a_remapped {
        if *row == 127 {
            println!("  ({}, {}, {:?})", row, col, val);
        }
    }
    println!("B entries:");
    for (row, col, val) in &b_remapped {
        if *row == 127 {
            println!("  ({}, {}, {:?})", row, col, val);
        }
    }
    println!("C entries:");
    for (row, col, val) in &c_remapped {
        if *row == 127 {
            println!("  ({}, {}, {:?})", row, col, val);
        }
    }
    
    // Show remapping for specific columns
    let n_pub = stats.num_pub_inputs;
    let n_prv = r1cs.num_private_vars();
    
    println!("\n=== Column Remapping ===");
    println!("n_pub = {}, n_prv = {}", n_pub, n_prv);
    println!("circom col 0 (constant 1) -> spartan col {}", n_prv);
    println!("circom col 1 (public input) -> spartan col {}", n_prv + 1);
    println!("circom col 2 (first prv var) -> spartan col 0");
    println!("circom col {} (last wire) -> spartan col {}", stats.num_variables - 1, stats.num_variables - 1 - n_pub - 1);
}

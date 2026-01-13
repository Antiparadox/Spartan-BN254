//! Simple example demonstrating Spartan-BN254 basic functionality

use spartan_bn254::{DensePolynomial, R1CSShape, Scalar};

fn main() {
    println!("Spartan-BN254 Simple Example");
    println!("=============================\n");

    // Test scalar operations
    let a = Scalar::from_u64(5);
    let b = Scalar::from_u64(7);
    let c = a * b;
    println!("Scalar arithmetic: 5 * 7 = {:?}", c);

    // Test dense polynomial
    let Z = vec![
        Scalar::from_u64(1),
        Scalar::from_u64(2),
        Scalar::from_u64(3),
        Scalar::from_u64(4),
    ];
    let poly = DensePolynomial::new(Z);
    println!("\nDense polynomial:");
    println!("  - Length: {}", poly.len());
    println!("  - Num vars: {}", poly.get_num_vars());

    // Evaluate polynomial at a random point
    let r = vec![Scalar::from_u64(2), Scalar::from_u64(3)];
    let eval = poly.evaluate(&r);
    println!("  - Evaluation at (2, 3): {:?}", eval);

    // Test R1CS
    println!("\nR1CS test:");
    
    // Simple constraint: x0 * x0 = x1 (x0^2 = x1)
    // With x0 = 3, x1 = 9, this should be satisfied
    let A = vec![(0, 0, Scalar::one())]; // A[0,0] = 1 (selects x0)
    let B = vec![(0, 0, Scalar::one())]; // B[0,0] = 1 (selects x0)
    let C = vec![(0, 1, Scalar::one())]; // C[0,1] = 1 (selects x1)

    // Need power-of-2 dimensions
    let num_cons: usize = 1; // Will be padded to 2
    let num_vars: usize = 2; // Will be padded to 2
    let num_inputs: usize = 0;

    // Pad to power of 2
    let num_cons_padded = num_cons.next_power_of_two().max(2);
    let num_vars_padded = num_vars.next_power_of_two();

    let shape = R1CSShape::new(num_cons_padded, num_vars_padded, num_inputs, &A, &B, &C);
    
    println!("  - Num constraints (padded): {}", shape.get_num_cons());
    println!("  - Num variables (padded): {}", shape.get_num_vars());

    // Test satisfaction with x0=3, x1=9
    let mut vars = vec![Scalar::zero(); num_vars_padded];
    vars[0] = Scalar::from_u64(3);
    vars[1] = Scalar::from_u64(9);
    
    let inputs: Vec<Scalar> = vec![];
    let is_sat = shape.is_sat(&vars, &inputs);
    println!("  - Constraint x0^2 = x1 with x0=3, x1=9: {}", 
             if is_sat { "SATISFIED ✓" } else { "NOT SATISFIED ✗" });

    // Test with wrong values
    vars[1] = Scalar::from_u64(10); // Wrong!
    let is_sat = shape.is_sat(&vars, &inputs);
    println!("  - Constraint x0^2 = x1 with x0=3, x1=10: {}", 
             if is_sat { "SATISFIED ✓" } else { "NOT SATISFIED ✗" });

    println!("\n✅ Spartan-BN254 basic functionality works!");
}

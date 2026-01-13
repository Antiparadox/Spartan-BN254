//! Spartan zkSNARK adapted for BN254 curve
//!
//! This is a port of Microsoft's Spartan to use the BN254 curve,
//! making it compatible with circom R1CS circuits.
//!
//! # SNARK Mode
//! Use `SNARKGens`, `SNARK`, and `SNARK::encode`/`prove`/`verify` for
//! succinct verification with preprocessing.
//!
//! # NIZK Mode  
//! Use `NIZKGens` and `NIZK` for non-interactive zero-knowledge proofs
//! without preprocessing (verifier needs full R1CS).

#![allow(non_snake_case)]

pub mod commitments;
pub mod dense_mlpoly;
pub mod errors;
pub mod group;
pub mod math;
pub mod nizk;
pub mod product_tree;
pub mod r1cs;
pub mod r1cs_reader;
pub mod r1csproof;
pub mod random;
pub mod scalar;
pub mod snark;
pub mod sparse_mlpoly;
pub mod sparse_mlpoly_full;
pub mod sumcheck;
pub mod timer;
pub mod transcript;
pub mod unipoly;

// Re-exports
pub use commitments::{Commitments, MultiCommitGens};
pub use dense_mlpoly::{DensePolynomial, EqPolynomial};
pub use errors::{ProofVerifyError, R1CSError};
pub use group::{CompressedGroup, GroupElement};
pub use nizk::{DotProductProofGens, DotProductProofLog, EqualityProof, KnowledgeProof, ProductProof};
pub use r1cs::{R1CSInstance, R1CSShape};
pub use r1cs_reader::R1CS;
pub use r1csproof::{R1CSGens, R1CSProof, PolyCommitmentGens};
pub use random::RandomTape;
pub use scalar::Scalar;
pub use snark::{Assignment, InputsAssignment, Instance, NIZK, NIZKGens, SNARK, SNARKGens, VarsAssignment};
pub use sumcheck::{SumcheckInstanceProof, ZKSumcheckInstanceProof};
pub use unipoly::{CompressedUniPoly, UniPoly};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_scalar_operations() {
        let a = Scalar::from_u64(5);
        let b = Scalar::from_u64(7);
        assert_eq!(a + b, Scalar::from_u64(12));
        assert_eq!(a * b, Scalar::from_u64(35));
    }

    #[test]
    fn test_r1cs_reader_integration() {
        // This will be tested with actual .r1cs files
    }

    #[test]
    fn test_dense_poly() {
        let Z = vec![
            Scalar::from_u64(1),
            Scalar::from_u64(2),
            Scalar::from_u64(3),
            Scalar::from_u64(4),
        ];
        let poly = DensePolynomial::new(Z);
        assert_eq!(poly.len(), 4);
        assert_eq!(poly.get_num_vars(), 2);
    }
}

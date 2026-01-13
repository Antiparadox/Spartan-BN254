//! Pedersen commitments on BN254 G1
//! Port of Spartan's commitments.rs

use crate::group::{group_basepoint_compressed, CompressedGroup, GroupElement};
use crate::scalar::Scalar;
use serde::{Deserialize, Serialize};
use sha3::{Digest, Shake256, digest::{ExtendableOutput, XofReader}};

/// Generators for multi-scalar commitments
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MultiCommitGens {
    pub n: usize,
    pub G: Vec<GroupElement>,
    pub h: GroupElement,
}

impl MultiCommitGens {
    /// Create new generators by hashing a label
    pub fn new(n: usize, label: &[u8]) -> Self {
        use sha3::digest::Update;
        
        let mut shake = Shake256::default();
        shake.update(label);
        shake.update(group_basepoint_compressed().as_bytes());

        let mut reader = shake.finalize_xof();
        let mut gens: Vec<GroupElement> = Vec::new();
        let mut uniform_bytes = [0u8; 64];
        
        for _ in 0..n + 1 {
            reader.read(&mut uniform_bytes);
            gens.push(GroupElement::from_uniform_bytes(&uniform_bytes));
        }

        MultiCommitGens {
            n,
            G: gens[..n].to_vec(),
            h: gens[n],
        }
    }

    pub fn scale(&self, s: &Scalar) -> MultiCommitGens {
        MultiCommitGens {
            n: self.n,
            h: self.h,
            G: (0..self.n).map(|i| *s * &self.G[i]).collect(),
        }
    }

    pub fn split_at(&self, mid: usize) -> (MultiCommitGens, MultiCommitGens) {
        let (G1, G2) = self.G.split_at(mid);

        (
            MultiCommitGens {
                n: G1.len(),
                G: G1.to_vec(),
                h: self.h,
            },
            MultiCommitGens {
                n: G2.len(),
                G: G2.to_vec(),
                h: self.h,
            },
        )
    }
}

/// Trait for Pedersen commitments
pub trait Commitments {
    fn commit(&self, blind: &Scalar, gens_n: &MultiCommitGens) -> GroupElement;
}

impl Commitments for Scalar {
    fn commit(&self, blind: &Scalar, gens_n: &MultiCommitGens) -> GroupElement {
        assert_eq!(gens_n.n, 1);
        GroupElement::vartime_multiscalar_mul(&[*self, *blind], &[gens_n.G[0], gens_n.h])
    }
}

impl Commitments for Vec<Scalar> {
    fn commit(&self, blind: &Scalar, gens_n: &MultiCommitGens) -> GroupElement {
        assert_eq!(gens_n.n, self.len());
        GroupElement::vartime_multiscalar_mul(self, &gens_n.G) + *blind * gens_n.h
    }
}

impl Commitments for [Scalar] {
    fn commit(&self, blind: &Scalar, gens_n: &MultiCommitGens) -> GroupElement {
        assert_eq!(gens_n.n, self.len());
        GroupElement::vartime_multiscalar_mul(self, &gens_n.G) + *blind * gens_n.h
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_commitment_gens() {
        let gens = MultiCommitGens::new(4, b"test");
        assert_eq!(gens.n, 4);
        assert_eq!(gens.G.len(), 4);
    }

    #[test]
    fn test_scalar_commit() {
        let gens = MultiCommitGens::new(1, b"test");
        let val = Scalar::from_u64(42);
        let blind = Scalar::from_u64(123);
        let commit = val.commit(&blind, &gens);
        
        // Commitment should be non-trivial
        assert_ne!(commit, GroupElement::identity());
    }

    #[test]
    fn test_vector_commit() {
        let gens = MultiCommitGens::new(3, b"test");
        let vals = vec![Scalar::from_u64(1), Scalar::from_u64(2), Scalar::from_u64(3)];
        let blind = Scalar::from_u64(456);
        let commit = vals.commit(&blind, &gens);
        
        assert_ne!(commit, GroupElement::identity());
    }
}

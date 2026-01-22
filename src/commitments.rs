//! Pedersen commitments on BN254 G1
//! Port of Spartan's commitments.rs
//! 
//! Optimization: Stores generators as affine points to avoid repeated conversions during MSM.

use ark_bn254::G1Affine;
use ark_ec::CurveGroup;
use crate::group::{group_basepoint_compressed, CompressedGroup, GroupElement};
use crate::scalar::Scalar;
use serde::{Deserialize, Serialize};
use sha3::{Digest, Shake256, digest::{ExtendableOutput, XofReader}};

/// Generators for multi-scalar commitments
/// 
/// Stores generators in AFFINE form for faster MSM (avoids projective->affine conversion each call)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MultiCommitGens {
    pub n: usize,
    /// Generators stored as projective points (for compatibility)
    pub G: Vec<GroupElement>,
    pub h: GroupElement,
    /// Pre-computed affine generators for fast MSM (skips conversion)
    #[serde(skip)]
    pub G_affine: Vec<G1Affine>,
    #[serde(skip)]
    pub h_affine: G1Affine,
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

        let G = gens[..n].to_vec();
        let h = gens[n];
        
        // Pre-compute affine points for fast MSM
        let G_projective: Vec<_> = G.iter().map(|g| *g.inner()).collect();
        let G_affine = <ark_bn254::G1Projective as CurveGroup>::normalize_batch(&G_projective);
        let h_affine = h.inner().into_affine();

        MultiCommitGens {
            n,
            G,
            h,
            G_affine,
            h_affine,
        }
    }

    pub fn scale(&self, s: &Scalar) -> MultiCommitGens {
        let G: Vec<GroupElement> = (0..self.n).map(|i| *s * &self.G[i]).collect();
        let G_projective: Vec<_> = G.iter().map(|g| *g.inner()).collect();
        let G_affine = <ark_bn254::G1Projective as CurveGroup>::normalize_batch(&G_projective);
        
        MultiCommitGens {
            n: self.n,
            h: self.h,
            h_affine: self.h_affine,
            G,
            G_affine,
        }
    }

    pub fn split_at(&self, mid: usize) -> (MultiCommitGens, MultiCommitGens) {
        let (G1, G2) = self.G.split_at(mid);
        let (G1_affine, G2_affine) = self.G_affine.split_at(mid);

        (
            MultiCommitGens {
                n: G1.len(),
                G: G1.to_vec(),
                h: self.h,
                G_affine: G1_affine.to_vec(),
                h_affine: self.h_affine,
            },
            MultiCommitGens {
                n: G2.len(),
                G: G2.to_vec(),
                h: self.h,
                G_affine: G2_affine.to_vec(),
                h_affine: self.h_affine,
            },
        )
    }

    /// Create from existing generators (computes affine cache automatically)
    pub fn from_generators(G: Vec<GroupElement>, h: GroupElement) -> Self {
        let n = G.len();
        let G_projective: Vec<_> = G.iter().map(|g| *g.inner()).collect();
        let G_affine = <ark_bn254::G1Projective as CurveGroup>::normalize_batch(&G_projective);
        let h_affine = h.inner().into_affine();
        
        MultiCommitGens {
            n,
            G,
            h,
            G_affine,
            h_affine,
        }
    }
}

/// Trait for Pedersen commitments
pub trait Commitments {
    fn commit(&self, blind: &Scalar, gens_n: &MultiCommitGens) -> GroupElement;
}

impl Commitments for Scalar {
    fn commit(&self, blind: &Scalar, gens_n: &MultiCommitGens) -> GroupElement {
        assert_eq!(gens_n.n, 1);
        // Use affine generators for faster MSM
        let scalars = [*self, *blind];
        let points = [gens_n.G_affine[0], gens_n.h_affine];
        GroupElement::msm_affine(&scalars, &points)
    }
}

impl Commitments for Vec<Scalar> {
    fn commit(&self, blind: &Scalar, gens_n: &MultiCommitGens) -> GroupElement {
        assert_eq!(gens_n.n, self.len());
        // Use affine generators for faster MSM
        let mut scalars = self.clone();
        scalars.push(*blind);
        let mut points = gens_n.G_affine.clone();
        points.push(gens_n.h_affine);
        GroupElement::msm_affine(&scalars, &points)
    }
}

impl Commitments for [Scalar] {
    fn commit(&self, blind: &Scalar, gens_n: &MultiCommitGens) -> GroupElement {
        assert_eq!(gens_n.n, self.len());
        // Use affine generators for faster MSM
        let mut scalars = self.to_vec();
        scalars.push(*blind);
        let mut points = gens_n.G_affine.clone();
        points.push(gens_n.h_affine);
        GroupElement::msm_affine(&scalars, &points)
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

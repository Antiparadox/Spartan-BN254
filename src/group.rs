//! Group operations for BN254 G1
//! Port of Spartan's group.rs from Ristretto255 to BN254

use crate::scalar::Scalar;
use ark_bn254::{Fr, G1Affine, G1Projective};
use ark_ec::{AffineRepr, CurveGroup, PrimeGroup, VariableBaseMSM};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use core::borrow::Borrow;
use core::ops::{Add, Mul};
use serde::{Deserialize, Serialize};
use sha3::{Digest, Sha3_256};

/// Group element (BN254 G1 projective point)
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct GroupElement(pub G1Projective);

/// Compressed group element (BN254 G1 affine point serialized)
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CompressedGroup(pub Vec<u8>);

impl GroupElement {
    /// Identity element
    pub fn identity() -> Self {
        GroupElement(G1Projective::default())
    }

    /// Generator
    pub fn generator() -> Self {
        GroupElement(G1Projective::generator())
    }

    /// Create from uniform random bytes (hash-to-curve)
    pub fn from_uniform_bytes(bytes: &[u8; 64]) -> Self {
        // Simple hash-to-curve: hash bytes to scalar, multiply generator
        // Note: This is a simplified version; production should use proper hash-to-curve
        let mut hasher = Sha3_256::new();
        hasher.update(bytes);
        let hash = hasher.finalize();
        
        let mut scalar_bytes = [0u8; 32];
        scalar_bytes.copy_from_slice(&hash[..32]);
        
        if let Some(scalar) = Scalar::from_bytes(&scalar_bytes) {
            GroupElement(G1Projective::generator() * scalar.0)
        } else {
            // Fallback: use a different derivation
            let mut hasher2 = Sha3_256::new();
            hasher2.update(b"fallback");
            hasher2.update(bytes);
            let hash2 = hasher2.finalize();
            scalar_bytes.copy_from_slice(&hash2[..32]);
            let scalar = Scalar::from_bytes(&scalar_bytes).unwrap_or(Scalar::one());
            GroupElement(G1Projective::generator() * scalar.0)
        }
    }

    /// Compress to bytes
    pub fn compress(&self) -> CompressedGroup {
        let affine = self.0.into_affine();
        let mut bytes = Vec::new();
        affine.serialize_compressed(&mut bytes).unwrap();
        CompressedGroup(bytes)
    }

    /// Variable-time multi-scalar multiplication
    pub fn vartime_multiscalar_mul<I, J>(scalars: I, points: J) -> Self
    where
        I: IntoIterator,
        I::Item: Borrow<Scalar>,
        J: IntoIterator,
        J::Item: Borrow<Self>,
    {
        let scalars_vec: Vec<Fr> = scalars.into_iter().map(|s| s.borrow().0).collect();
        let points_vec: Vec<G1Affine> = points
            .into_iter()
            .map(|p| p.borrow().0.into_affine())
            .collect();

        let result = G1Projective::msm(&points_vec, &scalars_vec).unwrap_or_default();
        GroupElement(result)
    }
}

impl CompressedGroup {
    /// Create from raw bytes
    pub fn from_bytes(bytes: &[u8]) -> Self {
        CompressedGroup(bytes.to_vec())
    }
    
    /// Decompress to group element
    pub fn decompress(&self) -> Option<GroupElement> {
        G1Affine::deserialize_compressed(&self.0[..])
            .ok()
            .map(|affine| GroupElement(affine.into()))
    }

    /// Get bytes
    pub fn as_bytes(&self) -> &[u8] {
        &self.0
    }

    /// Convert to bytes
    pub fn to_bytes(&self) -> Vec<u8> {
        self.0.clone()
    }
}

impl Default for GroupElement {
    fn default() -> Self {
        Self::identity()
    }
}

impl Add for GroupElement {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        GroupElement(self.0 + rhs.0)
    }
}

impl<'a> Add<&'a GroupElement> for GroupElement {
    type Output = Self;
    fn add(self, rhs: &'a GroupElement) -> Self {
        GroupElement(self.0 + rhs.0)
    }
}

impl Mul<Scalar> for GroupElement {
    type Output = Self;
    fn mul(self, rhs: Scalar) -> Self {
        GroupElement(self.0 * rhs.0)
    }
}

impl<'a> Mul<&'a Scalar> for GroupElement {
    type Output = Self;
    fn mul(self, rhs: &'a Scalar) -> Self {
        GroupElement(self.0 * rhs.0)
    }
}

impl<'a> Mul<&'a Scalar> for &'a GroupElement {
    type Output = GroupElement;
    fn mul(self, rhs: &'a Scalar) -> GroupElement {
        GroupElement(self.0 * rhs.0)
    }
}

impl Mul<GroupElement> for Scalar {
    type Output = GroupElement;
    fn mul(self, rhs: GroupElement) -> GroupElement {
        GroupElement(rhs.0 * self.0)
    }
}

impl<'a> Mul<&'a GroupElement> for Scalar {
    type Output = GroupElement;
    fn mul(self, rhs: &'a GroupElement) -> GroupElement {
        GroupElement(rhs.0 * self.0)
    }
}

impl<'a> Mul<&'a GroupElement> for &'a Scalar {
    type Output = GroupElement;
    fn mul(self, rhs: &'a GroupElement) -> GroupElement {
        GroupElement(rhs.0 * self.0)
    }
}

// Serde support
impl Serialize for GroupElement {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let compressed = self.compress();
        compressed.0.serialize(serializer)
    }
}

impl<'de> Deserialize<'de> for GroupElement {
    fn deserialize<D: serde::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        let bytes: Vec<u8> = Deserialize::deserialize(deserializer)?;
        let compressed = CompressedGroup(bytes);
        compressed
            .decompress()
            .ok_or_else(|| serde::de::Error::custom("invalid group element"))
    }
}

impl Serialize for CompressedGroup {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        self.0.serialize(serializer)
    }
}

impl<'de> Deserialize<'de> for CompressedGroup {
    fn deserialize<D: serde::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        let bytes: Vec<u8> = Deserialize::deserialize(deserializer)?;
        Ok(CompressedGroup(bytes))
    }
}

/// Generator point (compressed)
pub fn group_basepoint_compressed() -> CompressedGroup {
    GroupElement::generator().compress()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_group_operations() {
        let g = GroupElement::generator();
        let two = Scalar::from_u64(2);
        let g2 = g * two;
        let g_plus_g = g + g;
        assert_eq!(g2, g_plus_g);
    }

    #[test]
    fn test_msm() {
        let g = GroupElement::generator();
        let scalars = vec![Scalar::from_u64(2), Scalar::from_u64(3)];
        let points = vec![g, g];
        let result = GroupElement::vartime_multiscalar_mul(&scalars, &points);
        let expected = g * Scalar::from_u64(5);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_compress_decompress() {
        let g = GroupElement::generator();
        let compressed = g.compress();
        let decompressed = compressed.decompress().unwrap();
        assert_eq!(g, decompressed);
    }
}

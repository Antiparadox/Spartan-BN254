//! Group operations for BN254 G1
//! Port of Spartan's group.rs from Ristretto255 to BN254
//!
//! Implements CanonicalSerialize/CanonicalDeserialize for cross-verification
//! compatibility with arkworks-spartan.

use crate::scalar::Scalar;
use ark_bn254::{Fr, G1Affine, G1Projective};
use ark_ec::{AffineRepr, CurveGroup, PrimeGroup, VariableBaseMSM};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Compress, Valid, Validate};
use core::borrow::Borrow;
use core::ops::{Add, Mul};
use serde::{Deserialize, Serialize};
use sha3::{Digest, Sha3_256};

/// Group element (BN254 G1 projective point)
/// 
/// Serializes as G1Projective for arkworks compatibility
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct GroupElement(pub G1Projective);

/// Compressed group element (BN254 G1 affine point serialized)
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CompressedGroup(pub Vec<u8>);

// CanonicalSerialize for GroupElement: serialize as G1Projective (matches arkworks-spartan)
impl CanonicalSerialize for GroupElement {
    fn serialize_with_mode<W: std::io::Write>(
        &self,
        writer: W,
        compress: Compress,
    ) -> Result<(), ark_serialize::SerializationError> {
        self.0.serialize_with_mode(writer, compress)
    }

    fn serialized_size(&self, compress: Compress) -> usize {
        self.0.serialized_size(compress)
    }
}

impl Valid for GroupElement {
    fn check(&self) -> Result<(), ark_serialize::SerializationError> {
        self.0.check()
    }
}

impl CanonicalDeserialize for GroupElement {
    fn deserialize_with_mode<R: std::io::Read>(
        reader: R,
        compress: Compress,
        validate: Validate,
    ) -> Result<Self, ark_serialize::SerializationError> {
        G1Projective::deserialize_with_mode(reader, compress, validate).map(GroupElement)
    }
}

// CanonicalSerialize for CompressedGroup: serialize as raw bytes (G1Affine compressed)
impl CanonicalSerialize for CompressedGroup {
    fn serialize_with_mode<W: std::io::Write>(
        &self,
        mut writer: W,
        _compress: Compress,
    ) -> Result<(), ark_serialize::SerializationError> {
        // Write length prefix then bytes (to match Vec<u8> serialization pattern)
        writer.write_all(&self.0)?;
        Ok(())
    }

    fn serialized_size(&self, _compress: Compress) -> usize {
        self.0.len()
    }
}

impl Valid for CompressedGroup {
    fn check(&self) -> Result<(), ark_serialize::SerializationError> {
        // Validate by trying to decompress
        if self.decompress().is_some() {
            Ok(())
        } else {
            Err(ark_serialize::SerializationError::InvalidData)
        }
    }
}

impl CanonicalDeserialize for CompressedGroup {
    fn deserialize_with_mode<R: std::io::Read>(
        mut reader: R,
        _compress: Compress,
        _validate: Validate,
    ) -> Result<Self, ark_serialize::SerializationError> {
        // Read compressed G1Affine (32 bytes for BN254)
        let mut bytes = vec![0u8; 32];
        reader.read_exact(&mut bytes)?;
        Ok(CompressedGroup(bytes))
    }
}

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

    /// Get the underlying projective point
    pub fn inner(&self) -> &G1Projective {
        &self.0
    }

    /// Create from affine point
    pub fn from_affine(p: G1Affine) -> Self {
        GroupElement(p.into())
    }

    /// Variable-time MSM with pre-converted affine points (faster for repeated MSMs)
    pub fn msm_affine(scalars: &[Scalar], points: &[G1Affine]) -> Self {
        let scalars_vec: Vec<Fr> = scalars.iter().map(|s| s.0).collect();
        let result = G1Projective::msm(points, &scalars_vec).unwrap_or_default();
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

//! BN254 scalar field wrapper for Spartan
//! Replaces curve25519-dalek's Scalar with ark_bn254::Fr

use ark_bn254::Fr;
use ark_ff::{Field, PrimeField, UniformRand, Zero, One};
use ark_serialize::{CanonicalSerialize, CanonicalDeserialize};
use ark_std::rand::RngCore;
use std::ops::{Add, Sub, Mul, Neg, AddAssign, SubAssign, MulAssign};

/// Wrapper around ark_bn254::Fr for Spartan compatibility
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub struct Scalar(pub Fr);

impl Scalar {
    /// The zero scalar
    pub fn zero() -> Self {
        Scalar(Fr::zero())
    }

    /// The one scalar
    pub fn one() -> Self {
        Scalar(Fr::one())
    }

    /// Generate a random scalar
    pub fn random<R: RngCore>(rng: &mut R) -> Self {
        Scalar(Fr::rand(rng))
    }

    /// Compute the multiplicative inverse
    pub fn invert(&self) -> Option<Self> {
        self.0.inverse().map(Scalar)
    }

    /// Convert to bytes (little-endian, 32 bytes)
    pub fn to_bytes(&self) -> [u8; 32] {
        let mut bytes = [0u8; 32];
        // Serialize in little-endian format
        let bigint = self.0.into_bigint();
        for (i, limb) in bigint.0.iter().enumerate() {
            let limb_bytes = limb.to_le_bytes();
            bytes[i * 8..(i + 1) * 8].copy_from_slice(&limb_bytes);
        }
        bytes
    }

    /// Create from bytes (little-endian, 32 bytes)
    pub fn from_bytes(bytes: &[u8; 32]) -> Option<Self> {
        // Reconstruct from little-endian bytes
        let mut limbs = [0u64; 4];
        for i in 0..4 {
            limbs[i] = u64::from_le_bytes(bytes[i * 8..(i + 1) * 8].try_into().unwrap());
        }
        let bigint = ark_ff::BigInt(limbs);
        Fr::from_bigint(bigint).map(Scalar)
    }

    /// Create from u64
    pub fn from_u64(val: u64) -> Self {
        Scalar(Fr::from(val))
    }

    /// Square the scalar
    pub fn square(&self) -> Self {
        Scalar(self.0.square())
    }

    /// Raise to a power
    pub fn pow(&self, exp: &[u64]) -> Self {
        Scalar(self.0.pow(exp))
    }
}

// Arithmetic operations
impl Add for Scalar {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Scalar(self.0 + rhs.0)
    }
}

impl Sub for Scalar {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Scalar(self.0 - rhs.0)
    }
}

impl Mul for Scalar {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        Scalar(self.0 * rhs.0)
    }
}

impl Neg for Scalar {
    type Output = Self;
    fn neg(self) -> Self {
        Scalar(-self.0)
    }
}

impl AddAssign for Scalar {
    fn add_assign(&mut self, rhs: Self) {
        self.0 += rhs.0;
    }
}

impl SubAssign for Scalar {
    fn sub_assign(&mut self, rhs: Self) {
        self.0 -= rhs.0;
    }
}

impl MulAssign for Scalar {
    fn mul_assign(&mut self, rhs: Self) {
        self.0 *= rhs.0;
    }
}

impl std::iter::Sum for Scalar {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Scalar::zero(), |acc, x| acc + x)
    }
}

impl<'a> std::iter::Sum<&'a Scalar> for Scalar {
    fn sum<I: Iterator<Item = &'a Scalar>>(iter: I) -> Self {
        iter.fold(Scalar::zero(), |acc, x| acc + *x)
    }
}

impl std::iter::Product for Scalar {
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Scalar::one(), |acc, x| acc * x)
    }
}

// Serde support
impl serde::Serialize for Scalar {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        self.to_bytes().serialize(serializer)
    }
}

impl<'de> serde::Deserialize<'de> for Scalar {
    fn deserialize<D: serde::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        let bytes: [u8; 32] = serde::Deserialize::deserialize(deserializer)?;
        Scalar::from_bytes(&bytes).ok_or_else(|| serde::de::Error::custom("invalid scalar"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_scalar_arithmetic() {
        let a = Scalar::from_u64(5);
        let b = Scalar::from_u64(7);
        
        assert_eq!(a + b, Scalar::from_u64(12));
        assert_eq!(b - a, Scalar::from_u64(2));
        assert_eq!(a * b, Scalar::from_u64(35));
    }

    #[test]
    fn test_scalar_bytes_roundtrip() {
        let mut rng = ark_std::test_rng();
        for _ in 0..100 {
            let s = Scalar::random(&mut rng);
            let bytes = s.to_bytes();
            let recovered = Scalar::from_bytes(&bytes).unwrap();
            assert_eq!(s, recovered);
        }
    }

    #[test]
    fn test_scalar_invert() {
        let a = Scalar::from_u64(7);
        let a_inv = a.invert().unwrap();
        assert_eq!(a * a_inv, Scalar::one());
    }
}

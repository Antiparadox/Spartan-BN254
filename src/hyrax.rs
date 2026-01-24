//! Dense multilinear polynomial representation

use ark_serialize::{CanonicalSerialize, CanonicalDeserialize};
use crate::commitments::{Commitments, MultiCommitGens};
use crate::errors::ProofVerifyError;
use crate::group::GroupElement;
use crate::math::Math;
use crate::nizk::{DotProductProofGens, DotProductProofLog};
use crate::random::RandomTape;
use crate::scalar::Scalar;
use crate::transcript::{AppendToTranscript, ProofTranscript};
use core::ops::Index;
use merlin::Transcript;
use serde::{Deserialize, Serialize};

#[cfg(feature = "parallel")]
use rayon::prelude::*;

/// Polynomial commitment generators
#[derive(Clone, Serialize, Deserialize)]
pub struct PolyCommitmentGens {
    pub gens: DotProductProofGens,
}

impl PolyCommitmentGens {
    pub fn new(num_vars: usize, label: &'static [u8]) -> PolyCommitmentGens {
        let (_left, right) = EqPolynomial::compute_factored_lens(num_vars);
        let gens = DotProductProofGens::new(right.pow2(), label);
        PolyCommitmentGens { gens }
    }
}

/// Blinds for polynomial commitment
pub struct PolyCommitmentBlinds {
    pub blinds: Vec<Scalar>,
}

/// Polynomial commitment
#[derive(Debug, Serialize, Deserialize, Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct PolyCommitment {
    pub C: Vec<GroupElement>,
}

impl AppendToTranscript for PolyCommitment {
    fn append_to_transcript(&self, label: &'static [u8], transcript: &mut Transcript) {
        transcript.append_message(label, b"poly_commitment_begin");
        for c in &self.C {
            c.append_to_transcript(b"poly_commitment_share", transcript);
        }
        transcript.append_message(label, b"poly_commitment_end");
    }
}

/// Polynomial evaluation proof
#[derive(Debug, Clone, Serialize, Deserialize, CanonicalSerialize, CanonicalDeserialize)]
pub struct PolyEvalProof {
    proof: DotProductProofLog,
}

impl PolyEvalProof {
    fn protocol_name() -> &'static [u8] {
        b"polynomial evaluation proof"
    }

    pub fn prove(
        poly: &DensePolynomial,
        blinds_opt: Option<&PolyCommitmentBlinds>,
        r: &[Scalar],
        Zr: &Scalar,
        blind_Zr_opt: Option<&Scalar>,
        gens: &PolyCommitmentGens,
        transcript: &mut Transcript,
        random_tape: &mut RandomTape,
    ) -> (PolyEvalProof, GroupElement) {
        transcript.append_protocol_name(PolyEvalProof::protocol_name());

        assert_eq!(poly.get_num_vars(), r.len());

        let (left_num_vars, right_num_vars) = EqPolynomial::compute_factored_lens(r.len());
        let L_size = left_num_vars.pow2();
        let R_size = right_num_vars.pow2();

        let default_blinds = PolyCommitmentBlinds {
            blinds: vec![Scalar::zero(); L_size],
        };
        let blinds = blinds_opt.map_or(&default_blinds, |p| p);

        assert_eq!(blinds.blinds.len(), L_size);

        let zero = Scalar::zero();
        let blind_Zr = blind_Zr_opt.map_or(&zero, |p| p);

        // compute the L and R vectors
        let eq = EqPolynomial::new(r.to_vec());
        let (L, R) = eq.compute_factored_evals();
        assert_eq!(L.len(), L_size);
        assert_eq!(R.len(), R_size);

        // compute the vector underneath L*Z and the L*blinds
        let LZ = poly.bound(&L);
        let LZ_blind: Scalar = (0..L.len()).map(|i| blinds.blinds[i] * L[i]).sum();

        // a dot product proof of size R_size
        let (proof, _C_LR, C_Zr_prime) = DotProductProofLog::prove(
            &gens.gens,
            transcript,
            random_tape,
            &LZ,
            &LZ_blind,
            &R,
            Zr,
            blind_Zr,
        );

        (PolyEvalProof { proof }, C_Zr_prime)
    }

    pub fn verify(
        &self,
        gens: &PolyCommitmentGens,
        transcript: &mut Transcript,
        r: &[Scalar],
        C_Zr: &GroupElement,
        comm: &PolyCommitment,
    ) -> Result<(), ProofVerifyError> {
        transcript.append_protocol_name(PolyEvalProof::protocol_name());

        // compute L and R
        let eq = EqPolynomial::new(r.to_vec());
        let (L, R) = eq.compute_factored_evals();

        // compute a weighted sum of commitments and L
        let C_LZ = GroupElement::vartime_multiscalar_mul(&L, &comm.C);

        self.proof
            .verify(R.len(), &gens.gens, transcript, &R, &C_LZ, C_Zr)
    }

    pub fn verify_plain(
        &self,
        gens: &PolyCommitmentGens,
        transcript: &mut Transcript,
        r: &[Scalar],
        Zr: &Scalar,
        comm: &PolyCommitment,
    ) -> Result<(), ProofVerifyError> {
        // compute a commitment to Zr with a blind of zero
        let C_Zr = Zr.commit(&Scalar::zero(), &gens.gens.gens_1);

        self.verify(gens, transcript, r, &C_Zr, comm)
    }
}

/// Dense multilinear polynomial
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DensePolynomial {
    num_vars: usize,
    len: usize,
    Z: Vec<Scalar>,
}

impl DensePolynomial {
    pub fn new(Z: Vec<Scalar>) -> Self {
        let len = Z.len();
        let num_vars = if len > 0 { len.log_2() } else { 0 };
        DensePolynomial { num_vars, len, Z }
    }

    pub fn get_num_vars(&self) -> usize {
        self.num_vars
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Return the evaluation vector
    pub fn evals(&self) -> &[Scalar] {
        &self.Z
    }

    pub fn split(&self, idx: usize) -> (DensePolynomial, DensePolynomial) {
        assert!(idx < self.len());
        (
            DensePolynomial::new(self.Z[..idx].to_vec()),
            DensePolynomial::new(self.Z[idx..2 * idx].to_vec()),
        )
    }

    /// Bind the polynomial's top variable to r
    pub fn bound_poly_var_top(&mut self, r: &Scalar) {
        let n = self.len() / 2;
        for i in 0..n {
            self.Z[i] = self.Z[i] + *r * (self.Z[i + n] - self.Z[i]);
        }
        self.Z.truncate(n);
        self.num_vars -= 1;
        self.len = n;
    }

    /// Bind the polynomial's bottom variable to r
    pub fn bound_poly_var_bot(&mut self, r: &Scalar) {
        let n = self.len() / 2;
        for i in 0..n {
            self.Z[i] = self.Z[2 * i] + *r * (self.Z[2 * i + 1] - self.Z[2 * i]);
        }
        self.Z.truncate(n);
        self.num_vars -= 1;
        self.len = n;
    }

    /// Evaluate the polynomial at point r
    pub fn evaluate(&self, r: &[Scalar]) -> Scalar {
        assert_eq!(r.len(), self.get_num_vars());
        let chis = EqPolynomial::new(r.to_vec()).evals();
        assert_eq!(chis.len(), self.Z.len());
        compute_dotproduct(&self.Z, &chis)
    }

    pub fn vec(&self) -> &Vec<Scalar> {
        &self.Z
    }

    pub fn extend(&mut self, other: &DensePolynomial) {
        assert_eq!(self.Z.len(), self.len);
        let other_vec = other.vec();
        assert_eq!(other_vec.len(), self.len);
        self.Z.extend(other_vec);
        self.num_vars += 1;
        self.len *= 2;
    }

    pub fn merge<'a, I>(polys: I) -> DensePolynomial
    where
        I: IntoIterator<Item = &'a DensePolynomial>,
    {
        let mut Z: Vec<Scalar> = Vec::new();
        for poly in polys.into_iter() {
            Z.extend(poly.vec());
        }
        Z.resize(Z.len().next_power_of_two(), Scalar::zero());
        DensePolynomial::new(Z)
    }

    pub fn from_usize(Z: &[usize]) -> Self {
        DensePolynomial::new(Z.iter().map(|&x| Scalar::from_u64(x as u64)).collect())
    }

    #[cfg(feature = "parallel")]
    fn commit_inner(&self, blinds: &[Scalar], gens: &MultiCommitGens) -> PolyCommitment {
        use rayon::prelude::*;
        let L_size = blinds.len();
        let R_size = self.Z.len() / L_size;
        assert_eq!(L_size * R_size, self.Z.len());
        let C = (0..L_size)
            .into_par_iter()
            .map(|i| {
                self.Z[R_size * i..R_size * (i + 1)]
                    .commit(&blinds[i], gens)
            })
            .collect();
        PolyCommitment { C }
    }

    #[cfg(not(feature = "parallel"))]
    fn commit_inner(&self, blinds: &[Scalar], gens: &MultiCommitGens) -> PolyCommitment {
        let L_size = blinds.len();
        let R_size = self.Z.len() / L_size;
        assert_eq!(L_size * R_size, self.Z.len());
        let C = (0..L_size)
            .map(|i| {
                self.Z[R_size * i..R_size * (i + 1)]
                    .commit(&blinds[i], gens)
            })
            .collect();
        PolyCommitment { C }
    }

    pub fn commit(
        &self,
        gens: &PolyCommitmentGens,
        random_tape: Option<&mut RandomTape>,
    ) -> (PolyCommitment, PolyCommitmentBlinds) {
        let n = self.Z.len();
        let ell = self.get_num_vars();
        assert_eq!(n, ell.pow2());

        let (left_num_vars, right_num_vars) = EqPolynomial::compute_factored_lens(ell);
        let L_size = left_num_vars.pow2();
        let R_size = right_num_vars.pow2();
        assert_eq!(L_size * R_size, n);

        let blinds = if let Some(t) = random_tape {
            PolyCommitmentBlinds {
                blinds: t.random_vector(b"poly_blinds", L_size),
            }
        } else {
            PolyCommitmentBlinds {
                blinds: vec![Scalar::zero(); L_size],
            }
        };

        (self.commit_inner(&blinds.blinds, &gens.gens.gens_n), blinds)
    }

    /// Compute L * Z where Z is viewed as a matrix
    pub fn bound(&self, L: &[Scalar]) -> Vec<Scalar> {
        let (left_num_vars, right_num_vars) = EqPolynomial::compute_factored_lens(self.get_num_vars());
        let L_size = left_num_vars.pow2();
        let R_size = right_num_vars.pow2();

        #[cfg(feature = "parallel")]
        let iter = (0..R_size).into_par_iter();
        #[cfg(not(feature = "parallel"))]
        let iter = 0..R_size;

        iter.map(|i| {
            (0..L_size).map(|j| L[j] * self.Z[j * R_size + i]).sum()
        }).collect()
    }
}

impl Index<usize> for DensePolynomial {
    type Output = Scalar;

    #[inline(always)]
    fn index(&self, index: usize) -> &Scalar {
        &self.Z[index]
    }
}

/// Equality polynomial: eq(x, r) = prod_i (r_i * x_i + (1 - r_i) * (1 - x_i))
pub struct EqPolynomial {
    r: Vec<Scalar>,
}

impl EqPolynomial {
    pub fn new(r: Vec<Scalar>) -> Self {
        EqPolynomial { r }
    }

    /// Evaluate eq(rx, ry)
    pub fn evaluate(&self, rx: &[Scalar]) -> Scalar {
        assert_eq!(self.r.len(), rx.len());
        (0..rx.len())
            .map(|i| self.r[i] * rx[i] + (Scalar::one() - self.r[i]) * (Scalar::one() - rx[i]))
            .product()
    }

    /// Compute all evaluations of eq(r, x) for x in {0,1}^n
    pub fn evals(&self) -> Vec<Scalar> {
        let ell = self.r.len();
        let mut evals: Vec<Scalar> = vec![Scalar::one(); ell.pow2()];
        let mut size = 1;

        for j in 0..ell {
            size *= 2;
            for i in (0..size).rev().step_by(2) {
                let scalar = evals[i / 2];
                evals[i] = scalar * self.r[j];
                evals[i - 1] = scalar - evals[i];
            }
        }
        evals
    }

    pub fn compute_factored_lens(ell: usize) -> (usize, usize) {
        (ell / 2, ell - ell / 2)
    }

    pub fn compute_factored_evals(&self) -> (Vec<Scalar>, Vec<Scalar>) {
        let ell = self.r.len();
        let (left_num_vars, _right_num_vars) = Self::compute_factored_lens(ell);

        let L = EqPolynomial::new(self.r[..left_num_vars].to_vec()).evals();
        let R = EqPolynomial::new(self.r[left_num_vars..ell].to_vec()).evals();

        (L, R)
    }
}

/// Identity polynomial
pub struct IdentityPolynomial {
    size_point: usize,
}

impl IdentityPolynomial {
    pub fn new(size_point: usize) -> Self {
        IdentityPolynomial { size_point }
    }

    pub fn evaluate(&self, r: &[Scalar]) -> Scalar {
        let len = r.len();
        assert_eq!(len, self.size_point);
        (0..len)
            .map(|i| Scalar::from_u64((len - i - 1).pow2() as u64) * r[i])
            .sum()
    }
}

/// Compute dot product of two vectors
pub fn compute_dotproduct(a: &[Scalar], b: &[Scalar]) -> Scalar {
    assert_eq!(a.len(), b.len());

    #[cfg(feature = "parallel")]
    {
        a.par_iter()
            .zip(b.par_iter())
            .map(|(a_i, b_i)| *a_i * *b_i)
            .sum()
    }

    #[cfg(not(feature = "parallel"))]
    {
        a.iter().zip(b.iter()).map(|(a_i, b_i)| *a_i * *b_i).sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_eq_polynomial() {
        let r = vec![Scalar::from_u64(2), Scalar::from_u64(3)];
        let eq = EqPolynomial::new(r);
        let evals = eq.evals();
        assert_eq!(evals.len(), 4);
    }

    #[test]
    fn test_dense_polynomial_evaluate() {
        // Z = [1, 2, 3, 4] represents f(x0, x1) where f(0,0)=1, f(0,1)=2, f(1,0)=3, f(1,1)=4
        let Z = vec![
            Scalar::from_u64(1),
            Scalar::from_u64(2),
            Scalar::from_u64(3),
            Scalar::from_u64(4),
        ];
        let poly = DensePolynomial::new(Z);
        
        // Evaluate at (0, 0) should give 1
        let r = vec![Scalar::zero(), Scalar::zero()];
        assert_eq!(poly.evaluate(&r), Scalar::from_u64(1));
        
        // Evaluate at (1, 1) should give 4
        let r = vec![Scalar::one(), Scalar::one()];
        assert_eq!(poly.evaluate(&r), Scalar::from_u64(4));
    }

    #[test]
    fn test_bound_poly_var_top() {
        let Z = vec![
            Scalar::from_u64(1),
            Scalar::from_u64(2),
            Scalar::from_u64(3),
            Scalar::from_u64(4),
        ];
        let mut poly = DensePolynomial::new(Z);
        
        // Bind first variable to 0
        poly.bound_poly_var_top(&Scalar::zero());
        assert_eq!(poly.len(), 2);
        assert_eq!(poly[0], Scalar::from_u64(1));
        assert_eq!(poly[1], Scalar::from_u64(2));
    }
}

//! KZG Polynomial Commitment Scheme
//!
//! This module provides KZG commitments as an alternative to Hyrax,
//! offering O(1) proof size instead of O(√n).
//!
//! Key differences from Hyrax:
//! - Commitment: 1 G1 point (vs √n points)
//! - Opening proof: 1 G1 point (vs log(n) points + scalars)
//! - Requires trusted setup (powers of tau)
//! - Uses pairings for verification

use ark_bn254::{Bn254, Fr, G1Affine, G1Projective, G2Affine, G2Projective};
use ark_ec::{pairing::Pairing, CurveGroup, PrimeGroup, VariableBaseMSM};
use ark_ff::{UniformRand, One};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::rand::RngCore;
use merlin::Transcript;

use crate::scalar::Scalar;
use crate::transcript::AppendToTranscript;

/// KZG Structured Reference String (SRS)
/// Contains powers of tau: [τ^0]G1, [τ^1]G1, ..., [τ^n]G1, [τ]G2
#[derive(Clone, Debug, CanonicalSerialize, CanonicalDeserialize)]
pub struct KZGSrs {
    /// Powers of tau in G1: [τ^i]G1 for i = 0..max_degree
    pub powers_g1: Vec<G1Affine>,
    /// [τ]G2 for pairing verification
    pub tau_g2: G2Affine,
    /// [1]G2 (generator)
    pub g2: G2Affine,
}

impl KZGSrs {
    /// Generate a new SRS with random tau (FOR TESTING ONLY - not secure!)
    /// In production, use a ceremony like powers of tau.
    pub fn setup<R: RngCore>(max_degree: usize, rng: &mut R) -> Self {
        let tau = Fr::rand(rng);
        let g1 = G1Projective::generator();
        let g2 = G2Projective::generator();
        
        // Compute powers of tau in G1
        let mut powers_g1 = Vec::with_capacity(max_degree + 1);
        let mut tau_power = Fr::one();
        for _ in 0..=max_degree {
            powers_g1.push((g1 * tau_power).into_affine());
            tau_power *= tau;
        }
        
        KZGSrs {
            powers_g1,
            tau_g2: (g2 * tau).into_affine(),
            g2: g2.into_affine(),
        }
    }
    
    /// Generate a deterministic SRS from a seed (FOR TESTING)
    pub fn setup_from_seed(max_degree: usize, seed: u64) -> Self {
        use ark_std::rand::SeedableRng;
        use ark_std::rand::rngs::StdRng;
        let mut rng = StdRng::seed_from_u64(seed);
        Self::setup(max_degree, &mut rng)
    }
    
    /// Save SRS to file
    pub fn save_to_file(&self, path: &str) -> std::io::Result<()> {
        use std::fs::File;
        use std::io::Write;
        
        let mut bytes = Vec::new();
        self.serialize_compressed(&mut bytes)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e.to_string()))?;
        
        let mut file = File::create(path)?;
        file.write_all(&bytes)?;
        Ok(())
    }
    
    /// Load SRS from file
    pub fn load_from_file(path: &str) -> std::io::Result<Self> {
        use std::fs::File;
        use std::io::Read;
        
        let mut file = File::open(path)?;
        let mut bytes = Vec::new();
        file.read_to_end(&mut bytes)?;
        
        Self::deserialize_compressed(&bytes[..])
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e.to_string()))
    }
    
    /// Load SRS from file, or generate and save if not found
    pub fn load_or_generate(path: &str, max_degree: usize, seed: u64) -> std::io::Result<Self> {
        match Self::load_from_file(path) {
            Ok(srs) => {
                // Check if loaded SRS is large enough
                if srs.powers_g1.len() > max_degree {
                    Ok(srs)
                } else {
                    // Need to regenerate with larger size
                    println!("  Loaded SRS too small ({} < {}), regenerating...", 
                             srs.powers_g1.len(), max_degree);
                    let new_srs = Self::setup_from_seed(max_degree, seed);
                    new_srs.save_to_file(path)?;
                    Ok(new_srs)
                }
            }
            Err(_) => {
                // Generate and save
                let srs = Self::setup_from_seed(max_degree, seed);
                srs.save_to_file(path)?;
                Ok(srs)
            }
        }
    }
    
    /// Get the maximum degree this SRS can support
    pub fn max_degree(&self) -> usize {
        self.powers_g1.len().saturating_sub(1)
    }
}

/// KZG Polynomial Commitment (just one G1 point!)
#[derive(Clone, Debug, CanonicalSerialize, CanonicalDeserialize)]
pub struct KZGCommitment {
    pub commitment: G1Affine,
}

impl KZGCommitment {
    /// Commit to a polynomial using KZG
    /// C = Σ coeffs[i] * [τ^i]G1
    pub fn commit(coeffs: &[Scalar], srs: &KZGSrs) -> Self {
        assert!(coeffs.len() <= srs.powers_g1.len(), 
                "Polynomial degree {} exceeds SRS size {}", 
                coeffs.len(), srs.powers_g1.len());
        
        // Convert Scalars to Fr
        let fr_coeffs: Vec<Fr> = coeffs.iter().map(|s| s.inner()).collect();
        
        // MSM: C = Σ coeffs[i] * powers_g1[i]
        let commitment = G1Projective::msm(&srs.powers_g1[..coeffs.len()], &fr_coeffs)
            .expect("MSM failed")
            .into_affine();
        
        KZGCommitment { commitment }
    }
    
    /// Commit to evaluations of a multilinear polynomial (converts to univariate first)
    pub fn commit_multilinear(evals: &[Scalar], srs: &KZGSrs) -> Self {
        // For multilinear polynomials, we can commit directly to evaluations
        // using Lagrange basis, or convert to coefficient form.
        // For simplicity, we commit to evaluations directly as coefficients.
        Self::commit(evals, srs)
    }
}

impl AppendToTranscript for KZGCommitment {
    fn append_to_transcript(&self, label: &'static [u8], transcript: &mut Transcript) {
        let mut bytes = Vec::new();
        self.commitment.serialize_compressed(&mut bytes).unwrap();
        transcript.append_message(label, &bytes);
    }
}

/// KZG Opening Proof (just one G1 point!)
#[derive(Clone, Debug, CanonicalSerialize, CanonicalDeserialize)]
pub struct KZGProof {
    pub proof: G1Affine,
}

impl KZGProof {
    /// Create an opening proof for p(z) = y
    /// Proof: π = [q(τ)]G1 where q(x) = (p(x) - y) / (x - z)
    pub fn prove(coeffs: &[Scalar], point: &Scalar, srs: &KZGSrs) -> (Self, Scalar) {
        // Evaluate p(z)
        let eval = Self::evaluate_poly(coeffs, point);
        
        // Compute quotient polynomial q(x) = (p(x) - y) / (x - z)
        let quotient_coeffs = Self::compute_quotient(coeffs, point, &eval);
        
        // Commit to quotient
        let fr_coeffs: Vec<Fr> = quotient_coeffs.iter().map(|s| s.inner()).collect();
        let proof = if fr_coeffs.is_empty() {
            G1Affine::identity()
        } else {
            G1Projective::msm(&srs.powers_g1[..fr_coeffs.len()], &fr_coeffs)
                .expect("MSM failed")
                .into_affine()
        };
        
        (KZGProof { proof }, eval)
    }
    
    /// Verify an opening proof
    /// Check: e(C - [y]G1, G2) = e(π, [τ]G2 - [z]G2)
    pub fn verify(
        &self,
        commitment: &KZGCommitment,
        point: &Scalar,
        eval: &Scalar,
        srs: &KZGSrs,
    ) -> bool {
        let g1 = G1Projective::generator();
        let g2 = G2Projective::generator();
        
        // C - [y]G1
        let lhs_g1: G1Affine = (G1Projective::from(commitment.commitment) - g1 * eval.inner()).into_affine();
        
        // [τ]G2 - [z]G2
        let rhs_g2: G2Affine = (G2Projective::from(srs.tau_g2) - g2 * point.inner()).into_affine();
        
        // Check pairing: e(C - [y]G1, G2) = e(π, [τ - z]G2)
        let lhs = Bn254::pairing(lhs_g1, srs.g2);
        let rhs = Bn254::pairing(self.proof, rhs_g2);
        
        lhs == rhs
    }
    
    /// Evaluate polynomial at a point
    fn evaluate_poly(coeffs: &[Scalar], point: &Scalar) -> Scalar {
        let mut result = Scalar::zero();
        let mut power = Scalar::one();
        for coeff in coeffs {
            result = result + *coeff * power;
            power = power * *point;
        }
        result
    }
    
    /// Compute quotient polynomial (p(x) - y) / (x - z)
    fn compute_quotient(coeffs: &[Scalar], z: &Scalar, y: &Scalar) -> Vec<Scalar> {
        if coeffs.is_empty() {
            return vec![];
        }
        
        // p(x) - y as coefficients
        let mut p_minus_y: Vec<Scalar> = coeffs.to_vec();
        p_minus_y[0] = p_minus_y[0] - *y;
        
        // Synthetic division by (x - z)
        let n = p_minus_y.len();
        if n <= 1 {
            return vec![];
        }
        
        let mut quotient = vec![Scalar::zero(); n - 1];
        let mut remainder = p_minus_y[n - 1];
        
        for i in (0..n - 1).rev() {
            quotient[i] = remainder;
            remainder = p_minus_y[i] + remainder * *z;
        }
        
        // remainder should be 0 if y = p(z)
        quotient
    }
}

/// Batch KZG opening for multiple polynomials at the same point
#[derive(Clone, Debug, CanonicalSerialize, CanonicalDeserialize)]
pub struct KZGBatchProof {
    pub proof: G1Affine,
}

impl KZGBatchProof {
    /// Batch prove: open multiple polynomials at the same point
    /// Uses random linear combination
    pub fn batch_prove(
        polys: &[Vec<Scalar>],
        point: &Scalar,
        evals: &[Scalar],
        srs: &KZGSrs,
        transcript: &mut Transcript,
    ) -> Self {
        // Get random challenge for linear combination
        let mut challenge_bytes = [0u8; 32];
        transcript.challenge_bytes(b"batch_challenge", &mut challenge_bytes);
        let gamma = Scalar::from_bytes(&challenge_bytes).unwrap_or(Scalar::one());
        
        // Compute linear combination of polynomials
        let mut combined_coeffs = vec![Scalar::zero(); polys.iter().map(|p| p.len()).max().unwrap_or(0)];
        let mut gamma_power = Scalar::one();
        
        for poly in polys {
            for (i, coeff) in poly.iter().enumerate() {
                combined_coeffs[i] = combined_coeffs[i] + *coeff * gamma_power;
            }
            gamma_power = gamma_power * gamma;
        }
        
        // Compute combined evaluation
        let mut combined_eval = Scalar::zero();
        gamma_power = Scalar::one();
        for eval in evals {
            combined_eval = combined_eval + *eval * gamma_power;
            gamma_power = gamma_power * gamma;
        }
        
        // Prove the combined polynomial
        let quotient_coeffs = KZGProof::compute_quotient(&combined_coeffs, point, &combined_eval);
        
        let fr_coeffs: Vec<Fr> = quotient_coeffs.iter().map(|s| s.inner()).collect();
        let proof = if fr_coeffs.is_empty() {
            G1Affine::identity()
        } else {
            G1Projective::msm(&srs.powers_g1[..fr_coeffs.len()], &fr_coeffs)
                .expect("MSM failed")
                .into_affine()
        };
        
        KZGBatchProof { proof }
    }
    
    /// Verify batch opening
    pub fn batch_verify(
        &self,
        commitments: &[KZGCommitment],
        point: &Scalar,
        evals: &[Scalar],
        srs: &KZGSrs,
        transcript: &mut Transcript,
    ) -> bool {
        // Get same random challenge
        let mut challenge_bytes = [0u8; 32];
        transcript.challenge_bytes(b"batch_challenge", &mut challenge_bytes);
        let gamma = Scalar::from_bytes(&challenge_bytes).unwrap_or(Scalar::one());
        
        // Compute linear combination of commitments
        let mut combined_comm = G1Projective::from(G1Affine::identity());
        let mut gamma_power = Scalar::one();
        
        for comm in commitments {
            combined_comm = combined_comm + G1Projective::from(comm.commitment) * gamma_power.inner();
            gamma_power = gamma_power * gamma;
        }
        
        // Compute combined evaluation
        let mut combined_eval = Scalar::zero();
        gamma_power = Scalar::one();
        for eval in evals {
            combined_eval = combined_eval + *eval * gamma_power;
            gamma_power = gamma_power * gamma;
        }
        
        // Verify the combined proof
        let combined_commitment = KZGCommitment {
            commitment: combined_comm.into_affine(),
        };
        
        let single_proof = KZGProof { proof: self.proof };
        single_proof.verify(&combined_commitment, point, &combined_eval, srs)
    }
}

// ============================================================================
// KZG Polynomial Commitment for Dense Polynomials (replacement for Hyrax)
// ============================================================================

/// KZG-based polynomial commitment generators (replaces PolyCommitmentGens)
#[derive(Clone)]
pub struct KZGPolyCommitmentGens {
    pub srs: KZGSrs,
}

impl KZGPolyCommitmentGens {
    /// Create new generators from SRS
    pub fn new(srs: KZGSrs) -> Self {
        KZGPolyCommitmentGens { srs }
    }
    
    /// Create from file or generate
    pub fn from_file_or_generate(path: &str, max_degree: usize, seed: u64) -> std::io::Result<Self> {
        let srs = KZGSrs::load_or_generate(path, max_degree, seed)?;
        Ok(Self::new(srs))
    }
}

/// KZG polynomial commitment (single G1 point, replaces PolyCommitment)
#[derive(Debug, Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct KZGPolyCommitment {
    pub commitment: G1Affine,
}

impl KZGPolyCommitment {
    /// Commit to polynomial evaluations
    pub fn commit(evals: &[Scalar], gens: &KZGPolyCommitmentGens) -> Self {
        let fr_evals: Vec<Fr> = evals.iter().map(|s| s.inner()).collect();
        let n = fr_evals.len().min(gens.srs.powers_g1.len());
        
        let commitment = G1Projective::msm(&gens.srs.powers_g1[..n], &fr_evals[..n])
            .expect("MSM failed")
            .into_affine();
        
        KZGPolyCommitment { commitment }
    }
}

impl AppendToTranscript for KZGPolyCommitment {
    fn append_to_transcript(&self, label: &'static [u8], transcript: &mut Transcript) {
        let mut bytes = Vec::new();
        self.commitment.serialize_compressed(&mut bytes).unwrap();
        transcript.append_message(label, &bytes);
    }
}

/// KZG polynomial evaluation proof (single G1 point, replaces PolyEvalProof)
#[derive(Debug, Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct KZGPolyEvalProof {
    pub proof: G1Affine,
    pub eval: Scalar,
}

impl KZGPolyEvalProof {
    /// Prove evaluation of polynomial at a point
    pub fn prove(
        evals: &[Scalar],
        point: &Scalar,
        gens: &KZGPolyCommitmentGens,
    ) -> Self {
        let (kzg_proof, eval) = KZGProof::prove(evals, point, &gens.srs);
        KZGPolyEvalProof {
            proof: kzg_proof.proof,
            eval,
        }
    }
    
    /// Verify evaluation proof
    pub fn verify(
        &self,
        comm: &KZGPolyCommitment,
        point: &Scalar,
        gens: &KZGPolyCommitmentGens,
    ) -> bool {
        let kzg_proof = KZGProof { proof: self.proof };
        let kzg_comm = KZGCommitment { commitment: comm.commitment };
        kzg_proof.verify(&kzg_comm, point, &self.eval, &gens.srs)
    }
}

/// Batched KZG commitment for multiple polynomials
#[derive(Debug, Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct KZGBatchedCommitment {
    pub commitments: Vec<G1Affine>,
}

impl KZGBatchedCommitment {
    /// Commit to multiple polynomials
    pub fn commit(polys: &[&[Scalar]], gens: &KZGPolyCommitmentGens) -> Self {
        let commitments: Vec<G1Affine> = polys
            .iter()
            .map(|poly| KZGPolyCommitment::commit(poly, gens).commitment)
            .collect();
        KZGBatchedCommitment { commitments }
    }
}

impl AppendToTranscript for KZGBatchedCommitment {
    fn append_to_transcript(&self, label: &'static [u8], transcript: &mut Transcript) {
        transcript.append_message(label, b"kzg_batch_begin");
        for (i, c) in self.commitments.iter().enumerate() {
            let mut bytes = Vec::new();
            c.serialize_compressed(&mut bytes).unwrap();
            transcript.append_message(b"kzg_batch_elem", &bytes);
        }
        transcript.append_message(label, b"kzg_batch_end");
    }
}

/// Batched KZG evaluation proof for multiple polynomials at same point
#[derive(Debug, Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct KZGBatchedEvalProof {
    pub proof: G1Affine,
    pub evals: Vec<Scalar>,
}

impl KZGBatchedEvalProof {
    /// Prove evaluation of multiple polynomials at same point
    pub fn prove(
        polys: &[&[Scalar]],
        point: &Scalar,
        gens: &KZGPolyCommitmentGens,
        transcript: &mut Transcript,
    ) -> Self {
        // Evaluate each polynomial at the point
        let evals: Vec<Scalar> = polys
            .iter()
            .map(|poly| KZGProof::evaluate_poly(poly, point))
            .collect();
        
        // Convert to owned vecs for batch prove
        let owned_polys: Vec<Vec<Scalar>> = polys.iter().map(|p| p.to_vec()).collect();
        
        // Batch prove
        let batch_proof = KZGBatchProof::batch_prove(&owned_polys, point, &evals, &gens.srs, transcript);
        
        KZGBatchedEvalProof {
            proof: batch_proof.proof,
            evals,
        }
    }
    
    /// Verify batched evaluation proof
    pub fn verify(
        &self,
        comm: &KZGBatchedCommitment,
        point: &Scalar,
        gens: &KZGPolyCommitmentGens,
        transcript: &mut Transcript,
    ) -> bool {
        let kzg_comms: Vec<KZGCommitment> = comm.commitments
            .iter()
            .map(|c| KZGCommitment { commitment: *c })
            .collect();
        
        let batch_proof = KZGBatchProof { proof: self.proof };
        batch_proof.batch_verify(&kzg_comms, point, &self.evals, &gens.srs, transcript)
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_kzg_commit_open() {
        // Setup SRS
        let srs = KZGSrs::setup_from_seed(1024, 12345);
        
        // Create a polynomial: p(x) = 1 + 2x + 3x²
        let coeffs = vec![
            Scalar::from_u64(1),
            Scalar::from_u64(2),
            Scalar::from_u64(3),
        ];
        
        // Commit
        let commitment = KZGCommitment::commit(&coeffs, &srs);
        
        // Open at point z = 5
        let z = Scalar::from_u64(5);
        let (proof, eval) = KZGProof::prove(&coeffs, &z, &srs);
        
        // Expected: p(5) = 1 + 2*5 + 3*25 = 1 + 10 + 75 = 86
        assert_eq!(eval, Scalar::from_u64(86));
        
        // Verify
        assert!(proof.verify(&commitment, &z, &eval, &srs));
        
        // Verify fails with wrong eval
        let wrong_eval = Scalar::from_u64(100);
        assert!(!proof.verify(&commitment, &z, &wrong_eval, &srs));
    }
    
    #[test]
    fn test_kzg_proof_size() {
        let srs = KZGSrs::setup_from_seed(1024, 12345);
        
        // Large polynomial
        let coeffs: Vec<Scalar> = (0..1000).map(|i| Scalar::from_u64(i)).collect();
        let commitment = KZGCommitment::commit(&coeffs, &srs);
        let z = Scalar::from_u64(42);
        let (proof, _) = KZGProof::prove(&coeffs, &z, &srs);
        
        // Measure sizes
        let mut comm_bytes = Vec::new();
        commitment.serialize_compressed(&mut comm_bytes).unwrap();
        
        let mut proof_bytes = Vec::new();
        proof.serialize_compressed(&mut proof_bytes).unwrap();
        
        println!("Commitment size: {} bytes", comm_bytes.len());
        println!("Proof size: {} bytes", proof_bytes.len());
        
        // Both should be 32 bytes (1 compressed G1 point)
        assert_eq!(comm_bytes.len(), 32);
        assert_eq!(proof_bytes.len(), 32);
    }
}

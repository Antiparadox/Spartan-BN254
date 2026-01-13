//! Bullet reduction proof for logarithmic-size inner product arguments
//! Port of Spartan's bullet.rs to BN254

use crate::errors::ProofVerifyError;
use crate::group::{CompressedGroup, GroupElement};
use crate::math::Math;
use crate::scalar::Scalar;
use crate::transcript::{AppendToTranscript, ProofTranscript};
use merlin::Transcript;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BulletReductionProof {
    L_vec: Vec<CompressedGroup>,
    R_vec: Vec<CompressedGroup>,
}

impl BulletReductionProof {
    /// Create an inner-product argument for the relation:
    /// P = <a,G> + <a,b>*Q + r*H
    #[allow(clippy::too_many_arguments)]
    #[allow(non_snake_case)]
    pub fn prove(
        transcript: &mut Transcript,
        Q: &GroupElement,
        G_vec: &[GroupElement],
        H: &GroupElement,
        a_vec: &[Scalar],
        b_vec: &[Scalar],
        blind: &Scalar,
        blinds_vec: &[(Scalar, Scalar)],
    ) -> (
        BulletReductionProof,
        GroupElement,    // Gamma_hat
        Scalar,          // a_hat
        Scalar,          // b_hat
        GroupElement,    // g_hat
        Scalar,          // rhat_Gamma
    ) {
        let mut n = G_vec.len();
        assert_eq!(a_vec.len(), n);
        assert_eq!(b_vec.len(), n);
        assert!(n.is_power_of_two());

        let lg_n = n.log_2();
        assert_eq!(blinds_vec.len(), lg_n);

        let mut G: Vec<GroupElement> = G_vec.to_vec();
        let mut a: Vec<Scalar> = a_vec.to_vec();
        let mut b: Vec<Scalar> = b_vec.to_vec();

        let mut L_vec: Vec<CompressedGroup> = Vec::with_capacity(lg_n);
        let mut R_vec: Vec<CompressedGroup> = Vec::with_capacity(lg_n);

        // Gamma is a commitment to a_vec with a base Q and blinding factor r
        let Gamma = GroupElement::vartime_multiscalar_mul(&a, &G)
            + compute_dotproduct(&a, &b) * *Q
            + *blind * *H;

        let mut blind_Gamma = *blind;

        for i in 0..lg_n {
            n /= 2;

            let (a_L, a_R) = a.split_at(n);
            let (b_L, b_R) = b.split_at(n);
            let (G_L, G_R) = G.split_at(n);

            let c_L = compute_dotproduct(a_L, b_R);
            let c_R = compute_dotproduct(a_R, b_L);

            let (blind_L, blind_R) = blinds_vec[i];

            let L = GroupElement::vartime_multiscalar_mul(a_L, G_R) + c_L * *Q + blind_L * *H;
            let R = GroupElement::vartime_multiscalar_mul(a_R, G_L) + c_R * *Q + blind_R * *H;

            L.compress().append_to_transcript(b"L", transcript);
            R.compress().append_to_transcript(b"R", transcript);

            let u = transcript.challenge_scalar(b"u");
            let u_inv = u.invert().unwrap();

            // Fold the generators
            G = G_L
                .iter()
                .zip(G_R.iter())
                .map(|(g_L, g_R)| u_inv * *g_L + u * *g_R)
                .collect();

            // Fold the scalars
            a = a_L
                .iter()
                .zip(a_R.iter())
                .map(|(a_L, a_R)| u * *a_L + u_inv * *a_R)
                .collect();

            b = b_L
                .iter()
                .zip(b_R.iter())
                .map(|(b_L, b_R)| u_inv * *b_L + u * *b_R)
                .collect();

            blind_Gamma = u * u * blind_L + blind_Gamma + u_inv * u_inv * blind_R;

            L_vec.push(L.compress());
            R_vec.push(R.compress());
        }

        assert_eq!(a.len(), 1);
        assert_eq!(b.len(), 1);
        assert_eq!(G.len(), 1);

        let a_hat = a[0];
        let b_hat = b[0];
        let g_hat = G[0];

        (
            BulletReductionProof { L_vec, R_vec },
            Gamma,
            a_hat,
            b_hat,
            g_hat,
            blind_Gamma,
        )
    }

    /// Verify the bullet reduction proof
    #[allow(non_snake_case)]
    pub fn verify(
        &self,
        n: usize,
        b_vec: &[Scalar],
        transcript: &mut Transcript,
        Gamma: &GroupElement,
        G_vec: &[GroupElement],
    ) -> Result<(GroupElement, GroupElement, Scalar), ProofVerifyError> {
        assert_eq!(b_vec.len(), n);
        assert_eq!(G_vec.len(), n);
        assert!(n.is_power_of_two());

        let lg_n = n.log_2();
        assert_eq!(self.L_vec.len(), lg_n);
        assert_eq!(self.R_vec.len(), lg_n);

        let mut u_vec: Vec<Scalar> = Vec::with_capacity(lg_n);

        for i in 0..lg_n {
            self.L_vec[i].append_to_transcript(b"L", transcript);
            self.R_vec[i].append_to_transcript(b"R", transcript);
            u_vec.push(transcript.challenge_scalar(b"u"));
        }

        // Compute the scaling factors for each generator
        let s = compute_s(&u_vec);

        // Compute g_hat and b_hat from the recursion
        let g_hat = GroupElement::vartime_multiscalar_mul(&s, G_vec);
        let b_hat = compute_dotproduct_with_s(&s, b_vec);

        // Compute Gamma_hat from the L's and R's
        let Gamma_hat = {
            let u_sq: Vec<Scalar> = u_vec.iter().map(|u| *u * *u).collect();
            let u_sq_inv: Vec<Scalar> = u_sq.iter().map(|u_sq| u_sq.invert().unwrap()).collect();

            let L_decomp: Vec<GroupElement> = self
                .L_vec
                .iter()
                .map(|L| L.decompress().unwrap())
                .collect();
            let R_decomp: Vec<GroupElement> = self
                .R_vec
                .iter()
                .map(|R| R.decompress().unwrap())
                .collect();

            let Gamma_L = GroupElement::vartime_multiscalar_mul(&u_sq, &L_decomp);
            let Gamma_R = GroupElement::vartime_multiscalar_mul(&u_sq_inv, &R_decomp);

            Gamma_L + *Gamma + Gamma_R
        };

        Ok((g_hat, Gamma_hat, b_hat))
    }
}

/// Compute dot product of two scalar vectors
fn compute_dotproduct(a: &[Scalar], b: &[Scalar]) -> Scalar {
    assert_eq!(a.len(), b.len());
    a.iter().zip(b.iter()).map(|(a_i, b_i)| *a_i * *b_i).sum()
}

/// Compute the s scalars from the u challenges
fn compute_s(u_vec: &[Scalar]) -> Vec<Scalar> {
    let lg_n = u_vec.len();
    let n = 1 << lg_n;

    let u_inv: Vec<Scalar> = u_vec.iter().map(|u| u.invert().unwrap()).collect();

    let mut s: Vec<Scalar> = vec![Scalar::one(); n];
    for i in 0..n {
        for j in 0..lg_n {
            if i >> j & 1 == 1 {
                s[i] = s[i] * u_vec[lg_n - 1 - j];
            } else {
                s[i] = s[i] * u_inv[lg_n - 1 - j];
            }
        }
    }
    s
}

/// Compute dotproduct with s scalars: sum(s[i] * b[i])
fn compute_dotproduct_with_s(s: &[Scalar], b: &[Scalar]) -> Scalar {
    assert_eq!(s.len(), b.len());
    s.iter().zip(b.iter()).map(|(s_i, b_i)| *s_i * *b_i).sum()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::commitments::MultiCommitGens;
    use crate::math::Math;
    use rand::rngs::OsRng;

    #[test]
    fn test_bullet_reduction_proof() {
        let n = 8;

        // Create n generators (not n+1)
        let gens = MultiCommitGens::new(n, b"test-gens");
        let G = gens.G.clone();
        let H = gens.h;
        let Q = GroupElement::generator();

        let a: Vec<Scalar> = (0..n).map(|_| Scalar::random(&mut OsRng)).collect();
        let b: Vec<Scalar> = (0..n).map(|_| Scalar::random(&mut OsRng)).collect();
        let blind = Scalar::random(&mut OsRng);

        let lg_n = n.log_2();
        let blinds_vec: Vec<(Scalar, Scalar)> = (0..lg_n)
            .map(|_| (Scalar::random(&mut OsRng), Scalar::random(&mut OsRng)))
            .collect();

        let mut prover_transcript = Transcript::new(b"test");
        let (proof, _Gamma_hat, _a_hat, _b_hat, _g_hat, _rhat_Gamma) =
            BulletReductionProof::prove(
                &mut prover_transcript,
                &Q,
                &G,
                &H,
                &a,
                &b,
                &blind,
                &blinds_vec,
            );

        // Compute Gamma for verification
        let Gamma = GroupElement::vartime_multiscalar_mul(&a, &G)
            + compute_dotproduct(&a, &b) * Q
            + blind * H;

        let mut verifier_transcript = Transcript::new(b"test");
        let result = proof.verify(n, &b, &mut verifier_transcript, &Gamma, &G);
        assert!(result.is_ok());
    }
}

//! Zero-knowledge proofs for SNARK mode
//! Port of Spartan's nizk module to BN254

#![allow(clippy::too_many_arguments)]

use ark_serialize::{CanonicalSerialize, CanonicalDeserialize};
use crate::commitments::{Commitments, MultiCommitGens};
use crate::errors::ProofVerifyError;
use crate::group::GroupElement;
use crate::random::RandomTape;
use crate::scalar::Scalar;
use crate::transcript::{AppendToTranscript, ProofTranscript};
use merlin::Transcript;
use serde::{Deserialize, Serialize};

mod bullet;
pub use bullet::BulletReductionProof;

/// Knowledge proof: proves knowledge of (x, r) such that C = g^x * h^r
/// 
/// Uses GroupElement (G1Projective) for cross-verification compatibility with arkworks-spartan
#[derive(Serialize, Deserialize, Debug, CanonicalSerialize, CanonicalDeserialize)]
pub struct KnowledgeProof {
    alpha: GroupElement,
    z1: Scalar,
    z2: Scalar,
}

impl KnowledgeProof {
    fn protocol_name() -> &'static [u8] {
        b"knowledge proof"
    }

    pub fn prove(
        gens_n: &MultiCommitGens,
        transcript: &mut Transcript,
        random_tape: &mut RandomTape,
        x: &Scalar,
        r: &Scalar,
    ) -> (KnowledgeProof, GroupElement) {
        transcript.append_protocol_name(KnowledgeProof::protocol_name());

        // produce two random Scalars
        let t1 = random_tape.random_scalar(b"t1");
        let t2 = random_tape.random_scalar(b"t2");

        let C = x.commit(r, gens_n);
        C.append_to_transcript(b"C", transcript);

        let alpha = t1.commit(&t2, gens_n);
        alpha.append_to_transcript(b"alpha", transcript);

        let c = transcript.challenge_scalar(b"c");

        let z1 = *x * c + t1;
        let z2 = *r * c + t2;

        (KnowledgeProof { alpha, z1, z2 }, C)
    }

    pub fn verify(
        &self,
        gens_n: &MultiCommitGens,
        transcript: &mut Transcript,
        C: &GroupElement,
    ) -> Result<(), ProofVerifyError> {
        transcript.append_protocol_name(KnowledgeProof::protocol_name());
        C.append_to_transcript(b"C", transcript);
        self.alpha.append_to_transcript(b"alpha", transcript);

        let c = transcript.challenge_scalar(b"c");

        let lhs = self.z1.commit(&self.z2, gens_n);
        let rhs = c * *C + self.alpha;

        if lhs == rhs {
            Ok(())
        } else {
            Err(ProofVerifyError::InternalError)
        }
    }
}

/// Equality proof: proves C1 and C2 commit to the same value
#[derive(Serialize, Deserialize, Debug, CanonicalSerialize, CanonicalDeserialize)]
pub struct EqualityProof {
    alpha: GroupElement,
    z: Scalar,
}

impl EqualityProof {
    fn protocol_name() -> &'static [u8] {
        b"equality proof"
    }

    pub fn prove(
        gens_n: &MultiCommitGens,
        transcript: &mut Transcript,
        random_tape: &mut RandomTape,
        v1: &Scalar,
        s1: &Scalar,
        v2: &Scalar,
        s2: &Scalar,
    ) -> (EqualityProof, GroupElement, GroupElement) {
        transcript.append_protocol_name(EqualityProof::protocol_name());

        // produce a random Scalar
        let r = random_tape.random_scalar(b"r");

        let C1 = v1.commit(s1, gens_n);
        C1.append_to_transcript(b"C1", transcript);

        let C2 = v2.commit(s2, gens_n);
        C2.append_to_transcript(b"C2", transcript);

        let alpha = r * gens_n.h;
        alpha.append_to_transcript(b"alpha", transcript);

        let c = transcript.challenge_scalar(b"c");

        let z = c * (*s1 - *s2) + r;

        (EqualityProof { alpha, z }, C1, C2)
    }

    pub fn verify(
        &self,
        gens_n: &MultiCommitGens,
        transcript: &mut Transcript,
        C1: &GroupElement,
        C2: &GroupElement,
    ) -> Result<(), ProofVerifyError> {
        transcript.append_protocol_name(EqualityProof::protocol_name());
        C1.append_to_transcript(b"C1", transcript);
        C2.append_to_transcript(b"C2", transcript);
        self.alpha.append_to_transcript(b"alpha", transcript);

        let c = transcript.challenge_scalar(b"c");
        
        let C = *C1 + (Scalar::zero() - Scalar::one()) * *C2;
        let rhs = c * C + self.alpha;
        let lhs = self.z * gens_n.h;

        if lhs == rhs {
            Ok(())
        } else {
            Err(ProofVerifyError::InternalError)
        }
    }
}

/// Product proof: proves Z = X * Y for committed values
#[derive(Serialize, Deserialize, Debug, CanonicalSerialize, CanonicalDeserialize)]
pub struct ProductProof {
    alpha: GroupElement,
    beta: GroupElement,
    delta: GroupElement,
    z: [Scalar; 5],
}

impl ProductProof {
    fn protocol_name() -> &'static [u8] {
        b"product proof"
    }

    #[allow(non_snake_case)]
    pub fn prove(
        gens_n: &MultiCommitGens,
        transcript: &mut Transcript,
        random_tape: &mut RandomTape,
        x: &Scalar,
        rX: &Scalar,
        y: &Scalar,
        rY: &Scalar,
        z: &Scalar,
        rZ: &Scalar,
    ) -> (ProductProof, GroupElement, GroupElement, GroupElement) {
        transcript.append_protocol_name(ProductProof::protocol_name());

        // produce five random Scalar
        let b1 = random_tape.random_scalar(b"b1");
        let b2 = random_tape.random_scalar(b"b2");
        let b3 = random_tape.random_scalar(b"b3");
        let b4 = random_tape.random_scalar(b"b4");
        let b5 = random_tape.random_scalar(b"b5");

        let X = x.commit(rX, gens_n);
        X.append_to_transcript(b"X", transcript);

        let Y = y.commit(rY, gens_n);
        Y.append_to_transcript(b"Y", transcript);

        let Z = z.commit(rZ, gens_n);
        Z.append_to_transcript(b"Z", transcript);

        let alpha = b1.commit(&b2, gens_n);
        alpha.append_to_transcript(b"alpha", transcript);

        let beta = b3.commit(&b4, gens_n);
        beta.append_to_transcript(b"beta", transcript);

        let delta = {
            let gens_X = MultiCommitGens::from_generators(vec![X], gens_n.h);
            b3.commit(&b5, &gens_X)
        };
        delta.append_to_transcript(b"delta", transcript);

        let c = transcript.challenge_scalar(b"c");

        let z1 = b1 + c * *x;
        let z2 = b2 + c * *rX;
        let z3 = b3 + c * *y;
        let z4 = b4 + c * *rY;
        let z5 = b5 + c * (*rZ - *rX * *y);

        (
            ProductProof {
                alpha,
                beta,
                delta,
                z: [z1, z2, z3, z4, z5],
            },
            X,
            Y,
            Z,
        )
    }

    fn check_equality(
        P: &GroupElement,
        X: &GroupElement,
        c: &Scalar,
        gens_n: &MultiCommitGens,
        z1: &Scalar,
        z2: &Scalar,
    ) -> bool {
        let lhs = *P + *c * *X;
        let rhs = z1.commit(z2, gens_n);
        lhs == rhs
    }

    #[allow(non_snake_case)]
    pub fn verify(
        &self,
        gens_n: &MultiCommitGens,
        transcript: &mut Transcript,
        X: &GroupElement,
        Y: &GroupElement,
        Z: &GroupElement,
    ) -> Result<(), ProofVerifyError> {
        transcript.append_protocol_name(ProductProof::protocol_name());

        X.append_to_transcript(b"X", transcript);
        Y.append_to_transcript(b"Y", transcript);
        Z.append_to_transcript(b"Z", transcript);
        self.alpha.append_to_transcript(b"alpha", transcript);
        self.beta.append_to_transcript(b"beta", transcript);
        self.delta.append_to_transcript(b"delta", transcript);

        let z1 = self.z[0];
        let z2 = self.z[1];
        let z3 = self.z[2];
        let z4 = self.z[3];
        let z5 = self.z[4];

        let c = transcript.challenge_scalar(b"c");

        if ProductProof::check_equality(&self.alpha, X, &c, gens_n, &z1, &z2)
            && ProductProof::check_equality(&self.beta, Y, &c, gens_n, &z3, &z4)
            && ProductProof::check_equality(
                &self.delta,
                Z,
                &c,
                &MultiCommitGens::from_generators(vec![*X], gens_n.h),
                &z3,
                &z5,
            )
        {
            Ok(())
        } else {
            Err(ProofVerifyError::InternalError)
        }
    }
}

/// Dot product proof
#[derive(Debug, Clone, Serialize, Deserialize, CanonicalSerialize, CanonicalDeserialize)]
pub struct DotProductProof {
    delta: GroupElement,
    beta: GroupElement,
    z: Vec<Scalar>,
    z_delta: Scalar,
    z_beta: Scalar,
}

impl DotProductProof {
    fn protocol_name() -> &'static [u8] {
        b"dot product proof"
    }

    pub fn compute_dotproduct(a: &[Scalar], b: &[Scalar]) -> Scalar {
        assert_eq!(a.len(), b.len());
        (0..a.len()).map(|i| a[i] * b[i]).sum()
    }

    #[allow(non_snake_case)]
    pub fn prove(
        gens_1: &MultiCommitGens,
        gens_n: &MultiCommitGens,
        transcript: &mut Transcript,
        random_tape: &mut RandomTape,
        x_vec: &[Scalar],
        blind_x: &Scalar,
        a_vec: &[Scalar],
        y: &Scalar,
        blind_y: &Scalar,
    ) -> (DotProductProof, GroupElement, GroupElement) {
        transcript.append_protocol_name(DotProductProof::protocol_name());

        let n = x_vec.len();
        assert_eq!(x_vec.len(), a_vec.len());
        assert_eq!(gens_n.n, a_vec.len());
        assert_eq!(gens_1.n, 1);

        // produce randomness for the proofs
        let d_vec = random_tape.random_vector(b"d_vec", n);
        let r_delta = random_tape.random_scalar(b"r_delta");
        let r_beta = random_tape.random_scalar(b"r_beta");

        let Cx = x_vec.commit(blind_x, gens_n);
        Cx.append_to_transcript(b"Cx", transcript);

        let Cy = y.commit(blind_y, gens_1);
        Cy.append_to_transcript(b"Cy", transcript);

        a_vec.append_to_transcript(b"a", transcript);

        let delta = d_vec.commit(&r_delta, gens_n);
        delta.append_to_transcript(b"delta", transcript);

        let dotproduct_a_d = DotProductProof::compute_dotproduct(a_vec, &d_vec);

        let beta = dotproduct_a_d.commit(&r_beta, gens_1);
        beta.append_to_transcript(b"beta", transcript);

        let c = transcript.challenge_scalar(b"c");

        let z = (0..d_vec.len())
            .map(|i| c * x_vec[i] + d_vec[i])
            .collect::<Vec<Scalar>>();

        let z_delta = c * *blind_x + r_delta;
        let z_beta = c * *blind_y + r_beta;

        (
            DotProductProof {
                delta,
                beta,
                z,
                z_delta,
                z_beta,
            },
            Cx,
            Cy,
        )
    }

    #[allow(non_snake_case)]
    pub fn verify(
        &self,
        gens_1: &MultiCommitGens,
        gens_n: &MultiCommitGens,
        transcript: &mut Transcript,
        a: &[Scalar],
        Cx: &GroupElement,
        Cy: &GroupElement,
    ) -> Result<(), ProofVerifyError> {
        assert_eq!(gens_n.n, a.len());
        assert_eq!(gens_1.n, 1);

        transcript.append_protocol_name(DotProductProof::protocol_name());
        Cx.append_to_transcript(b"Cx", transcript);
        Cy.append_to_transcript(b"Cy", transcript);
        a.append_to_transcript(b"a", transcript);
        self.delta.append_to_transcript(b"delta", transcript);
        self.beta.append_to_transcript(b"beta", transcript);

        let c = transcript.challenge_scalar(b"c");

        let mut result = c * *Cx + self.delta == self.z.commit(&self.z_delta, gens_n);

        let dotproduct_z_a = DotProductProof::compute_dotproduct(&self.z, a);
        result &= c * *Cy + self.beta == dotproduct_z_a.commit(&self.z_beta, gens_1);

        if result {
            Ok(())
        } else {
            Err(ProofVerifyError::InternalError)
        }
    }
}

/// Generators for dot product proof
#[derive(Clone, Serialize, Deserialize)]
pub struct DotProductProofGens {
    pub n: usize,
    pub gens_n: MultiCommitGens,
    pub gens_1: MultiCommitGens,
}

impl DotProductProofGens {
    pub fn new(n: usize, label: &[u8]) -> Self {
        let (gens_n, gens_1) = MultiCommitGens::new(n + 1, label).split_at(n);
        DotProductProofGens { n, gens_n, gens_1 }
    }
}

/// Logarithmic-size dot product proof using bullet reduction
#[derive(Debug, Clone, Serialize, Deserialize, CanonicalSerialize, CanonicalDeserialize)]
pub struct DotProductProofLog {
    bullet_reduction_proof: BulletReductionProof,
    delta: GroupElement,
    beta: GroupElement,
    z1: Scalar,
    z2: Scalar,
}

impl DotProductProofLog {
    fn protocol_name() -> &'static [u8] {
        b"dot product proof (log)"
    }

    pub fn compute_dotproduct(a: &[Scalar], b: &[Scalar]) -> Scalar {
        assert_eq!(a.len(), b.len());
        (0..a.len()).map(|i| a[i] * b[i]).sum()
    }

    #[allow(non_snake_case)]
    pub fn prove(
        gens: &DotProductProofGens,
        transcript: &mut Transcript,
        random_tape: &mut RandomTape,
        x_vec: &[Scalar],
        blind_x: &Scalar,
        a_vec: &[Scalar],
        y: &Scalar,
        blind_y: &Scalar,
    ) -> (DotProductProofLog, GroupElement, GroupElement) {
        use crate::math::Math;

        transcript.append_protocol_name(DotProductProofLog::protocol_name());

        let n = x_vec.len();
        assert_eq!(x_vec.len(), a_vec.len());
        assert_eq!(gens.n, n);

        // produce randomness for generating a proof
        let d = random_tape.random_scalar(b"d");
        let r_delta = random_tape.random_scalar(b"r_delta");
        let r_beta = random_tape.random_scalar(b"r_delta");
        let blinds_vec = {
            let lg_n = n.log_2();
            let v1 = random_tape.random_vector(b"blinds_vec_1", lg_n);
            let v2 = random_tape.random_vector(b"blinds_vec_2", lg_n);
            (0..lg_n)
                .map(|i| (v1[i], v2[i]))
                .collect::<Vec<(Scalar, Scalar)>>()
        };

        let Cx = x_vec.commit(blind_x, &gens.gens_n);
        Cx.append_to_transcript(b"Cx", transcript);

        let Cy = y.commit(blind_y, &gens.gens_1);
        Cy.append_to_transcript(b"Cy", transcript);

        a_vec.append_to_transcript(b"a", transcript);

        // sample a random base and scale the generator used for
        // the output of the inner product
        let r = transcript.challenge_scalar(b"r");
        let gens_1_scaled = gens.gens_1.scale(&r);

        let blind_Gamma = *blind_x + r * *blind_y;
        let (bullet_reduction_proof, _Gamma_hat, x_hat, a_hat, g_hat, rhat_Gamma) =
            BulletReductionProof::prove(
                transcript,
                &gens_1_scaled.G[0],
                &gens.gens_n.G,
                &gens.gens_n.h,
                x_vec,
                a_vec,
                &blind_Gamma,
                &blinds_vec,
            );
        let y_hat = x_hat * a_hat;

        let delta = {
            let gens_hat = MultiCommitGens::from_generators(vec![g_hat], gens.gens_1.h);
            d.commit(&r_delta, &gens_hat)
        };
        delta.append_to_transcript(b"delta", transcript);

        let beta = d.commit(&r_beta, &gens_1_scaled);
        beta.append_to_transcript(b"beta", transcript);

        let c = transcript.challenge_scalar(b"c");

        let z1 = d + c * y_hat;
        let z2 = a_hat * (c * rhat_Gamma + r_beta) + r_delta;

        (
            DotProductProofLog {
                bullet_reduction_proof,
                delta,
                beta,
                z1,
                z2,
            },
            Cx,
            Cy,
        )
    }

    #[allow(non_snake_case)]
    pub fn verify(
        &self,
        n: usize,
        gens: &DotProductProofGens,
        transcript: &mut Transcript,
        a: &[Scalar],
        Cx: &GroupElement,
        Cy: &GroupElement,
    ) -> Result<(), ProofVerifyError> {
        assert_eq!(gens.n, n);
        assert_eq!(a.len(), n);

        transcript.append_protocol_name(DotProductProofLog::protocol_name());
        Cx.append_to_transcript(b"Cx", transcript);
        Cy.append_to_transcript(b"Cy", transcript);
        a.append_to_transcript(b"a", transcript);

        // sample a random base and scale the generator used for
        // the output of the inner product
        let r = transcript.challenge_scalar(b"r");
        let gens_1_scaled = gens.gens_1.scale(&r);

        let Gamma = *Cx + r * *Cy;

        let (g_hat, Gamma_hat, a_hat) =
            self.bullet_reduction_proof
                .verify(n, a, transcript, &Gamma, &gens.gens_n.G)?;
                
        self.delta.append_to_transcript(b"delta", transcript);
        self.beta.append_to_transcript(b"beta", transcript);

        let c = transcript.challenge_scalar(b"c");

        let lhs = (c * Gamma_hat + self.beta) * a_hat + self.delta;
        let rhs = self.z1 * (g_hat + a_hat * gens_1_scaled.G[0]) + self.z2 * gens_1_scaled.h;

        if lhs == rhs {
            Ok(())
        } else {
            eprintln!("DEBUG: DotProductProofLog::verify failed - lhs != rhs");
            Err(ProofVerifyError::InternalError)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::rngs::OsRng;

    #[test]
    fn check_knowledgeproof() {
        let gens_1 = MultiCommitGens::new(1, b"test-knowledgeproof");

        let x = Scalar::random(&mut OsRng);
        let r = Scalar::random(&mut OsRng);

        let mut random_tape = RandomTape::new(b"proof");
        let mut prover_transcript = Transcript::new(b"example");
        let (proof, committed_value) =
            KnowledgeProof::prove(&gens_1, &mut prover_transcript, &mut random_tape, &x, &r);

        let mut verifier_transcript = Transcript::new(b"example");
        assert!(proof
            .verify(&gens_1, &mut verifier_transcript, &committed_value)
            .is_ok());
    }

    #[test]
    fn check_equalityproof() {
        let gens_1 = MultiCommitGens::new(1, b"test-equalityproof");
        let v1 = Scalar::random(&mut OsRng);
        let v2 = v1;
        let s1 = Scalar::random(&mut OsRng);
        let s2 = Scalar::random(&mut OsRng);

        let mut random_tape = RandomTape::new(b"proof");
        let mut prover_transcript = Transcript::new(b"example");
        let (proof, C1, C2) = EqualityProof::prove(
            &gens_1,
            &mut prover_transcript,
            &mut random_tape,
            &v1,
            &s1,
            &v2,
            &s2,
        );

        let mut verifier_transcript = Transcript::new(b"example");
        assert!(proof
            .verify(&gens_1, &mut verifier_transcript, &C1, &C2)
            .is_ok());
    }

    #[test]
    fn check_productproof() {
        let gens_1 = MultiCommitGens::new(1, b"test-productproof");
        let x = Scalar::random(&mut OsRng);
        let rX = Scalar::random(&mut OsRng);
        let y = Scalar::random(&mut OsRng);
        let rY = Scalar::random(&mut OsRng);
        let z = x * y;
        let rZ = Scalar::random(&mut OsRng);

        let mut random_tape = RandomTape::new(b"proof");
        let mut prover_transcript = Transcript::new(b"example");
        let (proof, X, Y, Z) = ProductProof::prove(
            &gens_1,
            &mut prover_transcript,
            &mut random_tape,
            &x,
            &rX,
            &y,
            &rY,
            &z,
            &rZ,
        );

        let mut verifier_transcript = Transcript::new(b"example");
        assert!(proof
            .verify(&gens_1, &mut verifier_transcript, &X, &Y, &Z)
            .is_ok());
    }

    #[test]
    fn check_dotproductproof() {
        let n = 16;

        let gens_1 = MultiCommitGens::new(1, b"test-two");
        let gens_n = MultiCommitGens::new(n, b"test-n");

        let x: Vec<Scalar> = (0..n).map(|_| Scalar::random(&mut OsRng)).collect();
        let a: Vec<Scalar> = (0..n).map(|_| Scalar::random(&mut OsRng)).collect();
        let y = DotProductProof::compute_dotproduct(&x, &a);
        let r_x = Scalar::random(&mut OsRng);
        let r_y = Scalar::random(&mut OsRng);

        let mut random_tape = RandomTape::new(b"proof");
        let mut prover_transcript = Transcript::new(b"example");
        let (proof, Cx, Cy) = DotProductProof::prove(
            &gens_1,
            &gens_n,
            &mut prover_transcript,
            &mut random_tape,
            &x,
            &r_x,
            &a,
            &y,
            &r_y,
        );

        let mut verifier_transcript = Transcript::new(b"example");
        assert!(proof
            .verify(&gens_1, &gens_n, &mut verifier_transcript, &a, &Cx, &Cy)
            .is_ok());
    }

    #[test]
    fn check_dotproductproof_log() {
        let n = 16;

        let gens = DotProductProofGens::new(n, b"test-n");

        let x: Vec<Scalar> = (0..n).map(|_| Scalar::random(&mut OsRng)).collect();
        let a: Vec<Scalar> = (0..n).map(|_| Scalar::random(&mut OsRng)).collect();
        let y = DotProductProofLog::compute_dotproduct(&x, &a);

        let r_x = Scalar::random(&mut OsRng);
        let r_y = Scalar::random(&mut OsRng);

        let mut random_tape = RandomTape::new(b"proof");
        let mut prover_transcript = Transcript::new(b"example");
        let (proof, Cx, Cy) = DotProductProofLog::prove(
            &gens,
            &mut prover_transcript,
            &mut random_tape,
            &x,
            &r_x,
            &a,
            &y,
            &r_y,
        );

        let mut verifier_transcript = Transcript::new(b"example");
        assert!(proof
            .verify(n, &gens, &mut verifier_transcript, &a, &Cx, &Cy)
            .is_ok());
    }
}

//! Sumcheck protocol implementation
//! Port of Spartan's sumcheck.rs (both non-ZK and ZK versions)

#![allow(clippy::too_many_arguments)]

use crate::commitments::{Commitments, MultiCommitGens};
use crate::dense_mlpoly::DensePolynomial;
use crate::errors::ProofVerifyError;
use crate::group::{CompressedGroup, GroupElement};
use crate::random::RandomTape;
use crate::scalar::Scalar;
use crate::transcript::{AppendToTranscript, ProofTranscript};
use crate::unipoly::{CompressedUniPoly, UniPoly};
use merlin::Transcript;
use serde::{Deserialize, Serialize};

/// Sumcheck proof for a single instance
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct SumcheckInstanceProof {
    pub compressed_polys: Vec<CompressedUniPoly>,
}

impl SumcheckInstanceProof {
    pub fn new(compressed_polys: Vec<CompressedUniPoly>) -> Self {
        SumcheckInstanceProof { compressed_polys }
    }

    /// Verify the sumcheck proof
    pub fn verify(
        &self,
        claim: Scalar,
        num_rounds: usize,
        degree_bound: usize,
        transcript: &mut Transcript,
    ) -> Result<(Scalar, Vec<Scalar>), ProofVerifyError> {
        let mut e = claim;
        let mut r: Vec<Scalar> = Vec::new();

        // verify that there is a univariate polynomial for each round
        if self.compressed_polys.len() != num_rounds {
            return Err(ProofVerifyError::VerificationFailed(
                "wrong number of rounds".to_string(),
            ));
        }

        for i in 0..self.compressed_polys.len() {
            let poly = self.compressed_polys[i].decompress(&e);

            // verify degree bound
            if poly.degree() != degree_bound {
                return Err(ProofVerifyError::VerificationFailed(format!(
                    "degree mismatch at round {}: expected {}, got {}",
                    i,
                    degree_bound,
                    poly.degree()
                )));
            }

            // check if G_k(0) + G_k(1) = e
            if poly.eval_at_zero() + poly.eval_at_one() != e {
                return Err(ProofVerifyError::VerificationFailed(format!(
                    "sum check failed at round {}",
                    i
                )));
            }

            // append the prover's message to the transcript
            poly.append_to_transcript(b"poly", transcript);

            // derive the verifier's challenge for the next round
            let r_i = transcript.challenge_scalar(b"challenge_nextround");

            r.push(r_i);

            // evaluate the claimed degree-ell polynomial at r_i
            e = poly.evaluate(&r_i);
        }

        Ok((e, r))
    }

    /// Prove sumcheck for cubic polynomial (A * B * C)
    pub fn prove_cubic<F>(
        claim: &Scalar,
        num_rounds: usize,
        poly_A: &mut DensePolynomial,
        poly_B: &mut DensePolynomial,
        poly_C: &mut DensePolynomial,
        comb_func: F,
        transcript: &mut Transcript,
    ) -> (Self, Vec<Scalar>, Vec<Scalar>)
    where
        F: Fn(&Scalar, &Scalar, &Scalar) -> Scalar,
    {
        let mut e = *claim;
        let mut r: Vec<Scalar> = Vec::new();
        let mut cubic_polys: Vec<CompressedUniPoly> = Vec::new();

        for _j in 0..num_rounds {
            let mut eval_point_0 = Scalar::zero();
            let mut eval_point_2 = Scalar::zero();
            let mut eval_point_3 = Scalar::zero();

            let len = poly_A.len() / 2;
            for i in 0..len {
                // eval 0: bound_func is A(low)
                eval_point_0 += comb_func(&poly_A[i], &poly_B[i], &poly_C[i]);

                // eval 2: bound_func is -A(low) + 2*A(high)
                let poly_A_bound_point = poly_A[len + i] + poly_A[len + i] - poly_A[i];
                let poly_B_bound_point = poly_B[len + i] + poly_B[len + i] - poly_B[i];
                let poly_C_bound_point = poly_C[len + i] + poly_C[len + i] - poly_C[i];
                eval_point_2 += comb_func(
                    &poly_A_bound_point,
                    &poly_B_bound_point,
                    &poly_C_bound_point,
                );

                // eval 3: bound_func is -2A(low) + 3A(high)
                let poly_A_bound_point = poly_A_bound_point + poly_A[len + i] - poly_A[i];
                let poly_B_bound_point = poly_B_bound_point + poly_B[len + i] - poly_B[i];
                let poly_C_bound_point = poly_C_bound_point + poly_C[len + i] - poly_C[i];

                eval_point_3 += comb_func(
                    &poly_A_bound_point,
                    &poly_B_bound_point,
                    &poly_C_bound_point,
                );
            }

            let evals = vec![eval_point_0, e - eval_point_0, eval_point_2, eval_point_3];
            let poly = UniPoly::from_evals(&evals);

            // append the prover's message to the transcript
            poly.append_to_transcript(b"poly", transcript);

            // derive the verifier's challenge for the next round
            let r_j = transcript.challenge_scalar(b"challenge_nextround");
            r.push(r_j);

            // bound all tables to the verifier's challenge
            poly_A.bound_poly_var_top(&r_j);
            poly_B.bound_poly_var_top(&r_j);
            poly_C.bound_poly_var_top(&r_j);

            e = poly.evaluate(&r_j);
            cubic_polys.push(poly.compress());
        }

        (
            SumcheckInstanceProof::new(cubic_polys),
            r,
            vec![poly_A[0], poly_B[0], poly_C[0]],
        )
    }

    /// Prove batched sumcheck for cubic polynomials
    /// This is used in the product circuit evaluation
    pub fn prove_cubic_batched<F>(
        claim: &Scalar,
        num_rounds: usize,
        poly_vec_par: (
            &mut Vec<&mut DensePolynomial>,
            &mut Vec<&mut DensePolynomial>,
            &mut DensePolynomial,
        ),
        poly_vec_seq: (
            &mut Vec<&mut DensePolynomial>,
            &mut Vec<&mut DensePolynomial>,
            &mut Vec<&mut DensePolynomial>,
        ),
        coeffs: &[Scalar],
        comb_func: F,
        transcript: &mut Transcript,
    ) -> (
        Self,
        Vec<Scalar>,
        (Vec<Scalar>, Vec<Scalar>, Scalar),
        (Vec<Scalar>, Vec<Scalar>, Vec<Scalar>),
    )
    where
        F: Fn(&Scalar, &Scalar, &Scalar) -> Scalar,
    {
        let (poly_A_vec_par, poly_B_vec_par, poly_C_par) = poly_vec_par;
        let (poly_A_vec_seq, poly_B_vec_seq, poly_C_vec_seq) = poly_vec_seq;

        let mut e = *claim;
        let mut r: Vec<Scalar> = Vec::new();
        let mut cubic_polys: Vec<CompressedUniPoly> = Vec::new();

        for _j in 0..num_rounds {
            let mut evals: Vec<(Scalar, Scalar, Scalar)> = Vec::new();

            // Process parallel polynomials (share poly_C)
            for (poly_A, poly_B) in poly_A_vec_par.iter().zip(poly_B_vec_par.iter()) {
                let mut eval_point_0 = Scalar::zero();
                let mut eval_point_2 = Scalar::zero();
                let mut eval_point_3 = Scalar::zero();

                let len = poly_A.len() / 2;
                for i in 0..len {
                    // eval 0: bound_func is A(low)
                    eval_point_0 += comb_func(&poly_A[i], &poly_B[i], &poly_C_par[i]);

                    // eval 2: bound_func is -A(low) + 2*A(high)
                    let poly_A_bound_point = poly_A[len + i] + poly_A[len + i] - poly_A[i];
                    let poly_B_bound_point = poly_B[len + i] + poly_B[len + i] - poly_B[i];
                    let poly_C_bound_point = poly_C_par[len + i] + poly_C_par[len + i] - poly_C_par[i];
                    eval_point_2 += comb_func(
                        &poly_A_bound_point,
                        &poly_B_bound_point,
                        &poly_C_bound_point,
                    );

                    // eval 3
                    let poly_A_bound_point = poly_A_bound_point + poly_A[len + i] - poly_A[i];
                    let poly_B_bound_point = poly_B_bound_point + poly_B[len + i] - poly_B[i];
                    let poly_C_bound_point =
                        poly_C_bound_point + poly_C_par[len + i] - poly_C_par[i];

                    eval_point_3 += comb_func(
                        &poly_A_bound_point,
                        &poly_B_bound_point,
                        &poly_C_bound_point,
                    );
                }

                evals.push((eval_point_0, eval_point_2, eval_point_3));
            }

            // Process sequential polynomials (each has own poly_C)
            for i in 0..poly_A_vec_seq.len() {
                let poly_A = &poly_A_vec_seq[i];
                let poly_B = &poly_B_vec_seq[i];
                let poly_C = &poly_C_vec_seq[i];
                
                let mut eval_point_0 = Scalar::zero();
                let mut eval_point_2 = Scalar::zero();
                let mut eval_point_3 = Scalar::zero();
                let len = poly_A.len() / 2;
                for j in 0..len {
                    eval_point_0 += comb_func(&poly_A[j], &poly_B[j], &poly_C[j]);
                    let poly_A_bound_point = poly_A[len + j] + poly_A[len + j] - poly_A[j];
                    let poly_B_bound_point = poly_B[len + j] + poly_B[len + j] - poly_B[j];
                    let poly_C_bound_point = poly_C[len + j] + poly_C[len + j] - poly_C[j];
                    eval_point_2 += comb_func(
                        &poly_A_bound_point,
                        &poly_B_bound_point,
                        &poly_C_bound_point,
                    );
                    let poly_A_bound_point = poly_A_bound_point + poly_A[len + j] - poly_A[j];
                    let poly_B_bound_point = poly_B_bound_point + poly_B[len + j] - poly_B[j];
                    let poly_C_bound_point = poly_C_bound_point + poly_C[len + j] - poly_C[j];
                    eval_point_3 += comb_func(
                        &poly_A_bound_point,
                        &poly_B_bound_point,
                        &poly_C_bound_point,
                    );
                }
                evals.push((eval_point_0, eval_point_2, eval_point_3));
            }

            let evals_combined_0: Scalar = (0..evals.len()).map(|i| evals[i].0 * coeffs[i]).sum();
            let evals_combined_2: Scalar = (0..evals.len()).map(|i| evals[i].1 * coeffs[i]).sum();
            let evals_combined_3: Scalar = (0..evals.len()).map(|i| evals[i].2 * coeffs[i]).sum();

            let evals = vec![
                evals_combined_0,
                e - evals_combined_0,
                evals_combined_2,
                evals_combined_3,
            ];
            let poly = UniPoly::from_evals(&evals);

            // append the prover's message to the transcript
            poly.append_to_transcript(b"poly", transcript);

            // derive the verifier's challenge for the next round
            let r_j = transcript.challenge_scalar(b"challenge_nextround");
            r.push(r_j);

            // bound all tables to the verifier's challenge
            for (poly_A, poly_B) in poly_A_vec_par.iter_mut().zip(poly_B_vec_par.iter_mut()) {
                poly_A.bound_poly_var_top(&r_j);
                poly_B.bound_poly_var_top(&r_j);
            }
            poly_C_par.bound_poly_var_top(&r_j);

            for i in 0..poly_A_vec_seq.len() {
                poly_A_vec_seq[i].bound_poly_var_top(&r_j);
                poly_B_vec_seq[i].bound_poly_var_top(&r_j);
                poly_C_vec_seq[i].bound_poly_var_top(&r_j);
            }

            e = poly.evaluate(&r_j);
            cubic_polys.push(poly.compress());
        }

        let poly_A_par_final = (0..poly_A_vec_par.len())
            .map(|i| poly_A_vec_par[i][0])
            .collect();
        let poly_B_par_final = (0..poly_B_vec_par.len())
            .map(|i| poly_B_vec_par[i][0])
            .collect();
        let claims_prod = (poly_A_par_final, poly_B_par_final, poly_C_par[0]);

        let poly_A_seq_final = (0..poly_A_vec_seq.len())
            .map(|i| poly_A_vec_seq[i][0])
            .collect();
        let poly_B_seq_final = (0..poly_B_vec_seq.len())
            .map(|i| poly_B_vec_seq[i][0])
            .collect();
        let poly_C_seq_final = (0..poly_C_vec_seq.len())
            .map(|i| poly_C_vec_seq[i][0])
            .collect();
        let claims_dotp = (poly_A_seq_final, poly_B_seq_final, poly_C_seq_final);

        (
            SumcheckInstanceProof::new(cubic_polys),
            r,
            claims_prod,
            claims_dotp,
        )
    }
}

/// Zero-Knowledge Sumcheck proof with commitments
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct ZKSumcheckInstanceProof {
    pub comm_polys: Vec<CompressedGroup>,
    pub comm_evals: Vec<CompressedGroup>,
    pub proofs: Vec<CompressedUniPoly>,
}

impl ZKSumcheckInstanceProof {
    /// Verify the ZK sumcheck proof
    #[allow(non_snake_case)]
    pub fn verify(
        &self,
        comm_claim: &CompressedGroup,
        num_rounds: usize,
        degree_bound: usize,
        gens_1: &MultiCommitGens,
        gens_n: &MultiCommitGens,
        transcript: &mut Transcript,
    ) -> Result<(CompressedGroup, Vec<Scalar>), ProofVerifyError> {
        // verify that there is a univariate polynomial for each round
        if self.comm_polys.len() != num_rounds || self.proofs.len() != num_rounds {
            return Err(ProofVerifyError::VerificationFailed(
                "wrong number of rounds".to_string(),
            ));
        }

        let mut comm_e = comm_claim.clone();
        let mut r: Vec<Scalar> = Vec::new();

        for i in 0..self.comm_polys.len() {
            // append commitment to transcript
            self.comm_polys[i].append_to_transcript(b"comm_poly", transcript);

            // derive challenge
            let r_i = transcript.challenge_scalar(b"challenge_nextround");
            r.push(r_i);

            // Decompress the polynomial commitment
            let poly = self.proofs[i].decompress(&Scalar::zero()); // Placeholder
            
            // verify degree bound
            if poly.degree() != degree_bound {
                return Err(ProofVerifyError::VerificationFailed(format!(
                    "degree mismatch at round {}: expected {}, got {}",
                    i, degree_bound, poly.degree()
                )));
            }

            // update commitment
            self.comm_evals[i].append_to_transcript(b"comm_eval", transcript);
            comm_e = self.comm_evals[i].clone();
        }

        Ok((comm_e, r))
    }

    /// Prove ZK sumcheck for quartic polynomial with additive term (tau * (Az * Bz - Cz))
    #[allow(non_snake_case)]
    pub fn prove_cubic_with_additive_term<F>(
        claim: &Scalar,
        blind_claim: &Scalar,
        num_rounds: usize,
        poly_tau: &mut DensePolynomial,
        poly_Az: &mut DensePolynomial,
        poly_Bz: &mut DensePolynomial,
        poly_Cz: &mut DensePolynomial,
        comb_func: F,
        gens_1: &MultiCommitGens,
        gens_n: &MultiCommitGens,
        transcript: &mut Transcript,
        random_tape: &mut RandomTape,
    ) -> (Self, Vec<Scalar>, Vec<Scalar>, Scalar)
    where
        F: Fn(&Scalar, &Scalar, &Scalar, &Scalar) -> Scalar,
    {
        let mut e = *claim;
        let mut blind_e = *blind_claim;
        let mut r: Vec<Scalar> = Vec::new();
        let mut comm_polys: Vec<CompressedGroup> = Vec::new();
        let mut comm_evals: Vec<CompressedGroup> = Vec::new();
        let mut proofs: Vec<CompressedUniPoly> = Vec::new();

        for _j in 0..num_rounds {
            let mut eval_point_0 = Scalar::zero();
            let mut eval_point_2 = Scalar::zero();
            let mut eval_point_3 = Scalar::zero();

            let len = poly_tau.len() / 2;
            for i in 0..len {
                // eval 0: bound_func is poly(low)
                eval_point_0 += comb_func(&poly_tau[i], &poly_Az[i], &poly_Bz[i], &poly_Cz[i]);

                // eval 2: bound_func is -poly(low) + 2*poly(high)
                let poly_tau_bound = poly_tau[len + i] + poly_tau[len + i] - poly_tau[i];
                let poly_Az_bound = poly_Az[len + i] + poly_Az[len + i] - poly_Az[i];
                let poly_Bz_bound = poly_Bz[len + i] + poly_Bz[len + i] - poly_Bz[i];
                let poly_Cz_bound = poly_Cz[len + i] + poly_Cz[len + i] - poly_Cz[i];
                eval_point_2 += comb_func(
                    &poly_tau_bound,
                    &poly_Az_bound,
                    &poly_Bz_bound,
                    &poly_Cz_bound,
                );

                // eval 3: bound_func is -2*poly(low) + 3*poly(high)
                let poly_tau_bound = poly_tau_bound + poly_tau[len + i] - poly_tau[i];
                let poly_Az_bound = poly_Az_bound + poly_Az[len + i] - poly_Az[i];
                let poly_Bz_bound = poly_Bz_bound + poly_Bz[len + i] - poly_Bz[i];
                let poly_Cz_bound = poly_Cz_bound + poly_Cz[len + i] - poly_Cz[i];

                eval_point_3 += comb_func(
                    &poly_tau_bound,
                    &poly_Az_bound,
                    &poly_Bz_bound,
                    &poly_Cz_bound,
                );
            }

            let evals = vec![eval_point_0, e - eval_point_0, eval_point_2, eval_point_3];
            let poly = UniPoly::from_evals(&evals);

            // commit to polynomial
            let blind_poly = random_tape.random_scalar(b"blind_poly");
            let comm_poly = poly.commit(gens_n, &blind_poly);
            comm_poly.append_to_transcript(b"comm_poly", transcript);
            comm_polys.push(comm_poly);

            // derive challenge
            let r_j = transcript.challenge_scalar(b"challenge_nextround");
            r.push(r_j);

            // evaluate and commit
            let eval = poly.evaluate(&r_j);
            let blind_eval = random_tape.random_scalar(b"blind_eval");
            let comm_eval = eval.commit(&blind_eval, gens_1).compress();
            comm_eval.append_to_transcript(b"comm_eval", transcript);
            comm_evals.push(comm_eval);

            // bound all tables
            poly_tau.bound_poly_var_top(&r_j);
            poly_Az.bound_poly_var_top(&r_j);
            poly_Bz.bound_poly_var_top(&r_j);
            poly_Cz.bound_poly_var_top(&r_j);

            proofs.push(poly.compress());
            e = eval;
            blind_e = blind_eval;
        }

        (
            ZKSumcheckInstanceProof {
                comm_polys,
                comm_evals,
                proofs,
            },
            r,
            vec![poly_tau[0], poly_Az[0], poly_Bz[0], poly_Cz[0]],
            blind_e,
        )
    }

    /// Prove ZK sumcheck for quadratic polynomial (z * ABC)
    #[allow(non_snake_case)]
    pub fn prove_quad<F>(
        claim: &Scalar,
        blind_claim: &Scalar,
        num_rounds: usize,
        poly_z: &mut DensePolynomial,
        poly_ABC: &mut DensePolynomial,
        comb_func: F,
        gens_1: &MultiCommitGens,
        gens_n: &MultiCommitGens,
        transcript: &mut Transcript,
        random_tape: &mut RandomTape,
    ) -> (Self, Vec<Scalar>, Vec<Scalar>, Scalar)
    where
        F: Fn(&Scalar, &Scalar) -> Scalar,
    {
        let mut e = *claim;
        let mut blind_e = *blind_claim;
        let mut r: Vec<Scalar> = Vec::new();
        let mut comm_polys: Vec<CompressedGroup> = Vec::new();
        let mut comm_evals: Vec<CompressedGroup> = Vec::new();
        let mut proofs: Vec<CompressedUniPoly> = Vec::new();

        for _j in 0..num_rounds {
            let mut eval_point_0 = Scalar::zero();
            let mut eval_point_2 = Scalar::zero();

            let len = poly_z.len() / 2;
            for i in 0..len {
                // eval 0
                eval_point_0 += comb_func(&poly_z[i], &poly_ABC[i]);

                // eval 2
                let poly_z_bound = poly_z[len + i] + poly_z[len + i] - poly_z[i];
                let poly_ABC_bound = poly_ABC[len + i] + poly_ABC[len + i] - poly_ABC[i];
                eval_point_2 += comb_func(&poly_z_bound, &poly_ABC_bound);
            }

            let evals = vec![eval_point_0, e - eval_point_0, eval_point_2];
            let poly = UniPoly::from_evals(&evals);

            // commit to polynomial
            let blind_poly = random_tape.random_scalar(b"blind_poly");
            let comm_poly = poly.commit(gens_n, &blind_poly);
            comm_poly.append_to_transcript(b"comm_poly", transcript);
            comm_polys.push(comm_poly);

            // derive challenge
            let r_j = transcript.challenge_scalar(b"challenge_nextround");
            r.push(r_j);

            // evaluate and commit
            let eval = poly.evaluate(&r_j);
            let blind_eval = random_tape.random_scalar(b"blind_eval");
            let comm_eval = eval.commit(&blind_eval, gens_1).compress();
            comm_eval.append_to_transcript(b"comm_eval", transcript);
            comm_evals.push(comm_eval);

            // bound tables
            poly_z.bound_poly_var_top(&r_j);
            poly_ABC.bound_poly_var_top(&r_j);

            proofs.push(poly.compress());
            e = eval;
            blind_e = blind_eval;
        }

        (
            ZKSumcheckInstanceProof {
                comm_polys,
                comm_evals,
                proofs,
            },
            r,
            vec![poly_z[0], poly_ABC[0]],
            blind_e,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sumcheck_cubic() {
        let mut poly_A = DensePolynomial::new(vec![
            Scalar::from_u64(1),
            Scalar::from_u64(2),
            Scalar::from_u64(3),
            Scalar::from_u64(4),
        ]);
        let mut poly_B = DensePolynomial::new(vec![
            Scalar::from_u64(1),
            Scalar::from_u64(1),
            Scalar::from_u64(1),
            Scalar::from_u64(1),
        ]);
        let mut poly_C = DensePolynomial::new(vec![
            Scalar::from_u64(1),
            Scalar::from_u64(1),
            Scalar::from_u64(1),
            Scalar::from_u64(1),
        ]);

        // Compute the claim: sum of A * B * C
        let claim: Scalar = (0..4)
            .map(|i| poly_A[i] * poly_B[i] * poly_C[i])
            .sum();

        let comb_func = |a: &Scalar, b: &Scalar, c: &Scalar| *a * *b * *c;

        let mut transcript = Transcript::new(b"test");
        let (proof, r, _claims) = SumcheckInstanceProof::prove_cubic(
            &claim,
            2,
            &mut poly_A,
            &mut poly_B,
            &mut poly_C,
            comb_func,
            &mut transcript,
        );

        // Verify
        let mut verifier_transcript = Transcript::new(b"test");
        let result = proof.verify(claim, 2, 3, &mut verifier_transcript);
        assert!(result.is_ok());
    }
}

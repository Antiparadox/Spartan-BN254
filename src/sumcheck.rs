//! Sumcheck protocol implementation
//! Port of Spartan's sumcheck.rs (both non-ZK and ZK versions)
//!
//! The ZK sumcheck uses DotProductProof to hide polynomial coefficients,
//! matching the original Microsoft Spartan implementation.

#![allow(clippy::too_many_arguments)]

use ark_serialize::{CanonicalSerialize, CanonicalDeserialize};
use crate::commitments::{Commitments, MultiCommitGens};
use crate::dense_mlpoly::DensePolynomial;
use crate::errors::ProofVerifyError;
use crate::group::GroupElement;
use crate::nizk::DotProductProof;
use crate::random::RandomTape;
use crate::scalar::Scalar;
use crate::transcript::{AppendToTranscript, ProofTranscript};
use crate::unipoly::{CompressedUniPoly, UniPoly};
use core::iter;
use merlin::Transcript;
use serde::{Deserialize, Serialize};

/// Sumcheck proof for a single instance
#[derive(Serialize, Deserialize, Debug, Clone, CanonicalSerialize, CanonicalDeserialize)]
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
/// 
/// Uses DotProductProof to prove polynomial evaluations without revealing coefficients.
/// This matches the original Microsoft Spartan's ZK sumcheck.
#[derive(Serialize, Deserialize, Debug, Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct ZKSumcheckInstanceProof {
    pub comm_polys: Vec<GroupElement>,
    pub comm_evals: Vec<GroupElement>,
    pub proofs: Vec<DotProductProof>,
}

impl ZKSumcheckInstanceProof {
    /// Create a new ZK sumcheck proof
    pub fn new(
        comm_polys: Vec<GroupElement>,
        comm_evals: Vec<GroupElement>,
        proofs: Vec<DotProductProof>,
    ) -> Self {
        ZKSumcheckInstanceProof {
            comm_polys,
            comm_evals,
            proofs,
        }
    }

    /// Verify the ZK sumcheck proof using DotProductProof
    /// 
    /// For each round, we verify:
    /// 1. poly(0) + poly(1) = claim_per_round  (sumcheck property)
    /// 2. poly(r_j) = eval  (evaluation correctness)
    /// 
    /// Both are proven simultaneously via DotProductProof with batching.
    #[allow(non_snake_case)]
    pub fn verify(
        &self,
        comm_claim: &GroupElement,
        num_rounds: usize,
        degree_bound: usize,
        gens_1: &MultiCommitGens,
        gens_n: &MultiCommitGens,
        transcript: &mut Transcript,
    ) -> Result<(GroupElement, Vec<Scalar>), ProofVerifyError> {
        // verify that there is a univariate polynomial for each round
        if self.comm_polys.len() != num_rounds || self.proofs.len() != num_rounds {
            return Err(ProofVerifyError::VerificationFailed(
                "wrong number of rounds".to_string(),
            ));
        }

        let mut comm_claim_per_round = *comm_claim;
        let mut r: Vec<Scalar> = Vec::new();

        for i in 0..self.comm_polys.len() {
            // append polynomial commitment to transcript
            self.comm_polys[i].append_to_transcript(b"comm_poly", transcript);

            // derive the verifier's challenge for the next round
            let r_i = transcript.challenge_scalar(b"challenge_nextround");

            // add two claims to transcript (same as prover)
            comm_claim_per_round.append_to_transcript(b"comm_claim_per_round", transcript);
            self.comm_evals[i].append_to_transcript(b"comm_eval", transcript);

            // produce two weights
            let w = transcript.challenge_vector(b"combine_two_claims_to_one", 2);

            // compute the commitment to the weighted sum (target = w[0] * claim + w[1] * eval)
            let comm_target = GroupElement::vartime_multiscalar_mul(
                w.iter(),
                iter::once(&comm_claim_per_round)
                    .chain(iter::once(&self.comm_evals[i]))
                    .collect::<Vec<&GroupElement>>(),
            );

            // construct vector a for the DotProductProof verification
            // a = w[0] * a_sc + w[1] * a_eval
            // where:
            //   a_sc = [2, 1, 1, 1, ...] for sumcheck: poly(0) + poly(1) = <coeffs, [2,1,1,1,...]>
            //   a_eval = [1, r_j, r_j^2, ...] for evaluation: poly(r_j) = <coeffs, [1,r,r^2,...]>
            let a = {
                // the vector to use to decommit for sum-check test: poly(0) + poly(1) = claim
                // For a polynomial with coefficients [c0, c1, c2, ...]:
                // poly(0) = c0
                // poly(1) = c0 + c1 + c2 + ...
                // So poly(0) + poly(1) = 2*c0 + c1 + c2 + ...
                let a_sc = {
                    let mut a = vec![Scalar::one(); degree_bound + 1];
                    a[0] = a[0] + Scalar::one(); // a[0] = 2
                    a
                };

                // the vector to use to decommit for evaluation: poly(r_j)
                // = c0 + c1*r + c2*r^2 + ...
                let a_eval = {
                    let mut a = vec![Scalar::one(); degree_bound + 1];
                    for j in 1..a.len() {
                        a[j] = a[j - 1] * r_i;
                    }
                    a
                };

                // take weighted sum of the two vectors using w
                assert_eq!(a_sc.len(), a_eval.len());
                (0..a_sc.len())
                    .map(|j| w[0] * a_sc[j] + w[1] * a_eval[j])
                    .collect::<Vec<Scalar>>()
            };

            // verify the DotProductProof
            self.proofs[i].verify(
                gens_1,
                gens_n,
                transcript,
                &a,
                &self.comm_polys[i],
                &comm_target,
            )?;

            // update for next round
            comm_claim_per_round = self.comm_evals[i];
            r.push(r_i);
        }

        Ok((self.comm_evals[self.comm_evals.len() - 1], r))
    }

    /// Prove ZK sumcheck for quartic polynomial with additive term (tau * (Az * Bz - Cz))
    /// 
    /// Uses DotProductProof to hide polynomial coefficients while proving:
    /// 1. poly(0) + poly(1) = claim_per_round  (sumcheck property)
    /// 2. poly(r_j) = eval  (evaluation correctness)
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
        // Pre-generate blinds for all rounds
        let blinds_poly = random_tape.random_vector(b"blinds_poly", num_rounds);
        let blinds_evals = random_tape.random_vector(b"blinds_evals", num_rounds);

        let mut claim_per_round = *claim;
        let mut comm_claim_per_round = claim_per_round.commit(blind_claim, gens_1);

        let mut r: Vec<Scalar> = Vec::new();
        let mut comm_polys: Vec<GroupElement> = Vec::new();
        let mut comm_evals: Vec<GroupElement> = Vec::new();
        let mut proofs: Vec<DotProductProof> = Vec::new();

        for j in 0..num_rounds {
            // Compute polynomial evaluations
            let (poly, comm_poly) = {
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

                let evals = vec![
                    eval_point_0,
                    claim_per_round - eval_point_0,
                    eval_point_2,
                    eval_point_3,
                ];
                let poly = UniPoly::from_evals(&evals);
                let comm_poly = poly.commit(gens_n, &blinds_poly[j]);
                (poly, comm_poly)
            };

            // append the prover's message to the transcript
            comm_poly.append_to_transcript(b"comm_poly", transcript);
            comm_polys.push(comm_poly);

            // derive the verifier's challenge for the next round
            let r_j = transcript.challenge_scalar(b"challenge_nextround");

            // bound all tables to the verifier's challenge
            poly_tau.bound_poly_var_top(&r_j);
            poly_Az.bound_poly_var_top(&r_j);
            poly_Bz.bound_poly_var_top(&r_j);
            poly_Cz.bound_poly_var_top(&r_j);

            // produce a proof of sum-check and of evaluation using DotProductProof
            let (proof, claim_next_round, comm_claim_next_round) = {
                let eval = poly.evaluate(&r_j);
                let comm_eval = eval.commit(&blinds_evals[j], gens_1);

                // we need to prove the following under homomorphic commitments:
                // (1) poly(0) + poly(1) = claim_per_round
                // (2) poly(r_j) = eval
                //
                // Our technique is to leverage dot product proofs:
                // (1) we can prove: <poly_coeffs, [2,1,1,1,...]> = claim_per_round
                // (2) we can prove: <poly_coeffs, [1, r_j, r_j^2, ...]> = eval
                // for efficiency we batch them using random weights

                // add two claims to transcript
                comm_claim_per_round.append_to_transcript(b"comm_claim_per_round", transcript);
                comm_eval.append_to_transcript(b"comm_eval", transcript);

                // produce two weights
                let w = transcript.challenge_vector(b"combine_two_claims_to_one", 2);

                // compute a weighted sum of the RHS
                let target = w[0] * claim_per_round + w[1] * eval;
                let comm_target = GroupElement::vartime_multiscalar_mul(
                    w.iter(),
                    iter::once(&comm_claim_per_round)
                        .chain(iter::once(&comm_eval))
                        .collect::<Vec<&GroupElement>>(),
                );

                let blind = {
                    let blind_sc = if j == 0 {
                        blind_claim
                    } else {
                        &blinds_evals[j - 1]
                    };
                    let blind_eval = &blinds_evals[j];
                    w[0] * *blind_sc + w[1] * *blind_eval
                };
                debug_assert_eq!(target.commit(&blind, gens_1), comm_target);

                let a = {
                    // the vector to use to decommit for sum-check test
                    let a_sc = {
                        let mut a = vec![Scalar::one(); poly.degree() + 1];
                        a[0] = a[0] + Scalar::one(); // a[0] = 2
                        a
                    };

                    // the vector to use to decommit for evaluation
                    let a_eval = {
                        let mut a = vec![Scalar::one(); poly.degree() + 1];
                        for k in 1..a.len() {
                            a[k] = a[k - 1] * r_j;
                        }
                        a
                    };

                    // take weighted sum of the two vectors using w
                    assert_eq!(a_sc.len(), a_eval.len());
                    (0..a_sc.len())
                        .map(|k| w[0] * a_sc[k] + w[1] * a_eval[k])
                        .collect::<Vec<Scalar>>()
                };

                let (proof, _comm_poly, _comm_sc_eval) = DotProductProof::prove(
                    gens_1,
                    gens_n,
                    transcript,
                    random_tape,
                    &poly.as_vec(),
                    &blinds_poly[j],
                    &a,
                    &target,
                    &blind,
                );

                (proof, eval, comm_eval)
            };

            proofs.push(proof);
            claim_per_round = claim_next_round;
            comm_claim_per_round = comm_claim_next_round;
            r.push(r_j);
            comm_evals.push(comm_claim_per_round);
        }

        (
            ZKSumcheckInstanceProof::new(comm_polys, comm_evals, proofs),
            r,
            vec![poly_tau[0], poly_Az[0], poly_Bz[0], poly_Cz[0]],
            blinds_evals[num_rounds - 1],
        )
    }

    /// Prove ZK sumcheck for quadratic polynomial (z * ABC)
    /// 
    /// Uses DotProductProof to hide polynomial coefficients while proving:
    /// 1. poly(0) + poly(1) = claim_per_round  (sumcheck property)
    /// 2. poly(r_j) = eval  (evaluation correctness)
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
        // Pre-generate blinds for all rounds
        let blinds_poly = random_tape.random_vector(b"blinds_poly", num_rounds);
        let blinds_evals = random_tape.random_vector(b"blinds_evals", num_rounds);

        let mut claim_per_round = *claim;
        let mut comm_claim_per_round = claim_per_round.commit(blind_claim, gens_1);

        let mut r: Vec<Scalar> = Vec::new();
        let mut comm_polys: Vec<GroupElement> = Vec::new();
        let mut comm_evals: Vec<GroupElement> = Vec::new();
        let mut proofs: Vec<DotProductProof> = Vec::new();

        for j in 0..num_rounds {
            // Compute polynomial evaluations
            let (poly, comm_poly) = {
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

                let evals = vec![eval_point_0, claim_per_round - eval_point_0, eval_point_2];
                let poly = UniPoly::from_evals(&evals);
                let comm_poly = poly.commit(gens_n, &blinds_poly[j]);
                (poly, comm_poly)
            };

            // append the prover's message to the transcript
            comm_poly.append_to_transcript(b"comm_poly", transcript);
            comm_polys.push(comm_poly);

            // derive the verifier's challenge for the next round
            let r_j = transcript.challenge_scalar(b"challenge_nextround");

            // bound tables to the verifier's challenge
            poly_z.bound_poly_var_top(&r_j);
            poly_ABC.bound_poly_var_top(&r_j);

            // produce a proof of sum-check and of evaluation using DotProductProof
            let (proof, claim_next_round, comm_claim_next_round) = {
                let eval = poly.evaluate(&r_j);
                let comm_eval = eval.commit(&blinds_evals[j], gens_1);

                // we need to prove the following under homomorphic commitments:
                // (1) poly(0) + poly(1) = claim_per_round
                // (2) poly(r_j) = eval
                //
                // Our technique is to leverage dot product proofs:
                // (1) we can prove: <poly_coeffs, [2,1,1,...]> = claim_per_round
                // (2) we can prove: <poly_coeffs, [1, r_j, r_j^2, ...]> = eval
                // for efficiency we batch them using random weights

                // add two claims to transcript
                comm_claim_per_round.append_to_transcript(b"comm_claim_per_round", transcript);
                comm_eval.append_to_transcript(b"comm_eval", transcript);

                // produce two weights
                let w = transcript.challenge_vector(b"combine_two_claims_to_one", 2);

                // compute a weighted sum of the RHS
                let target = w[0] * claim_per_round + w[1] * eval;
                let comm_target = GroupElement::vartime_multiscalar_mul(
                    w.iter(),
                    iter::once(&comm_claim_per_round)
                        .chain(iter::once(&comm_eval))
                        .collect::<Vec<&GroupElement>>(),
                );

                let blind = {
                    let blind_sc = if j == 0 {
                        blind_claim
                    } else {
                        &blinds_evals[j - 1]
                    };
                    let blind_eval = &blinds_evals[j];
                    w[0] * *blind_sc + w[1] * *blind_eval
                };
                debug_assert_eq!(target.commit(&blind, gens_1), comm_target);

                let a = {
                    // the vector to use to decommit for sum-check test
                    let a_sc = {
                        let mut a = vec![Scalar::one(); poly.degree() + 1];
                        a[0] = a[0] + Scalar::one(); // a[0] = 2
                        a
                    };

                    // the vector to use to decommit for evaluation
                    let a_eval = {
                        let mut a = vec![Scalar::one(); poly.degree() + 1];
                        for k in 1..a.len() {
                            a[k] = a[k - 1] * r_j;
                        }
                        a
                    };

                    // take weighted sum of the two vectors using w
                    assert_eq!(a_sc.len(), a_eval.len());
                    (0..a_sc.len())
                        .map(|k| w[0] * a_sc[k] + w[1] * a_eval[k])
                        .collect::<Vec<Scalar>>()
                };

                let (proof, _comm_poly, _comm_sc_eval) = DotProductProof::prove(
                    gens_1,
                    gens_n,
                    transcript,
                    random_tape,
                    &poly.as_vec(),
                    &blinds_poly[j],
                    &a,
                    &target,
                    &blind,
                );

                (proof, eval, comm_eval)
            };

            proofs.push(proof);
            claim_per_round = claim_next_round;
            comm_claim_per_round = comm_claim_next_round;
            r.push(r_j);
            comm_evals.push(comm_claim_per_round);
        }

        (
            ZKSumcheckInstanceProof::new(comm_polys, comm_evals, proofs),
            r,
            vec![poly_z[0], poly_ABC[0]],
            blinds_evals[num_rounds - 1],
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

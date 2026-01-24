//! R1CS Proof - combines sumcheck with R1CS
//! Port of Spartan's r1csproof.rs

#![allow(clippy::too_many_arguments)]

use ark_serialize::{CanonicalSerialize, CanonicalDeserialize};
use crate::commitments::{Commitments, MultiCommitGens};
use crate::hyrax::{DensePolynomial, EqPolynomial};
use crate::errors::ProofVerifyError;
use crate::group::GroupElement;
use crate::math::Math;
use crate::nizk::{DotProductProofGens, DotProductProofLog, EqualityProof, KnowledgeProof, ProductProof};
use crate::r1cs::R1CSShape;
use crate::random::RandomTape;
use crate::scalar::Scalar;
use crate::sumcheck::ZKSumcheckInstanceProof;
use crate::timer::Timer;
use crate::transcript::{AppendToTranscript, ProofTranscript};
use merlin::Transcript;
use serde::{Deserialize, Serialize};

/// Commitment to a polynomial
#[derive(Debug, Serialize, Deserialize, CanonicalSerialize, CanonicalDeserialize)]
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

/// Blinds for polynomial commitment
pub struct PolyCommitmentBlinds {
    pub blinds: Vec<Scalar>,
}

/// Generators for polynomial commitment
#[derive(Clone, Serialize, Deserialize)]
pub struct PolyCommitmentGens {
    pub gens: DotProductProofGens,
}

impl PolyCommitmentGens {
    pub fn new(num_vars: usize, label: &'static [u8]) -> Self {
        let (_left, right) = EqPolynomial::compute_factored_lens(num_vars);
        let gens = DotProductProofGens::new(right.pow2(), label);
        PolyCommitmentGens { gens }
    }
}

/// Polynomial evaluation proof
#[derive(Debug, Serialize, Deserialize, CanonicalSerialize, CanonicalDeserialize)]
pub struct PolyEvalProof {
    proof: DotProductProofLog,
}

impl PolyEvalProof {
    fn protocol_name() -> &'static [u8] {
        b"polynomial evaluation proof"
    }

    /// Prove evaluation of a polynomial at a point
    #[allow(non_snake_case)]
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
        let blinds = blinds_opt.unwrap_or(&default_blinds);

        assert_eq!(blinds.blinds.len(), L_size);

        let zero = Scalar::zero();
        let blind_Zr = blind_Zr_opt.unwrap_or(&zero);

        // compute L and R vectors
        let eq = EqPolynomial::new(r.to_vec());
        let (L, R) = eq.compute_factored_evals();
        assert_eq!(L.len(), L_size);
        assert_eq!(R.len(), R_size);

        // compute LZ and L*blinds
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

    /// Verify the polynomial evaluation proof
    #[allow(non_snake_case)]
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
}

/// Generators for R1CS sumcheck
#[derive(Serialize, Deserialize)]
pub struct R1CSSumcheckGens {
    pub gens_1: MultiCommitGens,
    pub gens_3: MultiCommitGens,
    pub gens_4: MultiCommitGens,
}

impl R1CSSumcheckGens {
    pub fn new(label: &'static [u8], gens_1_ref: &MultiCommitGens) -> Self {
        let gens_1 = gens_1_ref.clone();
        let gens_3 = MultiCommitGens::new(3, label);
        let gens_4 = MultiCommitGens::new(4, label);

        R1CSSumcheckGens {
            gens_1,
            gens_3,
            gens_4,
        }
    }
}

/// Generators for R1CS proof
#[derive(Serialize, Deserialize)]
pub struct R1CSGens {
    pub gens_sc: R1CSSumcheckGens,
    pub gens_pc: PolyCommitmentGens,
}

impl R1CSGens {
    pub fn new(label: &'static [u8], _num_cons: usize, num_vars: usize) -> Self {
        let num_poly_vars = num_vars.log_2();
        let gens_pc = PolyCommitmentGens::new(num_poly_vars, label);
        let gens_sc = R1CSSumcheckGens::new(label, &gens_pc.gens.gens_1);
        R1CSGens { gens_sc, gens_pc }
    }
}

/// R1CS Proof - proves satisfiability of an R1CS instance
#[derive(Serialize, Deserialize, Debug, CanonicalSerialize, CanonicalDeserialize)]
pub struct R1CSProof {
    comm_vars: PolyCommitment,
    sc_proof_phase1: ZKSumcheckInstanceProof,
    claims_phase2: (
        GroupElement,
        GroupElement,
        GroupElement,
        GroupElement,
    ),
    pok_claims_phase2: (KnowledgeProof, ProductProof),
    proof_eq_sc_phase1: EqualityProof,
    sc_proof_phase2: ZKSumcheckInstanceProof,
    comm_vars_at_ry: GroupElement,
    proof_eval_vars_at_ry: PolyEvalProof,
    proof_eq_sc_phase2: EqualityProof,
}

impl R1CSProof {
    fn protocol_name() -> &'static [u8] {
        b"R1CS proof"
    }

    /// Commit a polynomial and return commitment + blinds
    fn commit_poly(
        poly: &DensePolynomial,
        gens: &PolyCommitmentGens,
        random_tape: &mut RandomTape,
    ) -> (PolyCommitment, PolyCommitmentBlinds) {
        let n = poly.len();
        let ell = poly.get_num_vars();
        assert_eq!(n, ell.pow2());

        let (left_num_vars, right_num_vars) = EqPolynomial::compute_factored_lens(ell);
        let L_size = left_num_vars.pow2();
        let R_size = right_num_vars.pow2();
        assert_eq!(L_size * R_size, n);

        let blinds = PolyCommitmentBlinds {
            blinds: random_tape.random_vector(b"poly_blinds", L_size),
        };

        let Z = poly.vec();
        let C: Vec<GroupElement> = (0..L_size)
            .map(|i| {
                Z[R_size * i..R_size * (i + 1)]
                    .commit(&blinds.blinds[i], &gens.gens.gens_n)
            })
            .collect();

        (PolyCommitment { C }, blinds)
    }

    /// Prove R1CS satisfiability
    #[allow(non_snake_case)]
    pub fn prove(
        inst: &R1CSShape,
        vars: Vec<Scalar>,
        input: &[Scalar],
        gens: &R1CSGens,
        transcript: &mut Transcript,
        random_tape: &mut RandomTape,
    ) -> (R1CSProof, Vec<Scalar>, Vec<Scalar>) {
        let timer_prove = Timer::new("R1CSProof::prove");
        transcript.append_protocol_name(R1CSProof::protocol_name());

        // we require the number of inputs + 1 to be at most number of vars
        assert!(input.len() < vars.len());

        input.append_to_transcript(b"input", transcript);

        let timer_commit = Timer::new("polycommit");
        let (poly_vars, comm_vars, blinds_vars) = {
            let poly_vars = DensePolynomial::new(vars.clone());
            let (comm_vars, blinds_vars) = R1CSProof::commit_poly(&poly_vars, &gens.gens_pc, random_tape);
            comm_vars.append_to_transcript(b"poly_commitment", transcript);
            (poly_vars, comm_vars, blinds_vars)
        };
        timer_commit.stop();

        let timer_sc_proof_phase1 = Timer::new("prove_sc_phase_one");

        // Build z in Spartan format: [vars, 1, inputs, padding]
        let z = {
            let num_inputs = input.len();
            let num_vars = vars.len();
            let mut z = vars;
            z.push(Scalar::one()); // constant term
            z.extend(input);
            z.resize(2 * num_vars, Scalar::zero()); // pad with zeros
            z
        };

        // derive the verifier's challenge tau
        let (num_rounds_x, num_rounds_y) = (inst.get_num_cons().log_2(), z.len().log_2());
        let tau = transcript.challenge_vector(b"challenge_tau", num_rounds_x);
        
        // compute the initial evaluation table for R(tau, x)
        let mut poly_tau = DensePolynomial::new(EqPolynomial::new(tau).evals());
        let (mut poly_Az, mut poly_Bz, mut poly_Cz) =
            inst.multiply_vec(inst.get_num_cons(), z.len(), &z);

        let comb_func = |poly_A_comp: &Scalar,
                         poly_B_comp: &Scalar,
                         poly_C_comp: &Scalar,
                         poly_D_comp: &Scalar|
         -> Scalar { *poly_A_comp * (*poly_B_comp * *poly_C_comp - *poly_D_comp) };

        let (sc_proof_phase1, rx, _claims_phase1, blind_claim_postsc1) =
            ZKSumcheckInstanceProof::prove_cubic_with_additive_term(
                &Scalar::zero(), // claim is zero
                &Scalar::zero(), // blind for claim is also zero
                num_rounds_x,
                &mut poly_tau,
                &mut poly_Az,
                &mut poly_Bz,
                &mut poly_Cz,
                comb_func,
                &gens.gens_sc.gens_1,
                &gens.gens_sc.gens_4,
                transcript,
                random_tape,
            );
        assert_eq!(poly_tau.len(), 1);
        assert_eq!(poly_Az.len(), 1);
        assert_eq!(poly_Bz.len(), 1);
        assert_eq!(poly_Cz.len(), 1);
        timer_sc_proof_phase1.stop();

        let (tau_claim, Az_claim, Bz_claim, Cz_claim) =
            (&poly_tau[0], &poly_Az[0], &poly_Bz[0], &poly_Cz[0]);
        let (Az_blind, Bz_blind, Cz_blind, prod_Az_Bz_blind) = (
            random_tape.random_scalar(b"Az_blind"),
            random_tape.random_scalar(b"Bz_blind"),
            random_tape.random_scalar(b"Cz_blind"),
            random_tape.random_scalar(b"prod_Az_Bz_blind"),
        );

        let (pok_Cz_claim, comm_Cz_claim) = {
            KnowledgeProof::prove(
                &gens.gens_sc.gens_1,
                transcript,
                random_tape,
                Cz_claim,
                &Cz_blind,
            )
        };

        let (proof_prod, comm_Az_claim, comm_Bz_claim, comm_prod_Az_Bz_claims) = {
            let prod = *Az_claim * *Bz_claim;
            ProductProof::prove(
                &gens.gens_sc.gens_1,
                transcript,
                random_tape,
                Az_claim,
                &Az_blind,
                Bz_claim,
                &Bz_blind,
                &prod,
                &prod_Az_Bz_blind,
            )
        };

        comm_Az_claim.append_to_transcript(b"comm_Az_claim", transcript);
        comm_Bz_claim.append_to_transcript(b"comm_Bz_claim", transcript);
        comm_Cz_claim.append_to_transcript(b"comm_Cz_claim", transcript);
        comm_prod_Az_Bz_claims.append_to_transcript(b"comm_prod_Az_Bz_claims", transcript);

        // prove the final step of sum-check #1
        let taus_bound_rx = tau_claim;
        let blind_expected_claim_postsc1 = *taus_bound_rx * (prod_Az_Bz_blind - Cz_blind);
        let claim_post_phase1 = (*Az_claim * *Bz_claim - *Cz_claim) * *taus_bound_rx;
        let (proof_eq_sc_phase1, _C1, _C2) = EqualityProof::prove(
            &gens.gens_sc.gens_1,
            transcript,
            random_tape,
            &claim_post_phase1,
            &blind_expected_claim_postsc1,
            &claim_post_phase1,
            &blind_claim_postsc1,
        );

        let timer_sc_proof_phase2 = Timer::new("prove_sc_phase_two");
        // combine the three claims into a single claim
        let r_A = transcript.challenge_scalar(b"challenge_Az");
        let r_B = transcript.challenge_scalar(b"challenge_Bz");
        let r_C = transcript.challenge_scalar(b"challenge_Cz");
        let claim_phase2 = r_A * *Az_claim + r_B * *Bz_claim + r_C * *Cz_claim;
        let blind_claim_phase2 = r_A * Az_blind + r_B * Bz_blind + r_C * Cz_blind;

        let evals_ABC = {
            // compute the initial evaluation table for R(tau, x)
            let evals_rx = EqPolynomial::new(rx.clone()).evals();
            let (evals_A, evals_B, evals_C) =
                inst.compute_eval_table_sparse(inst.get_num_cons(), z.len(), &evals_rx);

            assert_eq!(evals_A.len(), evals_B.len());
            assert_eq!(evals_A.len(), evals_C.len());
            (0..evals_A.len())
                .map(|i| r_A * evals_A[i] + r_B * evals_B[i] + r_C * evals_C[i])
                .collect::<Vec<Scalar>>()
        };

        let comb_func2 =
            |poly_A_comp: &Scalar, poly_B_comp: &Scalar| -> Scalar { *poly_A_comp * *poly_B_comp };

        // another instance of the sum-check protocol
        let (sc_proof_phase2, ry, claims_phase2, blind_claim_postsc2) =
            ZKSumcheckInstanceProof::prove_quad(
                &claim_phase2,
                &blind_claim_phase2,
                num_rounds_y,
                &mut DensePolynomial::new(z),
                &mut DensePolynomial::new(evals_ABC),
                comb_func2,
                &gens.gens_sc.gens_1,
                &gens.gens_sc.gens_3,
                transcript,
                random_tape,
            );
        timer_sc_proof_phase2.stop();

        let timer_polyeval = Timer::new("polyeval");
        let eval_vars_at_ry = poly_vars.evaluate(&ry[1..]);
        let blind_eval = random_tape.random_scalar(b"blind_eval");
        let (proof_eval_vars_at_ry, comm_vars_at_ry) = PolyEvalProof::prove(
            &poly_vars,
            Some(&blinds_vars),
            &ry[1..],
            &eval_vars_at_ry,
            Some(&blind_eval),
            &gens.gens_pc,
            transcript,
            random_tape,
        );
        timer_polyeval.stop();

        // prove the final step of sum-check #2
        let blind_eval_Z_at_ry = (Scalar::one() - ry[0]) * blind_eval;
        let blind_expected_claim_postsc2 = claims_phase2[1] * blind_eval_Z_at_ry;
        let claim_post_phase2 = claims_phase2[0] * claims_phase2[1];
        let (proof_eq_sc_phase2, _C1, _C2) = EqualityProof::prove(
            &gens.gens_pc.gens.gens_1,
            transcript,
            random_tape,
            &claim_post_phase2,
            &blind_expected_claim_postsc2,
            &claim_post_phase2,
            &blind_claim_postsc2,
        );

        timer_prove.stop();

        (
            R1CSProof {
                comm_vars,
                sc_proof_phase1,
                claims_phase2: (
                    comm_Az_claim,
                    comm_Bz_claim,
                    comm_Cz_claim,
                    comm_prod_Az_Bz_claims,
                ),
                pok_claims_phase2: (pok_Cz_claim, proof_prod),
                proof_eq_sc_phase1,
                sc_proof_phase2,
                comm_vars_at_ry,
                proof_eval_vars_at_ry,
                proof_eq_sc_phase2,
            },
            rx,
            ry,
        )
    }

    /// Verify R1CS proof
    #[allow(non_snake_case)]
    pub fn verify(
        &self,
        num_vars: usize,
        num_cons: usize,
        input: &[Scalar],
        evals: &(Scalar, Scalar, Scalar),
        transcript: &mut Transcript,
        gens: &R1CSGens,
    ) -> Result<(Vec<Scalar>, Vec<Scalar>), ProofVerifyError> {
        transcript.append_protocol_name(R1CSProof::protocol_name());

        input.append_to_transcript(b"input", transcript);

        let n = num_vars;
        // add the commitment to the verifier's transcript
        self.comm_vars
            .append_to_transcript(b"poly_commitment", transcript);

        let (num_rounds_x, num_rounds_y) = (num_cons.log_2(), (2 * num_vars).log_2());

        // derive the verifier's challenge tau
        let tau = transcript.challenge_vector(b"challenge_tau", num_rounds_x);

        // verify the first sum-check instance
        let claim_phase1 = Scalar::zero()
            .commit(&Scalar::zero(), &gens.gens_sc.gens_1);
        let sc_result = self.sc_proof_phase1.verify(
            &claim_phase1,
            num_rounds_x,
            3,
            &gens.gens_sc.gens_1,
            &gens.gens_sc.gens_4,
            transcript,
        );
        let (comm_claim_post_phase1, rx) = match sc_result {
            Ok(v) => v,
            Err(e) => {
                eprintln!("sc_proof_phase1.verify failed: {:?}", e);
                return Err(e);
            }
        };
        
        // perform the intermediate sum-check test with claimed Az, Bz, and Cz
        let (comm_Az_claim, comm_Bz_claim, comm_Cz_claim, comm_prod_Az_Bz_claims) =
            &self.claims_phase2;
        let (pok_Cz_claim, proof_prod) = &self.pok_claims_phase2;

        pok_Cz_claim.verify(&gens.gens_sc.gens_1, transcript, comm_Cz_claim)?;
        proof_prod.verify(
            &gens.gens_sc.gens_1,
            transcript,
            comm_Az_claim,
            comm_Bz_claim,
            comm_prod_Az_Bz_claims,
        )?;

        comm_Az_claim.append_to_transcript(b"comm_Az_claim", transcript);
        comm_Bz_claim.append_to_transcript(b"comm_Bz_claim", transcript);
        comm_Cz_claim.append_to_transcript(b"comm_Cz_claim", transcript);
        comm_prod_Az_Bz_claims.append_to_transcript(b"comm_prod_Az_Bz_claims", transcript);

        let taus_bound_rx: Scalar = (0..rx.len())
            .map(|i| rx[i] * tau[i] + (Scalar::one() - rx[i]) * (Scalar::one() - tau[i]))
            .product();
        let expected_claim_post_phase1 = taus_bound_rx
            * (*comm_prod_Az_Bz_claims + (Scalar::zero() - Scalar::one()) * *comm_Cz_claim);

        // verify proof that expected_claim_post_phase1 == claim_post_phase1
        self.proof_eq_sc_phase1.verify(
            &gens.gens_sc.gens_1,
            transcript,
            &expected_claim_post_phase1,
            &comm_claim_post_phase1,
        )?;

        // derive three public challenges and then derive a joint claim
        let r_A = transcript.challenge_scalar(b"challenge_Az");
        let r_B = transcript.challenge_scalar(b"challenge_Bz");
        let r_C = transcript.challenge_scalar(b"challenge_Cz");

        // r_A * comm_Az_claim + r_B * comm_Bz_claim + r_C * comm_Cz_claim
        let comm_claim_phase2 = GroupElement::vartime_multiscalar_mul(
            &[r_A, r_B, r_C],
            &[*comm_Az_claim, *comm_Bz_claim, *comm_Cz_claim],
        );

        // verify the joint claim with a sum-check protocol
        let sc2_result = self.sc_proof_phase2.verify(
            &comm_claim_phase2,
            num_rounds_y,
            2,
            &gens.gens_sc.gens_1,
            &gens.gens_sc.gens_3,
            transcript,
        );
        let (comm_claim_post_phase2, ry) = match sc2_result {
            Ok(v) => v,
            Err(e) => {
                eprintln!("sc_proof_phase2.verify failed: {:?}", e);
                return Err(e);
            }
        };

        // verify Z(ry) proof against the initial commitment
        let eval_result = self.proof_eval_vars_at_ry.verify(
            &gens.gens_pc,
            transcript,
            &ry[1..],
            &self.comm_vars_at_ry,
            &self.comm_vars,
        );
        if let Err(e) = &eval_result {
            eprintln!("proof_eval_vars_at_ry.verify failed: {:?}", e);
            return Err(e.clone());
        }

        // Build the sparse input polynomial evaluation
        let poly_input_eval = {
            use crate::sparse_mlpoly::{SparseMatEntry, SparseMatPolynomial};
            use crate::hyrax::EqPolynomial;
            
            // constant term
            let mut entries = vec![SparseMatEntry::new(0, 0, Scalar::one())];
            // remaining inputs
            for (i, inp) in input.iter().enumerate() {
                entries.push(SparseMatEntry::new(0, i + 1, *inp));
            }
            let sparse_poly = SparseMatPolynomial::new(1, n.log_2(), entries);
            // Compute full equality polynomial evaluations for proper table lookup
            let ry_evals = EqPolynomial::new(ry[1..].to_vec()).evals();
            sparse_poly.evaluate_with_tables(&[Scalar::one()], &ry_evals)
        };

        // compute commitment to eval_Z_at_ry
        let comm_eval_Z_at_ry = GroupElement::vartime_multiscalar_mul(
            &[Scalar::one() - ry[0], ry[0]],
            &[
                self.comm_vars_at_ry,
                poly_input_eval.commit(&Scalar::zero(), &gens.gens_pc.gens.gens_1),
            ],
        );

        // perform the final check in the second sum-check protocol
        let (eval_A_r, eval_B_r, eval_C_r) = evals;
        let expected_claim_post_phase2 =
            (r_A * *eval_A_r + r_B * *eval_B_r + r_C * *eval_C_r) * comm_eval_Z_at_ry;
            
        // verify proof that expected_claim_post_phase1 == claim_post_phase1
        self.proof_eq_sc_phase2.verify(
            &gens.gens_sc.gens_1,
            transcript,
            &expected_claim_post_phase2,
            &comm_claim_post_phase2,
        )?;

        Ok((rx, ry))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::r1cs::R1CSShape;

    fn produce_tiny_r1cs() -> (R1CSShape, Vec<Scalar>, Vec<Scalar>) {
        // A simple constraint: x * x = x (satisfied by x=0 or x=1)
        let num_cons = 4; // power of 2
        let num_vars = 4; // power of 2
        let num_inputs = 1;

        let one = Scalar::one();
        
        // constraint 0: x * x = x
        let A = vec![(0, 0, one)];
        let B = vec![(0, 0, one)];
        let C = vec![(0, 0, one)];

        let inst = R1CSShape::new(num_cons, num_vars, num_inputs, &A, &B, &C);

        // x = 1 satisfies x * x = x
        let mut vars = vec![Scalar::zero(); num_vars];
        vars[0] = Scalar::one();

        let inputs = vec![Scalar::one()];

        (inst, vars, inputs)
    }

    #[test]
    fn test_r1cs_proof_tiny() {
        let (inst, vars, input) = produce_tiny_r1cs();
        
        let gens = R1CSGens::new(b"test-gens", inst.get_num_cons(), inst.get_num_vars());

        let mut random_tape = RandomTape::new(b"proof");
        let mut prover_transcript = Transcript::new(b"example");
        let (proof, rx, ry) = R1CSProof::prove(
            &inst,
            vars,
            &input,
            &gens,
            &mut prover_transcript,
            &mut random_tape,
        );

        let inst_evals = inst.evaluate(&rx, &ry);

        let mut verifier_transcript = Transcript::new(b"example");
        assert!(proof
            .verify(
                inst.get_num_vars(),
                inst.get_num_cons(),
                &input,
                &inst_evals,
                &mut verifier_transcript,
                &gens,
            )
            .is_ok());
    }
}

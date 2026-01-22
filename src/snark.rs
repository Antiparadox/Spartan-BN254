//! SNARK and NIZK high-level interfaces
//! Port of Spartan's top-level SNARK and NIZK from lib.rs

use ark_serialize::{CanonicalSerialize, CanonicalDeserialize};
use crate::errors::{ProofVerifyError, R1CSError};
use crate::r1cs::{R1CSCommitment, R1CSCommitmentGens, R1CSDecommitment, R1CSEvalProof, R1CSShape};
use crate::r1csproof::{R1CSGens, R1CSProof};
use crate::random::RandomTape;
use crate::scalar::Scalar;
use crate::timer::Timer;
use crate::transcript::{AppendToTranscript, ProofTranscript};
use merlin::Transcript;
#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};
use std::cmp::max;

/// Assignment of values to variables
#[derive(Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Assignment {
    pub assignment: Vec<Scalar>,
}

impl Assignment {
    /// Create from byte arrays (Spartan's original format)
    pub fn new(assignment: &[[u8; 32]]) -> Result<Assignment, R1CSError> {
        let mut vec_scalar: Vec<Scalar> = Vec::new();
        for v in assignment {
            let val = Scalar::from_bytes(v).ok_or(R1CSError::InvalidScalar)?;
            vec_scalar.push(val);
        }
        Ok(Assignment {
            assignment: vec_scalar,
        })
    }

    /// Create from scalars directly
    pub fn from_scalars(scalars: Vec<Scalar>) -> Self {
        Assignment { assignment: scalars }
    }

    /// Pad to specified length
    fn pad(&self, len: usize) -> VarsAssignment {
        assert!(len > self.assignment.len());
        let mut padded = self.assignment.clone();
        padded.resize(len, Scalar::zero());
        VarsAssignment {
            assignment: padded,
        }
    }
}

/// Variable assignment
pub type VarsAssignment = Assignment;
/// Input assignment
pub type InputsAssignment = Assignment;

/// Instance wrapper around R1CSShape with digest
pub struct Instance {
    pub inst: R1CSShape,
    pub digest: Vec<u8>,
}

impl Instance {
    /// Create from raw triplets (original Spartan API)
    pub fn new(
        num_cons: usize,
        num_vars: usize,
        num_inputs: usize,
        A: &[(usize, usize, [u8; 32])],
        B: &[(usize, usize, [u8; 32])],
        C: &[(usize, usize, [u8; 32])],
    ) -> Result<Instance, R1CSError> {
        let (num_vars_padded, num_cons_padded) = {
            let num_vars_padded = {
                let mut n = max(num_vars, num_inputs + 1);
                if n.next_power_of_two() != n {
                    n = n.next_power_of_two();
                }
                n
            };
            let num_cons_padded = {
                let mut n = max(num_cons, 2);
                if n.next_power_of_two() != n {
                    n = n.next_power_of_two();
                }
                n
            };
            (num_vars_padded, num_cons_padded)
        };

        let convert = |tups: &[(usize, usize, [u8; 32])]| -> Result<Vec<(usize, usize, Scalar)>, R1CSError> {
            let mut mat = Vec::new();
            for &(row, col, val_bytes) in tups {
                if row >= num_cons {
                    return Err(R1CSError::InvalidIndex);
                }
                if col >= num_vars + 1 + num_inputs {
                    return Err(R1CSError::InvalidIndex);
                }
                let val = Scalar::from_bytes(&val_bytes).ok_or(R1CSError::InvalidScalar)?;
                let adjusted_col = if col >= num_vars {
                    col + num_vars_padded - num_vars
                } else {
                    col
                };
                mat.push((row, adjusted_col, val));
            }
            Ok(mat)
        };

        let A_scalar = convert(A)?;
        let B_scalar = convert(B)?;
        let C_scalar = convert(C)?;

        let inst = R1CSShape::new(
            num_cons_padded,
            num_vars_padded,
            num_inputs,
            &A_scalar,
            &B_scalar,
            &C_scalar,
        );

        let digest = inst.get_digest();

        Ok(Instance { inst, digest })
    }

    /// Create from R1CSShape directly
    pub fn from_shape(shape: R1CSShape) -> Self {
        let digest = shape.get_digest();
        Instance { inst: shape, digest }
    }

    /// Check if instance is satisfied
    pub fn is_sat(
        &self,
        vars: &VarsAssignment,
        inputs: &InputsAssignment,
    ) -> Result<bool, R1CSError> {
        if vars.assignment.len() > self.inst.get_num_vars() {
            return Err(R1CSError::InvalidNumberOfInputs);
        }
        if inputs.assignment.len() != self.inst.get_num_inputs() {
            return Err(R1CSError::InvalidNumberOfInputs);
        }

        let padded_vars = {
            let num_padded = self.inst.get_num_vars();
            if num_padded > vars.assignment.len() {
                vars.pad(num_padded)
            } else {
                vars.clone()
            }
        };

        Ok(self.inst.is_sat(&padded_vars.assignment, &inputs.assignment))
    }
}

/// Generators for NIZK proofs (no preprocessing)
pub struct NIZKGens {
    pub gens_r1cs_sat: R1CSGens,
}

impl NIZKGens {
    /// Create generators for the given R1CS dimensions
    pub fn new(num_cons: usize, num_vars: usize, num_inputs: usize) -> Self {
        let num_vars_padded = {
            let mut n = max(num_vars, num_inputs + 1);
            if n.next_power_of_two() != n {
                n = n.next_power_of_two();
            }
            n
        };

        let gens_r1cs_sat = R1CSGens::new(b"gens_r1cs_sat", num_cons, num_vars_padded);
        NIZKGens { gens_r1cs_sat }
    }
}

/// NIZK proof - non-interactive zero-knowledge proof
/// 
/// Note: For NIZK (unlike SNARK), the `r` vectors ARE needed because:
/// 1. Verifier must evaluate R1CS at (rx, ry) BEFORE sumcheck verification
/// 2. These evaluations are inputs to the sumcheck
/// 3. Sumcheck then produces (rx', ry') which must match
#[derive(Debug, CanonicalSerialize, CanonicalDeserialize)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct NIZK {
    r1cs_sat_proof: R1CSProof,
    r: (Vec<Scalar>, Vec<Scalar>),
}

impl NIZK {
    fn protocol_name() -> &'static [u8] {
        b"Spartan NIZK proof"
    }

    /// Produce a NIZK proof of R1CS satisfiability
    pub fn prove(
        inst: &Instance,
        vars: VarsAssignment,
        input: &InputsAssignment,
        gens: &NIZKGens,
        transcript: &mut Transcript,
    ) -> Self {
        let timer_prove = Timer::new("NIZK::prove");
        let mut random_tape = RandomTape::new(b"proof");

        transcript.append_protocol_name(NIZK::protocol_name());
        transcript.append_message(b"R1CSShapeDigest", &inst.digest);

        let (r1cs_sat_proof, rx, ry) = {
            let padded_vars = {
                let num_padded = inst.inst.get_num_vars();
                if num_padded > vars.assignment.len() {
                    vars.pad(num_padded)
                } else {
                    vars
                }
            };

            R1CSProof::prove(
                &inst.inst,
                padded_vars.assignment,
                &input.assignment,
                &gens.gens_r1cs_sat,
                transcript,
                &mut random_tape,
            )
        };

        timer_prove.stop();
        NIZK {
            r1cs_sat_proof,
            r: (rx, ry),
        }
    }

    /// Verify a NIZK proof
    pub fn verify(
        &self,
        inst: &Instance,
        input: &InputsAssignment,
        transcript: &mut Transcript,
        gens: &NIZKGens,
    ) -> Result<(), ProofVerifyError> {
        let timer_verify = Timer::new("NIZK::verify");

        transcript.append_protocol_name(NIZK::protocol_name());
        transcript.append_message(b"R1CSShapeDigest", &inst.digest);

        // Compute evaluations of A, B, C at r = (rx, ry)
        let timer_eval = Timer::new("eval_sparse_polys");
        let (claimed_rx, claimed_ry) = &self.r;
        let inst_evals = inst.inst.evaluate(claimed_rx, claimed_ry);
        timer_eval.stop();

        let timer_sat_proof = Timer::new("verify_sat_proof");
        assert_eq!(input.assignment.len(), inst.inst.get_num_inputs());
        let verify_result = self.r1cs_sat_proof.verify(
            inst.inst.get_num_vars(),
            inst.inst.get_num_cons(),
            &input.assignment,
            &inst_evals,
            transcript,
            &gens.gens_r1cs_sat,
        );
        let (rx, ry) = match verify_result {
            Ok(v) => v,
            Err(e) => {
                eprintln!("R1CS sat proof verification failed: {:?}", e);
                return Err(e);
            }
        };

        // verify that claimed rx and ry are correct
        assert_eq!(rx, *claimed_rx);
        assert_eq!(ry, *claimed_ry);
        timer_sat_proof.stop();
        timer_verify.stop();

        Ok(())
    }
}

/// Generators for SNARK proofs (with preprocessing)
pub struct SNARKGens {
    pub gens_r1cs_sat: R1CSGens,
    pub gens_r1cs_eval: R1CSCommitmentGens,
}

impl SNARKGens {
    /// Create generators for the given R1CS dimensions
    /// 
    /// # Arguments
    /// * `num_cons` - Number of constraints (will be padded to power of 2)
    /// * `num_vars` - Number of variables (will be padded to power of 2)
    /// * `num_inputs` - Number of public inputs
    /// * `num_nz_entries` - Maximum number of non-zero entries across A, B, C matrices
    ///                      (i.e., max(nnz_a, nnz_b, nnz_c))
    pub fn new(num_cons: usize, num_vars: usize, num_inputs: usize, num_nz_entries: usize) -> Self {
        let num_vars_padded = {
            let mut n = max(num_vars, num_inputs + 1);
            if n.next_power_of_two() != n {
                n = n.next_power_of_two();
            }
            n
        };
        let num_cons_padded = {
            let mut n = max(num_cons, 2);
            if n.next_power_of_two() != n {
                n = n.next_power_of_two();
            }
            n
        };

        let gens_r1cs_sat = R1CSGens::new(b"gens_r1cs_sat", num_cons_padded, num_vars_padded);
        let gens_r1cs_eval = R1CSCommitmentGens::new(b"gens_r1cs_eval", num_cons_padded, num_vars_padded, num_nz_entries);

        SNARKGens {
            gens_r1cs_sat,
            gens_r1cs_eval,
        }
    }
}

/// SNARK proof - succinct verification, verifier doesn't need full R1CS
/// 
/// Note: The random points (rx, ry) are NOT stored in the proof because
/// the verifier recomputes them from the Fiat-Shamir transcript.
/// 
/// Serialization is OPTIONAL - only needed if you're sending the proof
/// over a network or saving to disk. For in-process verification,
/// just pass the SNARK struct directly.
/// 
/// Uses CanonicalSerialize for cross-verification compatibility with arkworks-spartan.
#[derive(Debug, CanonicalSerialize, CanonicalDeserialize)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SNARK {
    r1cs_sat_proof: R1CSProof,
    inst_evals: (Scalar, Scalar, Scalar),
    r1cs_eval_proof: R1CSEvalProof,
}

impl SNARK {
    fn protocol_name() -> &'static [u8] {
        b"Spartan SNARK proof"
    }

    /// Encode the R1CS instance (preprocessing)
    pub fn encode(
        inst: &Instance,
        gens: &SNARKGens,
    ) -> (R1CSCommitment, R1CSDecommitment) {
        let timer_encode = Timer::new("SNARK::encode");
        let (comm, decomm) = inst.inst.commit(&gens.gens_r1cs_eval);
        timer_encode.stop();
        (comm, decomm)
    }

    /// Produce a SNARK proof of R1CS satisfiability
    pub fn prove(
        inst: &Instance,
        comm: &R1CSCommitment,
        decomm: &R1CSDecommitment,
        vars: VarsAssignment,
        input: &InputsAssignment,
        gens: &SNARKGens,
        transcript: &mut Transcript,
    ) -> Self {
        let timer_prove = Timer::new("SNARK::prove");
        let mut random_tape = RandomTape::new(b"snark_proof");

        transcript.append_protocol_name(SNARK::protocol_name());
        // Use the commitment append_to_transcript for consistency with arkworks-spartan
        comm.append_to_transcript(b"comm", transcript);

        let (r1cs_sat_proof, rx, ry) = {
            let padded_vars = {
                let num_padded = inst.inst.get_num_vars();
                if num_padded > vars.assignment.len() {
                    vars.pad(num_padded)
                } else {
                    vars
                }
            };

            R1CSProof::prove(
                &inst.inst,
                padded_vars.assignment,
                &input.assignment,
                &gens.gens_r1cs_sat,
                transcript,
                &mut random_tape,
            )
        };

        // Compute evaluations of A, B, C at (rx, ry)
        let inst_evals = inst.inst.evaluate(&rx, &ry);

        // Produce proof of correct evaluations
        let r1cs_eval_proof = R1CSEvalProof::prove(
            decomm,
            &rx,
            &ry,
            &inst_evals,
            &gens.gens_r1cs_eval,
            transcript,
            &mut random_tape,
        );

        timer_prove.stop();
        SNARK {
            r1cs_sat_proof,
            inst_evals,
            r1cs_eval_proof,
        }
    }

    /// Verify a SNARK proof (succinct - doesn't need full R1CS)
    pub fn verify(
        &self,
        comm: &R1CSCommitment,
        input: &InputsAssignment,
        transcript: &mut Transcript,
        gens: &SNARKGens,
    ) -> Result<(), ProofVerifyError> {
        let timer_verify = Timer::new("SNARK::verify");

        transcript.append_protocol_name(SNARK::protocol_name());
        // Use the commitment append_to_transcript for consistency with arkworks-spartan
        comm.append_to_transcript(b"comm", transcript);

        let timer_sat_proof = Timer::new("verify_sat_proof");
        assert_eq!(input.assignment.len(), comm.num_inputs);
        
        // Verify R1CS satisfiability proof - this derives rx, ry from the transcript
        let (rx, ry) = self.r1cs_sat_proof.verify(
            comm.num_vars,
            comm.num_cons,
            &input.assignment,
            &self.inst_evals,
            transcript,
            &gens.gens_r1cs_sat,
        )?;
        timer_sat_proof.stop();

        // Verify the R1CS evaluation proof
        let timer_eval_proof = Timer::new("verify_eval_proof");
        self.r1cs_eval_proof.verify(
            comm,
            &rx,
            &ry,
            &self.inst_evals,
            &gens.gens_r1cs_eval,
            transcript,
        )?;
        timer_eval_proof.stop();

        timer_verify.stop();
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_nizk_simple() {
        // Create a simple R1CS: x * x = x (satisfied by x=0 or x=1)
        let num_cons = 4;
        let num_vars = 4;
        let num_inputs = 1;

        let one = Scalar::one();
        let A = vec![(0, 0, one)];
        let B = vec![(0, 0, one)];
        let C = vec![(0, 0, one)];

        let shape = R1CSShape::new(num_cons, num_vars, num_inputs, &A, &B, &C);
        let inst = Instance::from_shape(shape);

        // x = 1 satisfies x * x = x
        let mut vars = vec![Scalar::zero(); num_vars];
        vars[0] = Scalar::one();
        let vars_assignment = VarsAssignment::from_scalars(vars);
        let inputs = InputsAssignment::from_scalars(vec![Scalar::one()]);

        // Check satisfaction
        assert!(inst.is_sat(&vars_assignment, &inputs).unwrap());

        // Create generators and proof
        let gens = NIZKGens::new(num_cons, num_vars, num_inputs);

        let mut prover_transcript = Transcript::new(b"nizk_test");
        let proof = NIZK::prove(&inst, vars_assignment.clone(), &inputs, &gens, &mut prover_transcript);

        // Verify
        let mut verifier_transcript = Transcript::new(b"nizk_test");
        let result = proof.verify(&inst, &inputs, &mut verifier_transcript, &gens);
        if let Err(e) = &result {
            eprintln!("Verification failed with error: {:?}", e);
        }
        assert!(result.is_ok());
    }

    #[test]
    fn test_snark_simple() {
        // Create a simple R1CS: x * x = x (satisfied by x=0 or x=1)
        let num_cons = 4;
        let num_vars = 4;
        let num_inputs = 1;

        let one = Scalar::one();
        let A = vec![(0, 0, one)];
        let B = vec![(0, 0, one)];
        let C = vec![(0, 0, one)];

        let shape = R1CSShape::new(num_cons, num_vars, num_inputs, &A, &B, &C);
        let inst = Instance::from_shape(shape);

        // x = 1 satisfies x * x = x
        let mut vars = vec![Scalar::zero(); num_vars];
        vars[0] = Scalar::one();
        let vars_assignment = VarsAssignment::from_scalars(vars);
        let inputs = InputsAssignment::from_scalars(vec![Scalar::one()]);

        // Check satisfaction
        assert!(inst.is_sat(&vars_assignment, &inputs).unwrap());

        // Create generators
        let gens = SNARKGens::new(num_cons, num_vars, num_inputs);

        // Encode (preprocess)
        let (comm, decomm) = SNARK::encode(&inst, &gens);

        // Prove
        let mut prover_transcript = Transcript::new(b"snark_test");
        let proof = SNARK::prove(&inst, &comm, &decomm, vars_assignment.clone(), &inputs, &gens, &mut prover_transcript);

        // Verify (succinct - uses only commitment, not full R1CS)
        let mut verifier_transcript = Transcript::new(b"snark_test");
        let result = proof.verify(&comm, &inputs, &mut verifier_transcript, &gens);
        if let Err(e) = &result {
            eprintln!("SNARK Verification failed with error: {:?}", e);
        }
        assert!(result.is_ok());
    }
}

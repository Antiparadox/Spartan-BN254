//! Error types for Spartan

use thiserror::Error;

#[derive(Error, Debug, Clone, PartialEq, Eq)]
pub enum R1CSError {
    #[error("Invalid index in R1CS matrix")]
    InvalidIndex,
    #[error("Invalid scalar value")]
    InvalidScalar,
    #[error("Invalid number of inputs")]
    InvalidNumberOfInputs,
    #[error("Constraint system is not satisfiable")]
    NotSatisfiable,
    #[error("Dimensions must be powers of 2")]
    InvalidDimensions,
}

#[derive(Error, Debug, Clone)]
pub enum ProofVerifyError {
    #[error("Proof verification failed: {0}")]
    VerificationFailed(String),
    #[error("Invalid proof format")]
    InvalidProof,
    #[error("Commitment verification failed")]
    CommitmentError,
    #[error("Sumcheck verification failed at round {0}")]
    SumcheckError(usize),
    #[error("Internal error")]
    InternalError,
}

//! Random tape for prover randomness
//! Port of Spartan's random.rs

use crate::scalar::Scalar;
use crate::transcript::ProofTranscript;
use merlin::Transcript;
use rand::rngs::OsRng;

/// Random tape for generating prover randomness
pub struct RandomTape {
    tape: Transcript,
}

impl RandomTape {
    pub fn new(name: &'static [u8]) -> Self {
        let tape = {
            let mut rng = OsRng;
            let mut tape = Transcript::new(name);
            tape.append_scalar(b"init_randomness", &Scalar::random(&mut rng));
            tape
        };
        Self { tape }
    }

    pub fn random_scalar(&mut self, label: &'static [u8]) -> Scalar {
        self.tape.challenge_scalar(label)
    }

    pub fn random_vector(&mut self, label: &'static [u8], len: usize) -> Vec<Scalar> {
        self.tape.challenge_scalars(label, len)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_random_tape() {
        let mut tape = RandomTape::new(b"test");
        let s1 = tape.random_scalar(b"r1");
        let s2 = tape.random_scalar(b"r2");
        // Different labels should give different scalars
        assert_ne!(s1, s2);
    }

    #[test]
    fn test_random_vector() {
        let mut tape = RandomTape::new(b"test");
        let v = tape.random_vector(b"vec", 5);
        assert_eq!(v.len(), 5);
    }
}

//! Fiat-Shamir transcript using Merlin

use crate::group::CompressedGroup;
use crate::scalar::Scalar;
use merlin::Transcript;

/// Trait for appending data to a transcript
pub trait AppendToTranscript {
    fn append_to_transcript(&self, label: &'static [u8], transcript: &mut Transcript);
}

/// Extension trait for Transcript to generate Scalars
pub trait ProofTranscript {
    /// Append a protocol name to the transcript
    fn append_protocol_name(&mut self, protocol_name: &'static [u8]);
    
    /// Append a scalar to the transcript
    fn append_scalar(&mut self, label: &'static [u8], scalar: &Scalar);
    
    /// Append multiple scalars to the transcript
    fn append_scalars(&mut self, label: &'static [u8], scalars: &[Scalar]);
    
    /// Append a point (compressed group element) to the transcript
    fn append_point(&mut self, label: &'static [u8], point: &CompressedGroup);
    
    /// Get a challenge scalar from the transcript
    fn challenge_scalar(&mut self, label: &'static [u8]) -> Scalar;
    
    /// Get multiple challenge scalars from the transcript
    fn challenge_scalars(&mut self, label: &'static [u8], n: usize) -> Vec<Scalar>;
    
    /// Get a vector of challenges
    fn challenge_vector(&mut self, label: &'static [u8], n: usize) -> Vec<Scalar>;
}

impl ProofTranscript for Transcript {
    fn append_protocol_name(&mut self, protocol_name: &'static [u8]) {
        self.append_message(b"protocol-name", protocol_name);
    }

    fn append_scalar(&mut self, label: &'static [u8], scalar: &Scalar) {
        self.append_message(label, &scalar.to_bytes());
    }

    fn append_scalars(&mut self, label: &'static [u8], scalars: &[Scalar]) {
        for scalar in scalars {
            self.append_message(label, &scalar.to_bytes());
        }
    }

    fn append_point(&mut self, label: &'static [u8], point: &CompressedGroup) {
        self.append_message(label, point.as_bytes());
    }

    fn challenge_scalar(&mut self, label: &'static [u8]) -> Scalar {
        let mut buf = [0u8; 64];
        self.challenge_bytes(label, &mut buf);
        
        // Convert 64 bytes to a scalar by interpreting as a big integer and reducing mod p
        use ark_ff::PrimeField;
        use ark_bn254::Fr;
        
        // Reduce modulo the field order
        let scalar = Fr::from_le_bytes_mod_order(&buf);
        Scalar(scalar)
    }

    fn challenge_scalars(&mut self, label: &'static [u8], n: usize) -> Vec<Scalar> {
        (0..n).map(|_| self.challenge_scalar(label)).collect()
    }

    fn challenge_vector(&mut self, label: &'static [u8], n: usize) -> Vec<Scalar> {
        self.challenge_scalars(label, n)
    }
}

impl AppendToTranscript for Scalar {
    fn append_to_transcript(&self, label: &'static [u8], transcript: &mut Transcript) {
        transcript.append_scalar(label, self);
    }
}

impl AppendToTranscript for [Scalar] {
    fn append_to_transcript(&self, label: &'static [u8], transcript: &mut Transcript) {
        transcript.append_scalars(label, self);
    }
}

impl AppendToTranscript for Vec<Scalar> {
    fn append_to_transcript(&self, label: &'static [u8], transcript: &mut Transcript) {
        transcript.append_scalars(label, self);
    }
}

impl AppendToTranscript for CompressedGroup {
    fn append_to_transcript(&self, label: &'static [u8], transcript: &mut Transcript) {
        transcript.append_point(label, self);
    }
}

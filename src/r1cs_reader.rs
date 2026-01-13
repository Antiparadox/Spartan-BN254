//! Reader for circom .r1cs binary format
//! Parses R1CS files produced by circom compiler

use crate::scalar::Scalar;
use ark_bn254::Fr;
use ark_ff::PrimeField;
use byteorder::{LittleEndian, ReadBytesExt};
use std::io::{Read, Seek, SeekFrom};
use thiserror::Error;

#[derive(Error, Debug)]
pub enum R1CSError {
    #[error("Invalid magic number")]
    InvalidMagic,
    #[error("Unsupported version: {0}")]
    UnsupportedVersion(u32),
    #[error("IO error: {0}")]
    IoError(#[from] std::io::Error),
    #[error("Invalid field size: expected 32, got {0}")]
    InvalidFieldSize(u32),
    #[error("Section not found: {0}")]
    SectionNotFound(u32),
}

/// Parsed R1CS constraint system
#[derive(Debug, Clone)]
pub struct R1CS {
    /// Number of constraints
    pub num_constraints: usize,
    /// Number of variables (wires)
    pub num_variables: usize,
    /// Number of public inputs (includes output wires)
    pub num_pub_inputs: usize,
    /// Number of private inputs
    pub num_prv_inputs: usize,
    /// Number of labels
    pub num_labels: u64,
    /// Matrix A as (row, col, value) triplets
    pub a: Vec<(usize, usize, Scalar)>,
    /// Matrix B as (row, col, value) triplets
    pub b: Vec<(usize, usize, Scalar)>,
    /// Matrix C as (row, col, value) triplets
    pub c: Vec<(usize, usize, Scalar)>,
}

impl R1CS {
    /// Read R1CS from a file path
    pub fn from_file(path: &str) -> Result<Self, R1CSError> {
        let file = std::fs::File::open(path)?;
        let mut reader = std::io::BufReader::new(file);
        Self::from_reader(&mut reader)
    }

    /// Read R1CS from any reader
    pub fn from_reader<R: Read + Seek>(reader: &mut R) -> Result<Self, R1CSError> {
        // Check magic number: "r1cs"
        let mut magic = [0u8; 4];
        reader.read_exact(&mut magic)?;
        if &magic != b"r1cs" {
            return Err(R1CSError::InvalidMagic);
        }

        // Version (should be 1)
        let version = reader.read_u32::<LittleEndian>()?;
        if version != 1 {
            return Err(R1CSError::UnsupportedVersion(version));
        }

        // Number of sections
        let num_sections = reader.read_u32::<LittleEndian>()?;

        // Find and parse header section (type 1)
        let (num_constraints, num_variables, num_pub_inputs, num_prv_inputs, num_labels, field_size) = 
            Self::read_header_section(reader, num_sections)?;

        if field_size != 32 {
            return Err(R1CSError::InvalidFieldSize(field_size));
        }

        // Find and parse constraints section (type 2)
        let (a, b, c) = Self::read_constraints_section(reader, num_sections, num_constraints, field_size)?;

        Ok(R1CS {
            num_constraints,
            num_variables,
            num_pub_inputs,
            num_prv_inputs,
            num_labels,
            a,
            b,
            c,
        })
    }

    fn read_header_section<R: Read + Seek>(
        reader: &mut R,
        num_sections: u32,
    ) -> Result<(usize, usize, usize, usize, u64, u32), R1CSError> {
        reader.seek(SeekFrom::Start(12))?; // After magic + version + num_sections

        for _ in 0..num_sections {
            let section_type = reader.read_u32::<LittleEndian>()?;
            let section_size = reader.read_u64::<LittleEndian>()?;

            if section_type == 1 {
                // Header section
                let field_size = reader.read_u32::<LittleEndian>()?;
                
                // Skip the prime field bytes
                reader.seek(SeekFrom::Current(field_size as i64))?;

                let num_variables = reader.read_u32::<LittleEndian>()? as usize;
                let num_pub_outputs = reader.read_u32::<LittleEndian>()? as usize;
                let num_pub_inputs = reader.read_u32::<LittleEndian>()? as usize;
                let num_prv_inputs = reader.read_u32::<LittleEndian>()? as usize;
                let num_labels = reader.read_u64::<LittleEndian>()?;
                let num_constraints = reader.read_u32::<LittleEndian>()? as usize;

                // Total public inputs includes outputs
                let total_pub = num_pub_outputs + num_pub_inputs;

                return Ok((num_constraints, num_variables, total_pub, num_prv_inputs, num_labels, field_size));
            } else {
                reader.seek(SeekFrom::Current(section_size as i64))?;
            }
        }

        Err(R1CSError::SectionNotFound(1))
    }

    fn read_constraints_section<R: Read + Seek>(
        reader: &mut R,
        num_sections: u32,
        num_constraints: usize,
        field_size: u32,
    ) -> Result<(Vec<(usize, usize, Scalar)>, Vec<(usize, usize, Scalar)>, Vec<(usize, usize, Scalar)>), R1CSError> {
        reader.seek(SeekFrom::Start(12))?; // After magic + version + num_sections

        for _ in 0..num_sections {
            let section_type = reader.read_u32::<LittleEndian>()?;
            let section_size = reader.read_u64::<LittleEndian>()?;

            if section_type == 2 {
                // Constraints section
                let mut a = Vec::new();
                let mut b = Vec::new();
                let mut c = Vec::new();

                for constraint_idx in 0..num_constraints {
                    // Read A entries
                    let num_a = reader.read_u32::<LittleEndian>()? as usize;
                    for _ in 0..num_a {
                        let col = reader.read_u32::<LittleEndian>()? as usize;
                        let mut val_bytes = vec![0u8; field_size as usize];
                        reader.read_exact(&mut val_bytes)?;
                        if let Some(scalar) = bytes_to_scalar(&val_bytes) {
                            a.push((constraint_idx, col, scalar));
                        }
                    }

                    // Read B entries
                    let num_b = reader.read_u32::<LittleEndian>()? as usize;
                    for _ in 0..num_b {
                        let col = reader.read_u32::<LittleEndian>()? as usize;
                        let mut val_bytes = vec![0u8; field_size as usize];
                        reader.read_exact(&mut val_bytes)?;
                        if let Some(scalar) = bytes_to_scalar(&val_bytes) {
                            b.push((constraint_idx, col, scalar));
                        }
                    }

                    // Read C entries
                    let num_c = reader.read_u32::<LittleEndian>()? as usize;
                    for _ in 0..num_c {
                        let col = reader.read_u32::<LittleEndian>()? as usize;
                        let mut val_bytes = vec![0u8; field_size as usize];
                        reader.read_exact(&mut val_bytes)?;
                        if let Some(scalar) = bytes_to_scalar(&val_bytes) {
                            c.push((constraint_idx, col, scalar));
                        }
                    }
                }

                return Ok((a, b, c));
            } else {
                reader.seek(SeekFrom::Current(section_size as i64))?;
            }
        }

        Err(R1CSError::SectionNotFound(2))
    }

    /// Get statistics about the R1CS
    pub fn stats(&self) -> R1CSStats {
        R1CSStats {
            num_constraints: self.num_constraints,
            num_variables: self.num_variables,
            num_pub_inputs: self.num_pub_inputs,
            num_prv_inputs: self.num_prv_inputs,
            nnz_a: self.a.len(),
            nnz_b: self.b.len(),
            nnz_c: self.c.len(),
        }
    }
    
    /// Convert to sparse matrix format for Spartan with column remapping
    /// 
    /// IMPORTANT: The remapping uses PADDED num_vars for column indices:
    /// - Circom: col 0 = constant 1, cols 1..n_pub+1 = public inputs, cols n_pub+1.. = private vars
    /// - Spartan: cols 0..num_vars_padded = private vars, col num_vars_padded = constant 1, cols num_vars_padded+1.. = public inputs
    /// 
    /// Returns (A, B, C) as vectors of (row, col, value) triplets
    pub fn to_sparse_matrices_padded(&self, num_vars_padded: usize) -> (
        Vec<(usize, usize, Scalar)>,
        Vec<(usize, usize, Scalar)>,
        Vec<(usize, usize, Scalar)>,
    ) {
        let n_pub = self.num_pub_inputs;
        
        // Remap column indices from circom to Spartan format using PADDED dimensions
        // circom col 0 (constant 1) -> spartan col num_vars_padded
        // circom cols 1..n_pub+1 (public inputs) -> spartan cols num_vars_padded+1..
        // circom cols n_pub+1.. (private vars) -> spartan cols 0..
        let remap_col = |col: usize| -> usize {
            if col == 0 {
                // constant 1 -> position after all private vars (padded)
                num_vars_padded
            } else if col <= n_pub {
                // public inputs -> after constant
                num_vars_padded + col
            } else {
                // private vars -> at the beginning (col - 1 - n_pub)
                col - n_pub - 1
            }
        };
        
        let convert = |mat: &[(usize, usize, Scalar)]| -> Vec<(usize, usize, Scalar)> {
            mat.iter().map(|&(row, col, val)| (row, remap_col(col), val)).collect()
        };
        
        (convert(&self.a), convert(&self.b), convert(&self.c))
    }
    
    /// Convert to sparse matrix format (legacy - uses unpadded num_vars)
    /// Prefer to_sparse_matrices_padded for proper Spartan compatibility
    pub fn to_sparse_matrices(&self) -> (
        Vec<(usize, usize, Scalar)>,
        Vec<(usize, usize, Scalar)>,
        Vec<(usize, usize, Scalar)>,
    ) {
        self.to_sparse_matrices_padded(self.num_private_vars())
    }
    
    /// Get the number of private variables (for Spartan's num_vars)
    pub fn num_private_vars(&self) -> usize {
        self.num_variables - 1 - self.num_pub_inputs
    }
}

#[derive(Debug, Clone)]
pub struct R1CSStats {
    pub num_constraints: usize,
    pub num_variables: usize,
    pub num_pub_inputs: usize,
    pub num_prv_inputs: usize,
    pub nnz_a: usize,
    pub nnz_b: usize,
    pub nnz_c: usize,
}

impl std::fmt::Display for R1CSStats {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "R1CS Statistics:")?;
        writeln!(f, "  Constraints: {}", self.num_constraints)?;
        writeln!(f, "  Variables:   {}", self.num_variables)?;
        writeln!(f, "  Public inputs:  {}", self.num_pub_inputs)?;
        writeln!(f, "  Private inputs: {}", self.num_prv_inputs)?;
        writeln!(f, "  Non-zeros in A: {}", self.nnz_a)?;
        writeln!(f, "  Non-zeros in B: {}", self.nnz_b)?;
        writeln!(f, "  Non-zeros in C: {}", self.nnz_c)?;
        writeln!(f, "  Total non-zeros: {}", self.nnz_a + self.nnz_b + self.nnz_c)
    }
}

/// Convert little-endian bytes to BN254 scalar
fn bytes_to_scalar(bytes: &[u8]) -> Option<Scalar> {
    if bytes.len() != 32 {
        return None;
    }
    
    // Convert to limbs (little-endian)
    let mut limbs = [0u64; 4];
    for i in 0..4 {
        limbs[i] = u64::from_le_bytes(bytes[i * 8..(i + 1) * 8].try_into().ok()?);
    }
    
    let bigint = ark_ff::BigInt(limbs);
    Fr::from_bigint(bigint).map(Scalar)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bytes_to_scalar() {
        // Test with value 1
        let mut bytes = [0u8; 32];
        bytes[0] = 1;
        let s = bytes_to_scalar(&bytes).unwrap();
        assert_eq!(s, Scalar::one());

        // Test with value 0
        let bytes = [0u8; 32];
        let s = bytes_to_scalar(&bytes).unwrap();
        assert_eq!(s, Scalar::zero());
    }
}

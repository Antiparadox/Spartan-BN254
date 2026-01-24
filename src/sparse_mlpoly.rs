//! Sparse multilinear polynomial representation for R1CS matrices

use crate::hyrax::EqPolynomial;
use crate::math::Math;
use crate::scalar::Scalar;
use serde::{Deserialize, Serialize};

/// Sparse matrix entry (row, col, value)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SparseMatEntry {
    pub row: usize,
    pub col: usize,
    pub val: Scalar,
}

impl SparseMatEntry {
    pub fn new(row: usize, col: usize, val: Scalar) -> Self {
        SparseMatEntry { row, col, val }
    }
    
    pub fn row(&self) -> usize {
        self.row
    }
    
    pub fn col(&self) -> usize {
        self.col
    }
    
    pub fn val(&self) -> Scalar {
        self.val
    }
}

/// Sparse matrix polynomial representation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SparseMatPolynomial {
    pub num_vars_x: usize,
    pub num_vars_y: usize,
    pub M: Vec<SparseMatEntry>,
}

impl SparseMatPolynomial {
    pub fn new(num_vars_x: usize, num_vars_y: usize, M: Vec<SparseMatEntry>) -> Self {
        SparseMatPolynomial {
            num_vars_x,
            num_vars_y,
            M,
        }
    }

    pub fn num_entries(&self) -> usize {
        self.M.len()
    }
    
    pub fn entries(&self) -> &[SparseMatEntry] {
        &self.M
    }
    
    pub fn get_entries(&self) -> &[SparseMatEntry] {
        &self.M
    }
    
    pub fn num_vars_x(&self) -> usize {
        self.num_vars_x
    }
    
    pub fn num_vars_y(&self) -> usize {
        self.num_vars_y
    }

    pub fn get_num_nz_entries(&self) -> usize {
        self.M.len().next_power_of_two()
    }

    /// Multiply the sparse matrix by a dense vector z
    /// Returns M * z
    pub fn multiply_vec(&self, num_rows: usize, num_cols: usize, z: &[Scalar]) -> Vec<Scalar> {
        assert_eq!(z.len(), num_cols);
        
        let mut result = vec![Scalar::zero(); num_rows];
        for entry in &self.M {
            if entry.col < num_cols {
                result[entry.row] += entry.val * z[entry.col];
            }
        }
        result
    }

    /// Evaluate the multilinear extension at point (rx, ry)
    pub fn evaluate(&self, rx: &[Scalar], ry: &[Scalar]) -> Scalar {
        let eq_rx = EqPolynomial::new(rx.to_vec()).evals();
        let eq_ry = EqPolynomial::new(ry.to_vec()).evals();

        self.M
            .iter()
            .map(|entry| {
                let val_x = if entry.row < eq_rx.len() {
                    eq_rx[entry.row]
                } else {
                    Scalar::zero()
                };
                let val_y = if entry.col < eq_ry.len() {
                    eq_ry[entry.col]
                } else {
                    Scalar::zero()
                };
                entry.val * val_x * val_y
            })
            .sum()
    }

    /// Evaluate multiple sparse polynomials at the same point
    pub fn multi_evaluate(
        polys: &[&SparseMatPolynomial],
        rx: &[Scalar],
        ry: &[Scalar],
    ) -> Vec<Scalar> {
        let eq_rx = EqPolynomial::new(rx.to_vec()).evals();
        let eq_ry = EqPolynomial::new(ry.to_vec()).evals();

        polys
            .iter()
            .map(|poly| {
                poly.M
                    .iter()
                    .map(|entry| {
                        let val_x = if entry.row < eq_rx.len() {
                            eq_rx[entry.row]
                        } else {
                            Scalar::zero()
                        };
                        let val_y = if entry.col < eq_ry.len() {
                            eq_ry[entry.col]
                        } else {
                            Scalar::zero()
                        };
                        entry.val * val_x * val_y
                    })
                    .sum()
            })
            .collect()
    }

    /// Compute evaluation tables for the sparse polynomial
    pub fn compute_eval_table_sparse(
        &self,
        evals: &[Scalar],
        _num_rows: usize,
        num_cols: usize,
    ) -> Vec<Scalar> {
        let mut table = vec![Scalar::zero(); num_cols];
        
        for entry in &self.M {
            if entry.row < evals.len() && entry.col < num_cols {
                table[entry.col] += evals[entry.row] * entry.val;
            }
        }
        
        table
    }

    /// Evaluate using precomputed eq tables - flexible version
    pub fn evaluate_with_tables(&self, rx_evals: &[Scalar], ry_evals: &[Scalar]) -> Scalar {
        self.M
            .iter()
            .map(|entry| {
                let val_x = if entry.row < rx_evals.len() {
                    rx_evals[entry.row]
                } else {
                    Scalar::zero()
                };
                let val_y = if entry.col < ry_evals.len() {
                    ry_evals[entry.col]
                } else {
                    Scalar::zero()
                };
                entry.val * val_x * val_y
            })
            .sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sparse_mat_multiply() {
        // Create a 2x2 identity matrix
        let entries = vec![
            SparseMatEntry::new(0, 0, Scalar::one()),
            SparseMatEntry::new(1, 1, Scalar::one()),
        ];
        let poly = SparseMatPolynomial::new(1, 1, entries);
        
        let z = vec![Scalar::from_u64(3), Scalar::from_u64(5)];
        let result = poly.multiply_vec(2, 2, &z);
        
        assert_eq!(result[0], Scalar::from_u64(3));
        assert_eq!(result[1], Scalar::from_u64(5));
    }

    #[test]
    fn test_sparse_mat_evaluate() {
        // Single entry at (0, 0) with value 7
        let entries = vec![SparseMatEntry::new(0, 0, Scalar::from_u64(7))];
        let poly = SparseMatPolynomial::new(1, 1, entries);
        
        // Evaluate at (0, 0) -> should give 7
        let rx = vec![Scalar::zero()];
        let ry = vec![Scalar::zero()];
        let result = poly.evaluate(&rx, &ry);
        assert_eq!(result, Scalar::from_u64(7));
    }
}

//! Math utilities for Spartan

/// Extension trait for math operations on usize
pub trait Math {
    /// Returns the number of bits needed to represent this number
    fn log_2(self) -> usize;
    /// Returns 2^self
    fn pow2(self) -> usize;
}

impl Math for usize {
    fn log_2(self) -> usize {
        assert!(self > 0);
        (std::mem::size_of::<usize>() * 8) - (self.leading_zeros() as usize) - 1
    }

    fn pow2(self) -> usize {
        1 << self
    }
}

/// Convert a number to its binary representation as a vector of bools
pub fn to_bits(val: usize, num_bits: usize) -> Vec<bool> {
    (0..num_bits)
        .map(|shift_amount| (val & (1 << (num_bits - shift_amount - 1))) > 0)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_log_2() {
        assert_eq!(1usize.log_2(), 0);
        assert_eq!(2usize.log_2(), 1);
        assert_eq!(4usize.log_2(), 2);
        assert_eq!(8usize.log_2(), 3);
        assert_eq!(1024usize.log_2(), 10);
    }

    #[test]
    fn test_pow2() {
        assert_eq!(0usize.pow2(), 1);
        assert_eq!(1usize.pow2(), 2);
        assert_eq!(10usize.pow2(), 1024);
    }

    #[test]
    fn test_to_bits() {
        assert_eq!(to_bits(0, 4), vec![false, false, false, false]);
        assert_eq!(to_bits(1, 4), vec![false, false, false, true]);
        assert_eq!(to_bits(5, 4), vec![false, true, false, true]);
        assert_eq!(to_bits(15, 4), vec![true, true, true, true]);
    }
}

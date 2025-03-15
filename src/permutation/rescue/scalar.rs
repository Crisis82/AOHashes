use dusk_bls12_381::BlsScalar;
use dusk_safe::Safe;

use super::RescueHash;
use crate::rescue::{ALPHA, ALPHA_INV, MDS_MATRIX, ROUND_CONSTANTS, WIDTH};

#[derive(Default)]
pub(crate) struct ScalarPermutation();

impl ScalarPermutation {
    /// Constructs a new `ScalarPermutation`.
    pub fn new() -> Self {
        Self()
    }
}

impl Safe<BlsScalar, WIDTH> for ScalarPermutation {
    fn permute(&mut self, state: &mut [BlsScalar; WIDTH]) {
        self.perm(state);
    }

    fn tag(&mut self, input: &[u8]) -> BlsScalar {
        BlsScalar::hash_to_scalar(input.as_ref())
    }

    fn add(&mut self, right: &BlsScalar, left: &BlsScalar) -> BlsScalar {
        right + left
    }
}

impl RescueHash<BlsScalar> for ScalarPermutation {
    fn sbox_layer(&mut self, state: &mut [BlsScalar]) {
        state
            .iter_mut()
            .for_each(|scalar| *scalar = scalar.pow_vartime(&ALPHA));
    }

    fn sbox_inv_layer(&mut self, state: &mut [BlsScalar]) {
        state
            .iter_mut()
            .for_each(|scalar| *scalar = scalar.pow_vartime(&*ALPHA_INV));
    }

    fn linear_layer(&mut self, state: &mut [BlsScalar]) {
        let mut result = [BlsScalar::zero(); WIDTH];

        for (j, scalar) in state.iter().enumerate() {
            for k in 0..WIDTH {
                result[k] += MDS_MATRIX[k][j] * scalar;
            }
        }

        state.copy_from_slice(&result);
    }

    fn constant_injection(
        &mut self,
        state: &mut [BlsScalar],
        round: usize,
        offset: usize,
    ) {
        state.iter_mut().enumerate().for_each(|(i, scalar)| {
            *scalar += ROUND_CONSTANTS[round][i + offset]
        });
    }
}

#[cfg(feature = "encryption")]
impl dusk_safe::Encryption<BlsScalar, WIDTH> for ScalarPermutation {
    fn subtract(
        &mut self,
        minuend: &BlsScalar,
        subtrahend: &BlsScalar,
    ) -> BlsScalar {
        minuend - subtrahend
    }

    fn is_equal(&mut self, lhs: &BlsScalar, rhs: &BlsScalar) -> bool {
        lhs == rhs
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rescue_det() {
        let mut x = [BlsScalar::from(17u64); WIDTH];
        let mut y = [BlsScalar::from(17u64); WIDTH];
        let mut z = [BlsScalar::from(19u64); WIDTH];

        ScalarPermutation::new().permute(&mut x);
        ScalarPermutation::new().permute(&mut y);
        ScalarPermutation::new().permute(&mut z);

        assert_eq!(x, y);
        assert_ne!(x, z);
    }
}

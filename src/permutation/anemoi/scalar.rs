use dusk_bls12_381::BlsScalar;
use dusk_safe::Safe;

use super::AnemoiHash;
use crate::anemoi::{
    ALPHA_INV, G, G_INV, L, MATRIX_X, MATRIX_Y, ROUND_CONSTANTS, WIDTH,
};

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

impl AnemoiHash<BlsScalar> for ScalarPermutation {
    fn constant_addition(&mut self, state: &mut [BlsScalar], round: usize) {
        for (i, scalar) in state.iter_mut().enumerate() {
            *scalar += ROUND_CONSTANTS[round][i];
        }
    }

    fn linear_layer(&mut self, state: &mut [BlsScalar]) {
        let mut result = [BlsScalar::zero(); WIDTH];

        for i in 0..L {
            for j in 0..L {
                result[i] += MATRIX_X[i][j] * state[j];
                result[i + L] += MATRIX_Y[i][j] * state[j + L];
            }
        }

        state.copy_from_slice(&result);
    }

    fn pht(&mut self, state: &mut [BlsScalar]) {
        // Applies Y <- Y + X
        for i in L..WIDTH {
            state[i] += state[i - L];
        }
        // Applies X <- X + Y, with Y updated to Y+X from previous step
        // thus is equivalent to X <- 2*X + Y
        for i in 0..L {
            state[i] += state[i + L];
        }
    }

    fn sbox_layer(&mut self, state: &mut [BlsScalar]) {
        for i in 0..L {
            // x_i = x_i - g * y_i^2 - g_inv
            state[i] = state[i] - G * (state[i + L].square()) - *G_INV;
            // y_i = y_i - x_i^alpha_inv
            state[i + L] -= state[i].pow_vartime(&ALPHA_INV);
            // x_i = x_i + g * y_i^2
            state[i] += G * (state[i + L].square());
        }
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
    fn anemoi_det() {
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

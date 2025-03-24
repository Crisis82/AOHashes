use dusk_bls12_381::BlsScalar;
use dusk_safe::Safe;

use super::GriffinPi;
use crate::griffin::{ALPHA, BETA, D, D_INV, MATRIX, ROUND_CONSTANTS, WIDTH};

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

impl GriffinPi<BlsScalar> for ScalarPermutation {
    fn non_linear_layer(&mut self, state: &mut [BlsScalar]) {
        let mut result = [BlsScalar::zero(); WIDTH];

        for i in 0..WIDTH {
            if i == 0 {
                result[i] = state[i].pow_vartime(&*D_INV);
            } else if i == 1 {
                result[i] = state[i].pow_vartime(&D);
            } else {
                let l = if i == 2 {
                    // gamma_i = i -1 = 2 -1 = 1
                    // l_2 = gamma_i * z_0 + z_1 + 0 = z_0 + z_1
                    result[0] + result[1]
                } else {
                    // gamma_i = i -1
                    // l_i = gamma_i * z_0 + z_1 + state_(i-1)
                    BlsScalar::from((i - 1) as u64) * result[0]
                        + result[1]
                        + state[i - 1]
                };
                // z_i = state_i * (l^2 + alpha_i * l + beta_i)
                result[i] = state[i] * (l.square() + ALPHA[i] * l + BETA[i]);
            }
        }

        state.copy_from_slice(&result);
    }

    fn linear_layer(&mut self, state: &mut [BlsScalar], _round: usize) {
        let mut result = [BlsScalar::zero(); WIDTH];

        for i in 0..WIDTH {
            for j in 0..WIDTH {
                result[i] += MATRIX[i][j] * state[j];
            }
        }

        state.copy_from_slice(&result);
    }

    fn constant_addition(&mut self, state: &mut [BlsScalar], round: usize) {
        state.iter_mut().for_each(|scalar| {
            *scalar += ROUND_CONSTANTS[round];
        })
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
    fn griffin_det() {
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

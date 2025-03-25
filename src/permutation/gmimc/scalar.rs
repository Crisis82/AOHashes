use dusk_bls12_381::BlsScalar;
use dusk_safe::Safe;

use super::Gmimc;
use crate::gmimc::{ROUND_CONSTANTS, WIDTH};

/// An implementation of the GMiMC permutation for `BlsScalar` as input values.
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

impl Gmimc<BlsScalar> for ScalarPermutation {
    fn f(&mut self, round: usize, value: &mut BlsScalar) {
        // x + c_i
        *value = *value + ROUND_CONSTANTS[round];
        // ^5
        *value = value.square().square() * *value;
    }

    fn erf(&mut self, state: &mut [BlsScalar]) {
        let result = state[0];
        state
            .iter_mut()
            .skip(1)
            .for_each(|scalar| *scalar += result);
    }

    fn left_rot(&mut self, state: &mut [BlsScalar]) {
        state.rotate_left(1);
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
    fn gmimc_det() {
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

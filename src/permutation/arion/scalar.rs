use dusk_bls12_381::BlsScalar;
use dusk_safe::Safe;

use super::ArionHash;
use crate::arion::{
    AFFINE_CONSTANTS, ALPHA1, ALPHA2, BETA, D1, D2_INV, MATRIX, WIDTH,
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

impl ArionHash<BlsScalar> for ScalarPermutation {
    fn gtds(&mut self, round: usize, state: &mut [BlsScalar]) {
        let mut result = [BlsScalar::zero(); WIDTH];
        result.copy_from_slice(&state);

        // High degree inverse S-Box is inverse of x^D2
        result[WIDTH - 1].pow_vartime(&*D2_INV);

        let mut s = state[WIDTH - 1].clone();
        s += result[WIDTH - 1];

        for i in (0..WIDTH - 1).rev() {
            result[i] = result[i].pow(&D1);

            // Evaluate g and h
            let tmp = s.square();
            // add linear term
            let mut g = s.clone();
            g *= ALPHA1[round][i];
            let mut h = s.clone();
            h *= BETA[round][i];
            // add quadratic term
            g += tmp;
            h += tmp;
            // add constant_term
            g += ALPHA2[round][i];

            // Multply g and add h
            result[i] *= g;
            result[i] += h;

            s += result[i];
            s += state[i];
        }

        state.copy_from_slice(&result);
    }

    fn mul_matrix(&mut self, state: &mut [BlsScalar]) {
        let mut result = [BlsScalar::zero(); WIDTH];

        for i in 0..WIDTH {
            for j in 0..WIDTH {
                result[i] += MATRIX[i][j] * state[j];
            }
        }

        state.copy_from_slice(&result);
    }

    fn affine_layer(&mut self, round: usize, state: &mut [BlsScalar]) {
        self.mul_matrix(state);

        let mut result = [BlsScalar::zero(); WIDTH];
        for i in 0..WIDTH {
            result[i] += state[i] + AFFINE_CONSTANTS[round][i];
        }

        state.copy_from_slice(&result);
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
    use crate::arion::D2;
    use ff::Field;
    use rand::{SeedableRng, rngs::StdRng};

    // High degree inverse S-Box is inverse of `x^D2`
    fn high_degree_permutation(samples: usize) {
        let mut rng = StdRng::from_entropy();
        for _i in 1..samples {
            let val = BlsScalar::random(&mut rng);
            let val_2 = val.pow_vartime(&*D2_INV);

            assert_eq!(val_2.pow_vartime(&D2), val);
        }
    }

    #[test]
    fn high_degree_permutation_1000() {
        high_degree_permutation(1000);
    }

    #[test]
    fn arion_det() {
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

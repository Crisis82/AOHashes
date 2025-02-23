use crate::griffin::{ROUNDS, WIDTH};
use crate::permutation::{BlsScalar, Shake256PRNG};

use alloc::vec;
use ndarray::{Array2, s};

/// Round, alpha and beta constants Griffin's generation function based on
/// SHAKE256 PRNG with fixed seed.
pub(crate) fn gen_const() -> [BlsScalar; ROUNDS + 2 * WIDTH] {
    let mut constants = [BlsScalar::zero(); ROUNDS];

    // SHAKE256 PRNG with fixed seed
    let seed = b"griffin_seed";
    let mut prng = Shake256PRNG::new(seed);

    // Generating the ROUND CONSTANTS
    // last round constant is zero
    for i in 0..(ROUNDS - 1) {
        constants[i] = prng.next_scalar();
    }

    // Generating ALPHA and BETA CONSTANTS
    // The first two elements remains zero, but are included just for
    // indexing convenience. We could remove them and adjust indexing
    // in the non-linear layer of Griffin-pi permutation
    let mut alpha = [BlsScalar::zero(); WIDTH];
    let mut beta = [BlsScalar::zero(); WIDTH];

    alpha[2] = prng.next_scalar();
    beta[2] = prng.next_scalar();

    for i in 3..WIDTH {
        let coeff = BlsScalar::from((i - 1) as u64);
        alpha[i] = coeff * alpha[2];
        beta[i] = coeff.square() * beta[2];
    }

    // writing results in a single matrix, to return all values
    let mut result = [BlsScalar::zero(); ROUNDS + 2 * WIDTH];
    for i in 0..ROUNDS {
        result[i] = constants[i];
    }
    for i in 0..WIDTH {
        result[ROUNDS + i] = alpha[i];
        result[ROUNDS + WIDTH + i] = beta[i];
    }

    result
}

/// Griffin matrix generation function.
pub(crate) fn gen_matrix() -> [[BlsScalar; WIDTH]; WIDTH] {
    let mut m: Array2<u64>;

    if WIDTH == 3 {
        #[allow(unused_unsafe)]
        unsafe {
            m = Array2::from_shape_vec((3, 3), vec![2, 1, 1, 1, 2, 1, 1, 1, 2])
                .unwrap();
        }
    } else {
        #[allow(unused_unsafe)]
        unsafe {
            m = Array2::from_shape_vec((4, 4), vec![
                5, 7, 1, 3, 4, 6, 1, 1, 1, 3, 5, 7, 1, 1, 4, 6,
            ])
            .unwrap();
        }
        if WIDTH > 4 {
            // Construct M' diagonal matrix
            let mut m1: Array2<u64> = Array2::zeros((WIDTH, WIDTH));
            let mut i = 0;
            while i < (WIDTH - 3) {
                m1.slice_mut(s![i..(i + 4), i..(i + 4)]).assign(&m);
                i += 4;
            }

            // Construct M'' circular matrix
            let identity: Array2<u64> = Array2::eye(4);
            let mut m2: Array2<u64> = Array2::zeros((WIDTH, WIDTH));

            let mut i = 0;
            while i < (WIDTH - 3) {
                m2.slice_mut(s![i..(i + 4), i..(i + 4)])
                    .assign(&(2 * &identity));

                let mut j = 0;
                while j < (WIDTH - 3) {
                    if j == i {
                        j += 4;
                        continue;
                    }
                    m2.slice_mut(s![i..(i + 4), j..(j + 4)]).assign(&identity);
                    j += 4;
                }
                i += 4;
            }

            // Compute M = M' * M''
            m = m1.dot(&m2);
        }
    }

    let mut matrix = [[BlsScalar::zero(); WIDTH]; WIDTH];
    for i in 0..WIDTH {
        for j in 0..WIDTH {
            matrix[i][j] = BlsScalar::from(m[[i, j]]);
        }
    }

    matrix
}

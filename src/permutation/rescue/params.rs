use crate::permutation::{BlsScalar, Shake256PRNG};
use crate::rescue::{ROUNDS, WIDTH};

/// Rescue round constants generation function using SHAKE256 with a fixed
/// seed.
pub(crate) fn gen_const() -> [[BlsScalar; 2 * WIDTH]; ROUNDS] {
    let mut constants = [[BlsScalar::zero(); 2 * WIDTH]; ROUNDS];

    // fixed seed
    let seed = b"rescue_fixed_seed";

    // PRNG based on SHAKE256
    let mut prng = Shake256PRNG::new(seed);

    // generation of the scalar constants
    for row in constants.iter_mut() {
        for scalar in row.iter_mut() {
            *scalar = prng.next_scalar();
        }
    }

    constants
}

/// Reduced Row Echelon Form of the given matrix.
fn rref(matrix: &mut [[BlsScalar; 2 * WIDTH]; WIDTH]) {
    let mut lead = 0;

    for r in 0..WIDTH {
        if lead >= (2 * WIDTH) {
            return;
        }

        let mut i = r;
        while matrix[i][lead] == BlsScalar::zero() {
            i += 1;
            if i == WIDTH {
                i = r;
                lead += 1;
                if lead == (2 * WIDTH) {
                    return;
                }
            }
        }

        matrix.swap(i, r);

        let lead_val = matrix[r][lead];
        if lead_val != BlsScalar::zero() {
            for j in 0..(2 * WIDTH) {
                // divide by using multiplicative inverse
                let inv_lead_val = lead_val.invert().unwrap();
                matrix[r][j] = matrix[r][j] * inv_lead_val;
            }
        }

        for i in 0..WIDTH {
            if i != r {
                let factor = matrix[i][lead];
                for j in 0..(2 * WIDTH) {
                    matrix[i][j] -= factor * matrix[r][j];
                }
            }
        }
        lead += 1;
    }
}

/// Extract the right half of the matrix.
fn right_half(
    matrix: [[BlsScalar; 2 * WIDTH]; WIDTH],
) -> [[BlsScalar; WIDTH]; WIDTH] {
    let mut right_half = [[BlsScalar::zero(); WIDTH]; WIDTH];
    for i in 0..WIDTH {
        for j in 0..WIDTH {
            right_half[i][j] = matrix[i][j + WIDTH];
        }
    }
    right_half
}

/// Matrix transposition.
fn transpose(
    matrix: [[BlsScalar; WIDTH]; WIDTH],
) -> [[BlsScalar; WIDTH]; WIDTH] {
    let mut transposed = [[BlsScalar::zero(); WIDTH]; WIDTH];
    for i in 0..WIDTH {
        for j in 0..WIDTH {
            transposed[j][i] = matrix[i][j];
        }
    }
    transposed
}

/// Rescue MDS matrix generation function.
pub(crate) fn gen_mds() -> [[BlsScalar; WIDTH]; WIDTH] {
    let mut matrix = [[BlsScalar::zero(); 2 * WIDTH]; WIDTH];

    // Primitive generator for p
    let g = BlsScalar::from(7);

    // Vandermonde matrix
    for i in 0..WIDTH {
        for j in 0..(2 * WIDTH) {
            let exp: [u64; 4] = [(i * j) as u64, 0, 0, 0];
            matrix[i][j] = g.pow(&exp);
        }
    }

    // transform in echelon form
    rref(&mut matrix);

    // return the transpose of the right half
    transpose(right_half(matrix))
}

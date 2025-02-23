use crate::arion::{ROUNDS, WIDTH};
use crate::permutation::{BlsScalar, u256_to_u64};
use ethnum::u256;
use ff::Field;
use rand::SeedableRng;
use rand::rngs::StdRng;

// P - 1
const P1: BlsScalar = BlsScalar::from_raw([
    0xffff_ffff_0000_0000,
    0x53bd_a402_fffe_5bfe,
    0x3339_d808_09a1_d805,
    0x73ed_a753_299d_7d48,
]);

fn legendre_symbol(a1: BlsScalar, a2: BlsScalar, exp: [u64; 4]) -> i8 {
    let val = a1.pow(&[2, 0, 0, 0]) - BlsScalar::from(4) * a2;

    let result = val.pow(&exp);

    if result == P1 {
        return -1;
    } else {
        return 1;
    }
}

/// ALPHA constants random generation function.
pub(crate) fn gen_alpha() -> [[[BlsScalar; WIDTH - 1]; ROUNDS]; 2] {
    let mut alpha1 = [[BlsScalar::zero(); WIDTH - 1]; ROUNDS];
    let mut alpha2 = [[BlsScalar::zero(); WIDTH - 1]; ROUNDS];

    // (P-1) / 2
    let e = u256::from_str_hex(
        "0x39f6d3a994cebea4199cec0404d0ec02a9ded2017fff2dff7fffffff80000000",
    )
    .expect("Invalid hex string.");
    let exp = u256_to_u64(e);

    let mut rng = StdRng::from_entropy();

    for i in 0..ROUNDS {
        for j in 0..(WIDTH - 1) {
            let mut a1 = BlsScalar::random(&mut rng);
            let mut a2 = BlsScalar::random(&mut rng);

            while legendre_symbol(a1, a2, exp) != -1 {
                a1 = BlsScalar::random(&mut rng);
                a2 = BlsScalar::random(&mut rng);
            }

            alpha1[i][j] = a1;
            alpha2[i][j] = a2;
        }
    }

    [alpha1, alpha2]
}

/// BETA constants random generation function.
pub(crate) fn gen_beta() -> [[BlsScalar; WIDTH - 1]; ROUNDS] {
    let mut beta = [[BlsScalar::zero(); WIDTH - 1]; ROUNDS];
    let mut rng = StdRng::from_entropy();

    for i in 0..ROUNDS {
        for j in 0..(WIDTH - 1) {
            beta[i][j] = BlsScalar::random(&mut rng);
        }
    }

    beta
}

/// Circulant matrix generation function.
pub(crate) fn gen_matrix() -> [[BlsScalar; WIDTH]; WIDTH] {
    let mut matrix = [[BlsScalar::zero(); WIDTH]; WIDTH];

    let one = BlsScalar::one();

    let mut i = 1;
    while i < WIDTH {
        matrix[0][i] = matrix[0][i - 1].add(&one);
        i += 1;
    }

    i = 1;
    let mut j = 1;
    while i < WIDTH {
        matrix[i][0] = matrix[i - 1][WIDTH - 1];
        while j < WIDTH {
            matrix[i][j] = matrix[i - 1][j - 1];
            j += 1;
        }
        j = 1;
        i += 1;
    }

    matrix
}

/// Affine constants random generation function.
pub(crate) fn gen_affine() -> [[BlsScalar; WIDTH]; ROUNDS + 1] {
    let mut affine_constants = [[BlsScalar::zero(); WIDTH]; ROUNDS + 1];
    let mut rng = StdRng::from_entropy();

    for i in 0..ROUNDS {
        for j in 0..WIDTH {
            affine_constants[i][j] = BlsScalar::random(&mut rng);
        }
    }

    affine_constants
}

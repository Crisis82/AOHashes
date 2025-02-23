use crate::anemoi::{ALPHA, G, G_INV, L, ROUNDS};
use crate::permutation::{BlsScalar, cauchy_matrix};

use alloc::vec::Vec;

const ONE: BlsScalar = BlsScalar::from_raw([1, 0, 0, 0]);

/// C and D round constants Anemoi's generation function.
pub(crate) fn gen_const() -> [[BlsScalar; 2 * L]; ROUNDS] {
    let pi_0: Vec<char> = "1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679".chars().collect();
    let pi_1: Vec<char> = "8214808651328230664709384460955058223172535940812848111745028410270193852110555964462294895493038196".chars().collect();

    let mut constants = [[BlsScalar::zero(); 2 * L]; ROUNDS];
    for i in 0..ROUNDS {
        let p0i = BlsScalar::from(pi_0[i] as u64);
        for j in 0..L {
            let p1j = BlsScalar::from(pi_1[j] as u64);
            let sum = (p0i + p1j).pow(&ALPHA);
            // c constants
            constants[i][j] = G * p0i.square() + sum;
            // d constants
            constants[i][j + L] = G * p1j.square() + sum + *G_INV;
        }
    }
    constants
}

/// Anemoi matrix generation function.
pub(crate) fn gen_matrix() -> [[BlsScalar; L]; L] {
    let mut matrix = [[BlsScalar::zero(); L]; L];

    if L == 1 {
        matrix[0][0] = ONE;
    } else if L == 2 {
        matrix[0][0] = ONE;
        matrix[0][1] = G;
        matrix[1][0] = G;
        matrix[1][1] = G.square() + ONE;
    } else if L == 3 {
        matrix[0][0] = G + ONE;
        matrix[0][1] = ONE;
        matrix[0][2] = G + ONE;
        matrix[1][0] = ONE;
        matrix[1][1] = ONE;
        matrix[1][2] = G;
        matrix[2][0] = G;
        matrix[2][1] = ONE;
        matrix[2][2] = ONE;
    } else if L == 4 {
        matrix[0][0] = ONE;
        matrix[0][1] = G.square();
        matrix[0][2] = G.square();
        matrix[0][3] = G + ONE;
        matrix[1][0] = G + ONE;
        matrix[1][1] = G.square() + G;
        matrix[1][2] = G.square();
        matrix[1][3] = G.double() + ONE;
        matrix[2][0] = G;
        matrix[2][1] = G + ONE;
        matrix[2][2] = ONE;
        matrix[2][3] = G;
        matrix[3][0] = G;
        matrix[3][1] = G.double() + ONE;
        matrix[3][2] = G + ONE;
        matrix[3][3] = G + ONE;
    } else {
        cauchy_matrix(&mut matrix);
    }

    matrix
}

use crate::permutation::inverse_mod;
use dusk_bls12_381::BlsScalar;
use lazy_static::lazy_static;

/// State width.
/// Can have value 3, 4, 5, 6, 8
pub const WIDTH: usize = 3;

/// D1 s-box exponent: smallest positive integer coprime with Q−1.
/// Usually have value 3 or 5.
pub(crate) const D1: [u64; 4] = [5, 0, 0, 0];

/// Number of rounds.
/// Depends on the exponent D1 and the state WIDTH.
pub(crate) const ROUNDS: usize = {
    match D1[0] {
        3 => match WIDTH {
            3 | 4 => 6,
            5 | 6 => 5,
            8 => 4,
            _ => panic!("Invalid WIDTH"),
        },
        5 => match WIDTH {
            3 => 6,
            4 | 5 | 6 => 5,
            8 => 4,
            _ => panic!("Invalid WIDTH"),
        },
        _ => panic!("Invalid value for D1"),
    }
};

/// D2 s-box exponent: arbitrary integer coprime with Q−1.
// Because of our choice of modulus, the modular invertible value
// between 121 and 257 are: 125, 161, 193, 257.
// Default value is 257.
pub(crate) const D2: [u64; 4] = [257, 0, 0, 0];

/// ArionHash parameters generation.
mod params;
use params::{gen_affine, gen_alpha, gen_beta, gen_matrix};

lazy_static! {
    static ref ALPHAS: [[[BlsScalar; WIDTH-1]; ROUNDS]; 2] = gen_alpha();
    /// Alpha_1 constant of function G.
    pub(crate) static ref ALPHA1: [[BlsScalar; WIDTH-1]; ROUNDS] = ALPHAS[0];
    /// Alpha_2 constant of function G.
    pub(crate) static ref ALPHA2: [[BlsScalar; WIDTH-1]; ROUNDS] = ALPHAS[1];

    /// Beta constant of function H.
    pub(crate) static ref BETA: [[BlsScalar; WIDTH-1]; ROUNDS] = gen_beta();

    /// Circulant matrix used by the affine layer.
    pub(crate) static ref MATRIX: [[BlsScalar; WIDTH]; WIDTH] = gen_matrix();

    /// Affine layer constants.
    /// It has a size of ROUNDS+1 to include in its last position (ROUNDS) the
    /// affine constant set to 0, for the initial call to the affine function.
    pub(crate) static ref AFFINE_CONSTANTS: [[BlsScalar; WIDTH]; ROUNDS+1] = gen_affine();

    /// Multiplicative inverse of D2 modulo Q-1.
    pub(crate) static ref D2_INV: [u64; 4] = inverse_mod(D2[0]);
}

/// ArionHash permutation struct operating in a plonk-circuit.
#[cfg(feature = "zk")]
pub(crate) mod gadget;

/// ArionHash permutation struct operating on [`BlsScalar`].
pub(crate) mod scalar;

/// Defines the ArionHash permutation algorithm.
pub(crate) trait ArionHash<T> {
    /// Generalized Triangular Dynamical System.
    fn gtds(&mut self, round: usize, state: &mut [T]);

    /// Matrix-vector multiplication using the precomputed circulant matrix.
    fn mul_matrix(&mut self, state: &mut [T]);

    /// Affine layer: matrix multiplication followed by a constant addition.
    fn affine_layer(&mut self, round: usize, state: &mut [T]);

    /// Applies a round.
    ///
    /// One round consists of:
    /// - GTDS (Generalized Triangular Dynamical System):
    ///     - `GTDS(x_i) = x_i^D1 + G(delta_{i+1}) + H(delta_{i+1})`,
    ///     - `GTDS(x_n) = x_n^D2_INV`;
    /// - Affine layer L: `L(x) = CircM * x + c`.
    fn apply_round(&mut self, round: usize, state: &mut [T]) {
        self.gtds(round, state);
        self.affine_layer(round, state);
    }

    /// Applies one permutation.
    fn perm(&mut self, state: &mut [T]) {
        // We initialize the permutation by applying a single
        // affine layer with affine constant set to 0.
        self.affine_layer(ROUNDS, state);
        // Apply R rounds
        for round in 0..ROUNDS {
            self.apply_round(round, state);
        }
    }
}

use crate::permutation::inverse_mod;
use dusk_bls12_381::BlsScalar;
use lazy_static::lazy_static;

/// Single state width.
/// Can have value 1, 2, 3, 4, 6, 8.
pub(crate) const L: usize = 4;

/// Total state width, which corresponds to the
/// concatenation of two single states X and Y of width L.
pub const WIDTH: usize = 2 * L;

/// Alpha s-box exponent: smallest integer comprime with Q-1.
/// Can have value 5 or 7.
pub(crate) const ALPHA: [u64; 4] = [5, 0, 0, 0];

/// Number of rounds.
/// Depends of the value of ALPHA and L.
pub(crate) const ROUNDS: usize = {
    match ALPHA[0] {
        5 => match L {
            1 => 21,
            2 => 14,
            3 => 12,
            4 => 12,
            6 | 8 => 10,
            _ => panic!("Invalid value for L"),
        },
        7 => match L {
            1 => 20,
            2 => 13,
            3 => 12,
            4 => 11,
            6 => 10,
            8 => 9,
            _ => panic!("Invalid value for L"),
        },
        _ => panic!("Invalid value for ALPHA"),
    }
};

/// Generator G for Q.
pub(crate) const G: BlsScalar = BlsScalar::from_raw([7, 0, 0, 0]);

/// AnemoiHash parameters generation.
mod params;
use params::{gen_const, gen_matrix};

lazy_static! {
    /// Multiplicative inverse of G.
    pub(crate) static ref G_INV: BlsScalar = G.invert().unwrap();

    /// Modular inverse of alpha modulo Q-1.
    pub(crate) static ref ALPHA_INV: [u64; 4] = inverse_mod(ALPHA[0]);

    /// Round constants for both state X and Y.
    pub(crate) static ref ROUND_CONSTANTS: [[BlsScalar; WIDTH]; ROUNDS] = gen_const();

    /// Matrix for state X.
    pub(crate) static ref MATRIX_X: [[BlsScalar; L]; L] = gen_matrix();

    /// Matrix for state Y.
    /// Instead of permuting state Y and multiplying it with MATRIX_X, we use a new
    /// matrix which is the permutation of MATRIX_X: [row1, row2, row3, row0].
    pub(crate) static ref MATRIX_Y: [[BlsScalar; L]; L] = {
        let mut matrix_y = [[BlsScalar::zero(); L]; L];

        let mut i = 0;
        while i < L {
            let mut j = 0;
            while j < (L - 1) {
                matrix_y[i][j] = MATRIX_X[i][j + 1];
                j += 1;
            }
            matrix_y[i][j] = MATRIX_X[i][0];
            i += 1;
        }

        matrix_y
    };
}

/// AnemoiHash permutation struct operating in a plonk-circuit.
#[cfg(feature = "zk")]
pub(crate) mod gadget;

/// AnemoiHash permutation struct operating on [`BlsScalar`].
pub(crate) mod scalar;

/// Defines the AnemoiHash permutation algorithm.
pub(crate) trait AnemoiHash<T> {
    /// Adds the round constants to the state.
    fn constant_addition(&mut self, state: &mut [T], round: usize);

    /// Matrix-vector multiplication, first between matrix X and state X, then
    /// between matrix Y and state Y.
    fn linear_layer(&mut self, state: &mut [T]);

    /// Pseudo-Hadamard Transformation P.
    fn pht(&mut self, state: &mut [T]);

    /// Open Flystel function H.
    fn sbox_layer(&mut self, state: &mut [T]);

    /// Applies a round.
    ///
    /// One round consists of:
    /// - Constant addition: addition of constants c to state X, and of
    ///   constants d to state Y;
    /// - Linear layer: matrix-vector multiplication between MATRIX_X and state
    ///   X, and between MATRIX_Y and state Y;
    /// - PHT: Pseudo-Hadamard Transformation P;
    /// - S-Box layer: application of the open Flystel function H with
    ///   parameters alpha, beta, gamma and delta.
    fn apply_round(&mut self, round: usize, state: &mut [T]) {
        self.constant_addition(state, round);
        self.linear_layer(state);
        self.pht(state);
        self.sbox_layer(state);
    }

    /// Applies one permutation.
    fn perm(&mut self, state: &mut [T]) {
        // Apply R rounds
        for round in 0..ROUNDS {
            self.apply_round(round, state);
        }
        // Conclude the permutation with a final linear layer
        self.linear_layer(state);
    }
}

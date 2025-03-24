use crate::permutation::inverse_mod;
use dusk_bls12_381::BlsScalar;
use lazy_static::lazy_static;

/// State width t, which can be
/// 3, 4, 8, 12, 16, 20, 24.
pub const WIDTH: usize = 4;

/// Exponent D, which can be 5 or 7.
pub(crate) const D: [u64; 4] = [5, 0, 0, 0];

/// Number of rounds, which is dependent on D and WIDTH.
// pub(crate) const ROUNDS: usize = 10;
pub(crate) const ROUNDS: usize = {
    match D[0] {
        5 => match WIDTH {
            3 => 14,
            4 => 11,
            8 | 12 | 16 | 20 | 24 => 9,
            _ => panic!("Invalid value for WIDTH"),
        },
        7 => match WIDTH {
            3 => 11,
            4 => 10,
            8 | 12 | 16 | 20 | 24 => 8,
            _ => panic!("Invalid value for WIDTH"),
        },
        _ => panic!("Invalid value for D"),
    }
};

/// Griffin-pi parameters generation.
mod params;
use params::{gen_const, gen_matrix};

lazy_static! {
    /// Modular inverse of exponent D modulo Q-1.
    pub(crate) static ref D_INV: [u64; 4] = inverse_mod(D[0]);

    pub(crate) static ref CONSTANTS: [BlsScalar; ROUNDS + 2 * WIDTH] = gen_const();
    /// Round constants.
    pub(crate) static ref ROUND_CONSTANTS: [BlsScalar; ROUNDS] = {
        let mut t = [BlsScalar::zero(); ROUNDS];
        t.copy_from_slice(&CONSTANTS[..ROUNDS]);
        t
    };
    /// Alpha constants.
    pub(crate) static ref ALPHA: [BlsScalar; WIDTH] = {
        let mut t = [BlsScalar::zero(); WIDTH];
        t.copy_from_slice(&CONSTANTS[ROUNDS..(ROUNDS + WIDTH)]);
        t
    };
    /// Beta constants.
    pub(crate) static ref BETA: [BlsScalar; WIDTH] = {
        let mut t = [BlsScalar::zero(); WIDTH];
        t.copy_from_slice(&CONSTANTS[(ROUNDS + WIDTH)..]);
        t
    };

    /// Invertible matrix.
    pub(crate) static ref MATRIX: [[BlsScalar; WIDTH]; WIDTH] = gen_matrix();
}

/// Griffin-pi permutation struct operating in a plonk-circuit.
#[cfg(feature = "zk")]
pub(crate) mod gadget;

/// Griffin-pi permutation struct operating on [`BlsScalar`].
pub(crate) mod scalar;

/// Defines the Griffin-pi permutation algorithm.
pub(crate) trait GriffinPi<T> {
    /// Substitution layer S.
    fn non_linear_layer(&mut self, state: &mut [T]);

    /// Matrix-vector multiplication with matrix M.
    fn linear_layer(&mut self, state: &mut [T], round: usize);

    /// Adds the round constants to the state.
    fn constant_addition(&mut self, state: &mut [T], round: usize);

    /// Applies a round.
    ///
    /// One round consists of:
    /// - **Non-linear layer S**: applies a different substitution depending on
    ///   the index;
    /// - **Linear layer**: matrix-vector multiplication between invertible
    ///   matrix M and state;
    /// - **Round constants addition**: add a constant depending on the round.
    fn apply_round(&mut self, round: usize, state: &mut [T]) {
        self.non_linear_layer(state);
        self.linear_layer(state, round);
        self.constant_addition(state, round);
    }

    /// Applies one permutation.
    fn perm(&mut self, state: &mut [T]) {
        // Begin the permutation with a linear layer
        self.linear_layer(state, ROUNDS - 1);
        // Apply R-1 rounds
        for round in 0..(ROUNDS - 1) {
            self.apply_round(round, state);
        }
        // Apply the last round without constant addition
        // Note: it's a small optimization because the last round constant is 0
        self.non_linear_layer(state);
        self.linear_layer(state, ROUNDS - 1);
    }
}

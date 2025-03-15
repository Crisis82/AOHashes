use crate::permutation::inverse_mod;
use dusk_bls12_381::BlsScalar;
use lazy_static::lazy_static;

/// State width, with possible values:
/// 3, 4, 5, 6, 8, 12, 16, 20 or 24.
pub const WIDTH: usize = 4;

/// Number of rounds, depending on width.
pub(crate) const ROUNDS: usize = {
    match WIDTH {
        3 => 14,
        4 => 11,
        5 => 9,
        6 | 8 | 12 | 16 | 20 | 24 => 8,
        _ => panic!("Invalid width"),
    }
};

/// Alpha s-box exponent.
pub(crate) const ALPHA: [u64; 4] = [5, 0, 0, 0];

/// Rescue parameters generation.
mod params;
use params::{gen_const, gen_mds};

lazy_static! {
    /// Round constants.
    pub(crate) static ref ROUND_CONSTANTS: [[BlsScalar; 2 * WIDTH]; ROUNDS] = gen_const();

    /// MDS square matrix.
    pub(crate) static ref MDS_MATRIX: [[BlsScalar; WIDTH]; WIDTH] = gen_mds();

    /// Multiplicative inverse of Alpha modulo Q-1.
    pub(crate) static ref ALPHA_INV: [u64; 4] = inverse_mod(ALPHA[0]);
}

/// Rescue permutation struct operating in a plonk-circuit.
#[cfg(feature = "zk")]
pub(crate) mod gadget;

/// Rescue permutation struct operating on [`BlsScalar`].
pub(crate) mod scalar;

/// Defines the Rescue permutation algorithm.
pub(crate) trait RescueHash<T> {
    /// Exponentiation to alpha of the state.
    fn sbox_layer(&mut self, state: &mut [T]);

    /// Exponentiation to alpha_inv of the state.
    fn sbox_inv_layer(&mut self, state: &mut [T]);

    /// Matrix-vector multiplication with the MDS matrix.
    fn linear_layer(&mut self, state: &mut [T]);

    /// Adds the round constants to the state.
    fn constant_injection(
        &mut self,
        state: &mut [T],
        round: usize,
        offset: usize,
    );

    /// Applies a round.
    ///
    /// One round consists of:
    /// - _First half_:
    ///     - **Inverse S-Box layer**: power map with alpha^(-1);
    ///     - **Linear layer**: matrix-vector multiplication with MDS matrix;
    ///     - **Constant injection**: round constant addition of WIDTH elements.
    /// - _Second half_:
    ///     - **S-Box layer**: power map with alpha;
    ///     - **Linear layer**: matrix-vector multiplication with MDS matrix;
    ///     - **Constant injection**: round constant addition of the next WIDTH
    ///       elements (same round).
    fn apply_round(&mut self, round: usize, state: &mut [T]) {
        // first half
        self.sbox_inv_layer(state);
        self.linear_layer(state);
        self.constant_injection(state, round, 0);

        // second half
        self.sbox_layer(state);
        self.linear_layer(state);
        self.constant_injection(state, round, WIDTH);
    }

    /// Applies one permutation.
    fn perm(&mut self, state: &mut [T]) {
        // Apply R rounds
        for round in 0..ROUNDS {
            self.apply_round(round, state);
        }
    }
}

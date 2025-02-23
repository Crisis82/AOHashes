use dusk_bls12_381::BlsScalar;
use lazy_static::lazy_static;

/// State width t, where t >=3 and
/// usually is equal to 3, 4, 5, 6 or 8.
pub const WIDTH: usize = 4;

/// Number of rounds, depending on the width.
pub(crate) const ROUNDS: usize = {
    match WIDTH {
        3 => 328,
        4 => 330,
        5 => 332,
        6 => 334,
        8 => 338,
        _ => panic!("Invalid WIDTH"),
    }
};

/// GMiMC parameters generation.
mod params;
use params::gen_const;

lazy_static! {
    /// Round constants.
    pub(crate) static ref ROUND_CONSTANTS: [BlsScalar; ROUNDS] = gen_const();
}

/// GMiMC permutation struct operating in a plonk-circuit.
#[cfg(feature = "zk")]
pub(crate) mod gadget;

/// GMiMC permutation struct operating on [`BlsScalar`].
pub(crate) mod scalar;

/// Defines the GMiMC permutation algorithm.
pub(crate) trait Gmimc<T> {
    /// Computes `f(x) => (x + c)^3`.
    fn f(&mut self, round: usize, value: &mut T);

    /// Expanding Round Function.
    fn erf(&mut self, state: &mut [T]);

    /// Left Rotation of 1 index.
    fn left_rot(&mut self, state: &mut [T]);

    /// Applies a round.
    ///
    /// One round consists of:
    /// - Computation of the function **F** over `state[0]`;
    /// - **ERF** (Expanding Round Function): propagation of the function F
    ///   result over `state[1..]`;
    /// - **Left rotation** by 1 index of the state.
    fn apply_round(&mut self, round: usize, state: &mut [T]) {
        // Computation of F
        self.f(round, &mut state[0]);

        // Applying Expanding Round function
        self.erf(state);

        // Left rotation
        self.left_rot(state);
    }

    /// Applies one GMiMC permutation.
    fn perm(&mut self, state: &mut [T]) {
        // Apply R rounds
        for round in 0..ROUNDS {
            self.apply_round(round, state);
        }
    }
}

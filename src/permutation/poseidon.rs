// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

//! This module contains an implementation of the `Hades252` permutation
//! algorithm specifically designed to work outside of Rank 1 Constraint Systems
//! (R1CS) or other custom Constraint Systems such as Add/Mul/Custom plonk
//! gate-circuits.
//!
//! The inputs of the permutation function have to be explicitly over the
//! scalar Field of the bls12_381 curve so over a modulus
//! `p = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001`.

use dusk_bls12_381::BlsScalar;
use lazy_static::lazy_static;

/// State width t.
/// It can be 2, 3, 4, 5, 6, 8, 12, 16, 20, 24.
/// Default is 5.
pub const WIDTH: usize = 8;

/// Number of full rounds.
const FULL_ROUNDS: usize = 8;

/// Number of partial rounds.
const PARTIAL_ROUNDS: usize = {
    match WIDTH {
        2 | 3 => 57,
        4 | 5 | 6 | 8 | 12 | 16 | 20 | 24 => 60,
        _ => panic!("Invalid width"),
    }
};

/// Total number of rounds.
pub(crate) const ROUNDS: usize = FULL_ROUNDS + PARTIAL_ROUNDS;

/// Hades parameters generation.
mod params;
use params::{gen_const, gen_mds};

lazy_static! {
    /// Round constants.
    pub(crate) static ref ROUND_CONSTANTS: [[BlsScalar; WIDTH]; ROUNDS] = gen_const();

    /// MDS square matrix.
    pub(crate) static ref MDS_MATRIX: [[BlsScalar; WIDTH]; WIDTH] = gen_mds();
}

/// Hades permutation struct operating in a plonk-circuit.
#[cfg(feature = "zk")]
pub(crate) mod gadget;

/// Hades permutation struct operating on [`BlsScalar`].
pub(crate) mod scalar;

/// Defines the Hades252 permutation algorithm.
///
/// This permutation is a 3-step process that:
/// - Applies half of the `FULL_ROUNDS` (which can be understood as linear ops).
/// - Applies the `PARTIAL_ROUNDS` (which can be understood as non-linear ops).
/// - Applies the other half of the `FULL_ROUNDS`.
///
/// This structure allows to minimize the number of non-linear ops while
/// maintaining the security.
pub(crate) trait Hades<T> {
    /// Add round constants to the state.
    ///
    /// This constants addition, also known as `ARC`, is used to reach
    /// `Confusion and Diffusion` properties for the algorithm.
    ///
    /// Basically it allows to destroy any connection between the inputs and the
    /// outputs of the function.
    fn add_round_constants(&mut self, round: usize, state: &mut [T; WIDTH]);

    /// Computes `input ^ 5 (mod p)`
    ///
    /// The modulo depends on the input you use. In our case the modulo is done
    /// in respect of the scalar field of the bls12_381 curve
    /// `p = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001`.
    fn quintic_s_box(&mut self, value: &mut T);

    /// Multiply the MDS matrix with the state.
    fn mul_matrix(&mut self, state: &mut [T; WIDTH]);

    /// Applies a `Partial Round` also known as a `Partial S-Box layer` to a set
    /// of inputs.
    ///
    /// One partial round consists of 3 steps:
    /// - ARC: Add round constants to the elements of the state.
    /// - Sub Words: Apply `quintic S-Box` just to **the last element of the
    ///   state** generated from the first step.
    /// - Mix Layer: Multiplies the output state from the second step by the
    ///   `MDS_MATRIX`.
    fn apply_partial_round(&mut self, round: usize, state: &mut [T; WIDTH]) {
        // Add round constants to each state element
        self.add_round_constants(round, state);

        // Then apply quintic s-box to the last element of the state
        self.quintic_s_box(&mut state[WIDTH - 1]);

        // Multiply this result by the MDS matrix
        self.mul_matrix(state);
    }

    /// Applies a `Full Round` also known as a `Full S-Box layer` to a set of
    /// inputs.
    ///
    /// One full round consists of 3 steps:
    /// - ARC: Add round constants to the elements of the state.
    /// - Sub Words: Apply `quintic S-Box` to **all of the state-elements**
    ///   generated from the first step.
    /// - Mix Layer: Multiplies the output state from the second step by the
    ///   `MDS_MATRIX`.
    fn apply_full_round(&mut self, round: usize, state: &mut [T; WIDTH]) {
        // Add round constants to each state element
        self.add_round_constants(round, state);

        // Then apply quintic s-box to each element of the state
        state.iter_mut().for_each(|w| self.quintic_s_box(w));

        // Multiply this result by the MDS matrix
        self.mul_matrix(state);
    }

    /// Applies one Hades permutation.
    ///
    /// This permutation is a 3-step process that:
    /// - Applies half of the `FULL_ROUNDS` (which can be understood as linear
    ///   ops).
    /// - Applies the `PARTIAL_ROUNDS` (which can be understood as non-linear
    ///   ops).
    /// - Applies the other half of the `FULL_ROUNDS`.
    ///
    /// This structure allows to minimize the number of non-linear ops while
    /// maintaining the security.
    fn perm(&mut self, state: &mut [T; WIDTH]) {
        // Apply R_f full rounds
        for round in 0..FULL_ROUNDS / 2 {
            self.apply_full_round(round, state);
        }

        // Apply R_P partial rounds
        for round in 0..PARTIAL_ROUNDS {
            self.apply_partial_round(round + FULL_ROUNDS / 2, state);
        }

        // Apply R_f full rounds
        for round in 0..FULL_ROUNDS / 2 {
            self.apply_full_round(
                round + FULL_ROUNDS / 2 + PARTIAL_ROUNDS,
                state,
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use dusk_bls12_381::BlsScalar;

    use super::ROUND_CONSTANTS;

    #[test]
    fn round_constants() {
        // Check each element is non-zero
        let zero = BlsScalar::zero();
        let has_zero = ROUND_CONSTANTS.iter().flatten().any(|&x| x == zero);
        for ctant in ROUND_CONSTANTS.iter().flatten() {
            let bytes = ctant.to_bytes();
            assert!(&BlsScalar::from_bytes(&bytes).unwrap() == ctant);
        }
        assert!(!has_zero);
    }
}

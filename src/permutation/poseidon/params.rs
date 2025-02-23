use crate::permutation::{BlsScalar, cauchy_matrix};
use crate::poseidon::{ROUNDS, WIDTH};

use sha2::{Digest, Sha512};

/// Poseidon round constants generation function based on Sha512 with a fixed
/// seed.
pub(crate) fn gen_const() -> [[BlsScalar; WIDTH]; ROUNDS] {
    let mut constants = [[BlsScalar::zero(); WIDTH]; ROUNDS];
    let mut p = BlsScalar::one();
    let mut bytes = b"poseidon-for-plonk".to_vec();

    for i in 0..ROUNDS {
        for j in 0..WIDTH {
            let mut hasher = Sha512::new();
            Digest::update(&mut hasher, bytes.as_slice());
            bytes = hasher.finalize().to_vec();

            let mut v = [0x00u8; 64];
            v.copy_from_slice(&bytes[0..64]);

            constants[i][j] = BlsScalar::from_bytes_wide(&v) + p;
            p = constants[i][j];
        }
    }

    constants
}

/// Poseidon's Cauchy MDS matrix generation function.
pub(crate) fn gen_mds() -> [[BlsScalar; WIDTH]; WIDTH] {
    let mut mds = [[BlsScalar::zero(); WIDTH]; WIDTH];
    cauchy_matrix(&mut mds);
    mds
}

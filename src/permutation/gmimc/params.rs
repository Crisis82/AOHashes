use crate::gmimc::ROUNDS;
use crate::permutation::{BlsScalar, Shake256PRNG};

/// GMiMC Constant generation function using SHAKE256 with fixed seed.
pub(crate) fn gen_const() -> [BlsScalar; ROUNDS] {
    let mut constants = [BlsScalar::zero(); ROUNDS];

    // SHAKE256 PRNG with fixed seed
    let seed = b"gmimc_seed";
    let mut prng = Shake256PRNG::new(seed);

    // Generating the ROUND CONSTANTS
    for i in 0..ROUNDS {
        constants[i] = prng.next_scalar();
    }

    constants
}

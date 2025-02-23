use dusk_bls12_381::BlsScalar;

cfg_if::cfg_if! {
    if #[cfg(feature = "anemoi")] {
        /// Anemoi
        pub mod anemoi;
        pub use anemoi::WIDTH;
        pub(crate) use anemoi::scalar::ScalarPermutation;
        #[cfg(feature = "zk")]
        pub(crate) use anemoi::gadget::GadgetPermutation;
    } else if #[cfg(feature = "arion")] {
        /// Arion
        pub mod arion;
        pub use arion::WIDTH;
        pub(crate) use arion::scalar::ScalarPermutation;
        #[cfg(feature = "zk")]
        pub(crate) use arion::gadget::GadgetPermutation;
    } else if #[cfg(feature = "gmimc")] {
        /// GMiMC
        pub mod gmimc;
        pub use gmimc::WIDTH;
        pub(crate) use gmimc::scalar::ScalarPermutation;
        #[cfg(feature = "zk")]
        pub(crate) use gmimc::gadget::GadgetPermutation;
    }else if #[cfg(feature = "griffin")] {
        /// Griffin
        pub mod griffin;
        pub use griffin::WIDTH;
        pub(crate) use griffin::scalar::ScalarPermutation;
        #[cfg(feature = "zk")]
        pub(crate) use griffin::gadget::GadgetPermutation;
    } else if #[cfg(feature = "poseidon")] {
        /// Poseidon
        pub mod poseidon;
        pub use poseidon::WIDTH;
        pub(crate) use poseidon::scalar::ScalarPermutation;
        #[cfg(feature = "zk")]
        pub(crate) use poseidon::gadget::GadgetPermutation;
    }  else if #[cfg(feature = "rescue")] {
        /// Rescue
        pub mod rescue;
        pub use rescue::WIDTH;
        pub(crate) use rescue::scalar::ScalarPermutation;
        #[cfg(feature = "zk")]
        pub(crate) use rescue::gadget::GadgetPermutation;
    } else if #[cfg(feature = "rescue_prime")] {
        /// Rescue-Prime
        pub mod rescue_prime;
        pub use rescue_prime::WIDTH;
        pub(crate) use rescue_prime::scalar::ScalarPermutation;
        #[cfg(feature = "zk")]
        pub(crate) use rescue_prime::gadget::GadgetPermutation;
    }
}

cfg_if::cfg_if! {
    if #[cfg(any(feature = "anemoi", feature = "poseidon"))] {
        /// Produces a cauchy matrix
        fn cauchy_matrix<const WIDTH: usize>(matrix: &mut [[BlsScalar; WIDTH]; WIDTH]) {
            let mut xs = [BlsScalar::zero(); WIDTH];
            let mut ys = [BlsScalar::zero(); WIDTH];
            (0..WIDTH).for_each(|i| {
                xs[i] = BlsScalar::from(i as u64);
                ys[i] = BlsScalar::from((i + WIDTH) as u64);
            });

            let mut m = 0;
            (0..WIDTH).for_each(|i| {
                (0..WIDTH).for_each(|j| {
                    matrix[m][j] = (xs[i] + ys[j]).invert().unwrap();
                });
                m += 1;
            });
        }
    }
}

cfg_if::cfg_if! {
    if #[cfg(any(feature = "gmimc", feature = "griffin", feature = "rescue", feature = "rescue_prime"))] {
        use sha3::{
            Shake256,
            digest::{ExtendableOutput, Update, XofReader},
        };

        /// PRNG based on SHAKE256
        struct Shake256PRNG {
            reader: sha3::Shake256Reader,
        }

        impl Shake256PRNG {
            fn new(seed: &[u8]) -> Self {
                let mut hasher = Shake256::default();
                hasher.update(seed);
                let reader = hasher.finalize_xof();
                Shake256PRNG { reader }
            }

            fn next_scalar(&mut self) -> BlsScalar {
                let mut buffer = [0u8; 64];
                self.reader.read(&mut buffer);
                BlsScalar::from_bytes_wide(&buffer)
            }
        }
    }
}

cfg_if::cfg_if! {
    if #[cfg(any(feature = "anemoi", feature = "arion", feature = "griffin", feature = "rescue", feature = "rescue_prime"))] {
        use ethnum::{i256, u256};

        /// Extended Euclidean algorithm
        fn xgcd(a: i256, b: i256) -> (i256, i256, i256) {
            if a == i256::new(0i128) {
                return (b, i256::new(0i128), i256::new(1i128));
            }

            let (gcd, x1, y1) = xgcd(b % a, a);

            let x = y1 - (b / a) * x1;
            let y = x1;

            return (gcd, x, y);
        }

        /// u256 to u64 conversion
        fn u256_to_u64(mut n: u256) -> [u64; 4] {
            // 2^64
            let base = u256::new(1u128 << 64);

            let mut b_adic_expansion = [0u64; 4];
            b_adic_expansion[0] = (n % base).as_u64();
            n = (n - u256::new(b_adic_expansion[0] as u128)) / base;

            let mut i = 1;
            while n != u256::new(0u128) {
                b_adic_expansion[i] = (n % base).as_u64();
                n = (n - u256::new(b_adic_expansion[i] as u128)) / base;
                i += 1;
            }

            b_adic_expansion
        }

        /// Modular inverse of a u64 number returned as [u64; 4] representation
        fn inverse_mod(n: u64) -> [u64; 4] {
            let p: u256 = u256::from_str_hex(
                "0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001",
            )
            .expect("Invalid hex string");

            // cast n to i256
            let exp = i256::new(n as i128);

            // compute Exteded Euclidean Algorithm
            let (_gcd, x, _y) = xgcd(exp, (p - 1).as_i256());
            // cast result to u256
            let d2_inv = x.as_u256();

            // convert u256 to [u64; 4] internal representation
            u256_to_u64(d2_inv)
        }
    }
}

#[cfg(test)]
mod tests {
    use dusk_bls12_381::BlsScalar;
    use dusk_safe::{Call, Safe, Sponge};

    use crate::hash::{Domain, Hash};
    use crate::permutation::{ScalarPermutation, WIDTH};

    extern crate std;
    use std::vec;

    use ff::Field;
    use rand::rngs::OsRng;

    const TESTS: usize = 4;

    #[derive(Default, Debug, Clone, Copy, PartialEq)]
    struct Test();

    impl Safe<BlsScalar, { WIDTH }> for Test {
        // apply permutation
        fn permute(&mut self, state: &mut [BlsScalar; WIDTH]) {
            ScalarPermutation::new().permute(state);
        }

        fn tag(&mut self, input: &[u8]) -> BlsScalar {
            // let _ = input;
            // BlsScalar::zero()
            BlsScalar::hash_to_scalar(input)
        }

        fn add(&mut self, right: &BlsScalar, left: &BlsScalar) -> BlsScalar {
            right + left
        }
    }

    impl Test {
        pub fn new() -> Self {
            Self()
        }
    }

    fn create_hash(input: &[BlsScalar]) -> BlsScalar {
        // let iopattern = vec![Call::Absorb(input.len()), Call::Absorb(1),
        // Call::Squeeze(1)];
        #[allow(unused_unsafe)]
        let iopattern =
            unsafe { vec![Call::Absorb(input.len()), Call::Squeeze(1)] };

        let domain_sep = 0;
        let mut sponge = Sponge::start(Test::new(), iopattern, domain_sep)
            .expect("IO pattern should be valid");
        // absorb given input
        sponge
            .absorb(input.len(), input)
            .expect("Absorption of the input should work fine");

        // absorb padding of one BlsScalar::one()
        // sponge
        //     .absorb(1, &[BlsScalar::one()])
        //     .expect("Absorption of padding should work fine");

        // squeeze single output
        sponge.squeeze(1).expect("Squeezing should work fine");
        let output = sponge.finish().expect("Finish should work fine");
        output[0]
    }

    #[test]
    fn hash() {
        let mut inputs = [BlsScalar::zero(); TESTS];
        for scalar in inputs.iter_mut() {
            *scalar = BlsScalar::random(&mut OsRng);
        }

        let mut outputs = [BlsScalar::zero(); TESTS];
        for (i, scalar) in outputs.iter_mut().enumerate() {
            *scalar = Hash::digest(Domain::Other, &[inputs[i]])[0];
        }

        for i in 0..TESTS {
            assert_eq!(create_hash(&[inputs[i]]), outputs[i]);
        }
    }
}

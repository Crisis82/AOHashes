#![no_std]
#![doc = include_str!("../README.md")]
#![deny(missing_docs)]

extern crate alloc;

mod error;
pub use error::Error;

mod permutation;
pub use permutation::WIDTH;

cfg_if::cfg_if! {
    if #[cfg(feature = "anemoi")] {
        pub use permutation::anemoi;
    } else if #[cfg(feature = "arion")] {
        pub use permutation::arion;
    } else if #[cfg(feature = "gmimc")] {
        pub use permutation::gmimc;
    } else if #[cfg(feature = "griffin")] {
        pub use permutation::griffin;
    } else if #[cfg(feature = "poseidon")] {
        pub use permutation::poseidon;
    } else if #[cfg(feature = "rescue")] {
        pub use permutation::rescue;
    } else if #[cfg(feature = "rescue_prime")] {
        pub use permutation::rescue_prime;
    }
}

mod hash;
#[cfg(feature = "zk")]
pub use hash::gadget::HashGadget;
pub use hash::{Domain, Hash};

cfg_if::cfg_if! {
    if #[cfg(feature = "encryption")] {
        mod encryption;
        pub use encryption::{decrypt, encrypt};
        #[cfg(feature = "zk")]
        pub use encryption::gadget::{decrypt_gadget, encrypt_gadget};
    }
}

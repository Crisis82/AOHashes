use aohashes::{Domain, Hash};
use dusk_bls12_381::BlsScalar;

use ff::Field;
use rand::SeedableRng;
use rand::rngs::StdRng;

fn hash_example() {
    // generate random input
    let mut rng = StdRng::seed_from_u64(0xbeef);
    // let mut rng = StdRng::from_entropy();
    let mut input = [BlsScalar::zero(); 4];
    for scalar in input.iter_mut() {
        *scalar = BlsScalar::random(&mut rng);
    }
    // println!("Input: {:?}", input);

    // digest the input all at once
    let hash = Hash::digest(Domain::Other, &input);

    // update the input gradually
    let mut hasher = Hash::new(Domain::Other);
    hasher.update(&input[..1]);
    hasher.update(&input[1..]);
    assert_eq!(hash, hasher.finalize());
    println!("\nEqual hashes!");
    println!("Hash1: {:?}", hash);
    println!("Hash2: {:?}", hasher.finalize());

    // create a hash used for merkle tree hashing with arity = 4
    let merkle_hash = Hash::digest(Domain::Merkle4, &input);
    println!("\nMerkle Hash: {:?}", merkle_hash);

    // which is different when another domain is used
    assert_ne!(merkle_hash, Hash::digest(Domain::Other, &input));
}

fn main() {
    hash_example();
}

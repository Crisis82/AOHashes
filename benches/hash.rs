use aohashes::{Domain, Hash, HashGadget};
use core::time::Duration;
use criterion::{Criterion, black_box, criterion_group, criterion_main};
use dusk_plonk::prelude::*;
use ff::Field;
use rand::SeedableRng;
use rand::rngs::StdRng;

cfg_if::cfg_if! {
    if #[cfg(any(feature = "anemoi", feature = "arion", feature = "griffin", feature = "poseidon", feature = "rescue", feature = "rescue_prime"))] {
        const CAPACITY: usize = 11;
    } else if #[cfg(feature = "gmimc")] {
        const CAPACITY: usize = 16;
    }
}

#[derive(Default)]
struct SpongeCircuit {
    message: [BlsScalar; 4],
    output: BlsScalar,
}

impl SpongeCircuit {
    pub fn new(message: [BlsScalar; 4], output: BlsScalar) -> Self {
        SpongeCircuit { message, output }
    }
}

impl Circuit for SpongeCircuit {
    fn circuit(&self, composer: &mut Composer) -> Result<(), Error> {
        let mut w_message = [Composer::ZERO; 4];
        w_message
            .iter_mut()
            .zip(self.message)
            .for_each(|(witness, scalar)| {
                *witness = composer.append_witness(scalar);
            });

        let output_witness =
            HashGadget::digest(composer, Domain::Merkle4, &w_message);
        composer.assert_equal_constant(output_witness[0], 0, Some(self.output));

        Ok(())
    }
}

// Benchmark for running sponge on 5 BlsScalar, one permutation
fn bench_sponge(c: &mut Criterion) {
    // Prepare benchmarks and initialize variables
    let label = b"sponge benchmark";
    let mut rng = StdRng::seed_from_u64(0xc10d);
    let pp = PublicParameters::setup(1 << CAPACITY, &mut rng).unwrap();
    let (prover, verifier) = Compiler::compile::<SpongeCircuit>(&pp, label)
        .expect("Circuit should compile successfully");
    let mut proof = Proof::default();
    let message = [BlsScalar::random(&mut rng); 4];
    let public_inputs = Hash::digest(Domain::Merkle4, &message);
    let circuit = SpongeCircuit::new(message, public_inputs[0]);

    // Benchmark sponge native
    c.bench_function("hash 4 BlsScalar", |b| {
        b.iter(|| {
            let _ = Hash::digest(Domain::Merkle4, black_box(&circuit.message));
        })
    });

    // Benchmark proof creation
    c.bench_function("hash 4 BlsScalar proof generation", |b| {
        b.iter(|| {
            (proof, _) = prover
                .prove(&mut rng, black_box(&circuit))
                .expect("Proof generation should succeed");
        })
    });

    // Benchmark proof verification
    c.bench_function("hash 4 BlsScalar proof verification", |b| {
        b.iter(|| {
            verifier
                .verify(black_box(&proof), &public_inputs)
                .expect("Proof verification should succeed");
        })
    });
}

criterion_group! {
    name = benches;
    config = Criterion::default().sample_size(10);
    // usually 100s are fine, maybe gmimc needs 2000s
    // config = Criterion::default().sample_size(100).measurement_time(Duration::from_secs(100));
    targets = bench_sponge
}
criterion_main!(benches);

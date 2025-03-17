use aohashes::{decrypt, decrypt_gadget, encrypt};
// use core::time::Duration;
use criterion::{Criterion, black_box, criterion_group, criterion_main};
use dusk_bls12_381::BlsScalar;
use dusk_jubjub::{GENERATOR_EXTENDED, JubJubAffine, JubJubScalar};
use dusk_plonk::prelude::Error as PlonkError;
use dusk_plonk::prelude::*;
use ff::Field;
use once_cell::sync::Lazy;
use rand::SeedableRng;
use rand::rngs::StdRng;

cfg_if::cfg_if! {
    if #[cfg(any(feature = "anemoi", feature = "arion", feature = "griffin", feature = "poseidon", feature = "rescue", feature = "rescue_prime"))] {
        const CAPACITY: usize = 11;
    } else if #[cfg(feature = "gmimc")] {
        const CAPACITY: usize = 16;
    }
}

const MESSAGE_LEN: usize = 2;

static PUB_PARAMS: Lazy<PublicParameters> = Lazy::new(|| {
    let mut rng = StdRng::seed_from_u64(0xfab);

    PublicParameters::setup(1 << CAPACITY, &mut rng)
        .expect("Setup of public params should pass")
});
static LABEL: &[u8] = b"hash-gadget-tester";

#[derive(Debug)]
struct DecryptionCircuit {
    pub cipher: Vec<BlsScalar>,
    pub shared_secret: JubJubAffine,
    pub nonce: BlsScalar,
}

impl DecryptionCircuit {
    pub fn random(rng: &mut StdRng) -> Self {
        let mut message = [BlsScalar::zero(); MESSAGE_LEN];
        message
            .iter_mut()
            .for_each(|s| *s = BlsScalar::random(&mut *rng));
        let shared_secret =
            GENERATOR_EXTENDED * &JubJubScalar::random(&mut *rng);
        let shared_secret = shared_secret.into();
        let nonce = BlsScalar::random(&mut *rng);
        let cipher = encrypt(&message, &shared_secret, &nonce)
            .expect("encryption should not fail");

        Self {
            cipher,
            shared_secret,
            nonce,
        }
    }
}

impl Default for DecryptionCircuit {
    fn default() -> Self {
        let message = [BlsScalar::zero(); MESSAGE_LEN];
        let mut cipher = message.to_vec();
        cipher.push(BlsScalar::zero());
        let shared_secret = JubJubAffine::identity();
        let nonce = BlsScalar::zero();

        Self {
            cipher,
            shared_secret,
            nonce,
        }
    }
}

impl Circuit for DecryptionCircuit {
    fn circuit(&self, composer: &mut Composer) -> Result<(), PlonkError> {
        // append all variables to the circuit
        let mut cipher_wit = Vec::with_capacity(MESSAGE_LEN + 1);
        self.cipher
            .iter()
            .for_each(|c| cipher_wit.push(composer.append_witness(*c)));
        let secret_wit = composer.append_point(self.shared_secret);
        let nonce_wit = composer.append_witness(self.nonce);

        // decrypt the cipher with the gadget
        let _cipher_result =
            decrypt_gadget(composer, &cipher_wit, &secret_wit, &nonce_wit)
                .expect("decryption should pass");

        Ok(())
    }
}

fn bench_decryption(c: &mut Criterion) {
    let mut rng = StdRng::seed_from_u64(0x42424242);

    let (prover, verifier) =
        Compiler::compile::<DecryptionCircuit>(&PUB_PARAMS, LABEL)
            .expect("compilation should pass");

    let circuit: DecryptionCircuit = DecryptionCircuit::random(&mut rng);
    let public_inputs = Vec::new();
    let mut proof = Proof::default();

    // Benchmark native cipher decryption
    c.bench_function("decrypt 2 BlsScalar", |b| {
        b.iter(|| {
            _ = decrypt(
                black_box(&circuit.cipher),
                black_box(&circuit.shared_secret),
                black_box(&circuit.nonce),
            );
        })
    });

    // Benchmark proof creation
    c.bench_function("decrypt 2 BlsScalar proof generation", |b| {
        b.iter(|| {
            (proof, _) = prover
                .prove(&mut rng, black_box(&circuit))
                .expect("Proof generation should succeed");
        })
    });

    // Benchmark proof verification
    c.bench_function("decrypt 2 BlsScalar proof verification", |b| {
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
    // config = Criterion::default().sample_size(100).measurement_time(Duration::from_secs(100));
    targets = bench_decryption
}
criterion_main!(benches);

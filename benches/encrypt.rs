use aohashes::{encrypt, encrypt_gadget};
use core::time::Duration;
use criterion::{Criterion, black_box, criterion_group, criterion_main};
use dusk_bls12_381::BlsScalar;
use dusk_jubjub::{GENERATOR_EXTENDED, JubJubAffine, JubJubScalar};
use dusk_plonk::prelude::Error as PlonkError;
use dusk_plonk::prelude::*;
use ff::Field;
use once_cell::sync::Lazy;
use rand::SeedableRng;
use rand::rngs::StdRng;

// for gmimg and poseidon at least 12 if not more is needed
const CAPACITY: usize = 12;

const MESSAGE_LEN: usize = 2;

static PUB_PARAMS: Lazy<PublicParameters> = Lazy::new(|| {
    let mut rng = StdRng::seed_from_u64(0xfab);

    PublicParameters::setup(1 << CAPACITY, &mut rng)
        .expect("Setup of public params should pass")
});
static LABEL: &[u8] = b"hash-gadget-tester";

#[derive(Debug)]
struct EncryptionCircuit {
    pub message: [BlsScalar; MESSAGE_LEN],
    pub shared_secret: JubJubAffine,
    pub nonce: BlsScalar,
}

impl EncryptionCircuit {
    pub fn random(rng: &mut StdRng) -> Self {
        let mut message = [BlsScalar::zero(); MESSAGE_LEN];
        message
            .iter_mut()
            .for_each(|s| *s = BlsScalar::random(&mut *rng));
        let shared_secret =
            GENERATOR_EXTENDED * &JubJubScalar::random(&mut *rng);
        let nonce = BlsScalar::random(&mut *rng);

        Self {
            message,
            shared_secret: shared_secret.into(),
            nonce,
        }
    }
}

impl Default for EncryptionCircuit {
    fn default() -> Self {
        let message = [BlsScalar::zero(); MESSAGE_LEN];
        let shared_secret = JubJubAffine::identity();
        let nonce = BlsScalar::zero();

        Self {
            message,
            shared_secret,
            nonce,
        }
    }
}

impl Circuit for EncryptionCircuit {
    fn circuit(&self, composer: &mut Composer) -> Result<(), PlonkError> {
        // append all variables to the circuit
        let mut message_wit = [Composer::ZERO; MESSAGE_LEN];
        message_wit
            .iter_mut()
            .zip(self.message)
            .for_each(|(w, m)| *w = composer.append_witness(m));
        let secret_wit = composer.append_point(self.shared_secret);
        let nonce_wit = composer.append_witness(self.nonce);

        // encrypt the message with the gadget
        let _cipher_result =
            encrypt_gadget(composer, &message_wit, &secret_wit, &nonce_wit)
                .expect("encryption should pass");

        Ok(())
    }
}

fn bench_encryption(c: &mut Criterion) {
    let mut rng = StdRng::seed_from_u64(0x42424242);

    let (prover, verifier) =
        Compiler::compile::<EncryptionCircuit>(&PUB_PARAMS, LABEL)
            .expect("compilation should pass");

    let circuit: EncryptionCircuit = EncryptionCircuit::random(&mut rng);
    let public_inputs = Vec::new();
    let mut proof = Proof::default();

    // Benchmark native cipher decryption
    c.bench_function("encrypt 2 BlsScalar", |b| {
        b.iter(|| {
            let _ = encrypt(
                black_box(&circuit.message),
                black_box(&circuit.shared_secret),
                black_box(&circuit.nonce),
            );
        })
    });

    // Benchmark proof creation
    c.bench_function("encrypt 2 BlsScalar proof generation", |b| {
        b.iter(|| {
            (proof, _) = prover
                .prove(&mut rng, black_box(&circuit))
                .expect("Proof generation should succeed");
        })
    });

    // Benchmark proof verification
    c.bench_function("encrypt 2 BlsScalar proof verification", |b| {
        b.iter(|| {
            verifier
                .verify(black_box(&proof), &public_inputs)
                .expect("Proof verification should succeed");
        })
    });
}

criterion_group! {
    name = benches;
    // config = Criterion::default().sample_size(10);
    config = Criterion::default().sample_size(1000).measurement_time(Duration::from_secs(800));
    targets = bench_encryption
}
criterion_main!(benches);

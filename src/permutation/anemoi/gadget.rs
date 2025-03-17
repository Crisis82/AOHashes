use dusk_bls12_381::BlsScalar;
use dusk_plonk::prelude::*;
use dusk_safe::Safe;

use super::AnemoiHash;
use crate::anemoi::{
    ALPHA, ALPHA_INV, G, G_INV, L, MATRIX_X, MATRIX_Y, ROUND_CONSTANTS, WIDTH,
};

/// Requires a reference to a plonk circuit [`Composer`].
pub(crate) struct GadgetPermutation<'a> {
    /// A reference to the constraint system used by the gadgets
    composer: &'a mut Composer,
}

impl<'a> GadgetPermutation<'a> {
    /// Constructs a new `GadgetPermutation` with the constraint system.
    pub fn new(composer: &'a mut Composer) -> Self {
        Self { composer }
    }
}

impl<'a> Safe<Witness, WIDTH> for GadgetPermutation<'a> {
    fn permute(&mut self, state: &mut [Witness; WIDTH]) {
        self.perm(state);
    }

    fn tag(&mut self, input: &[u8]) -> Witness {
        let tag = BlsScalar::hash_to_scalar(input.as_ref());
        // append the tag as a constant
        self.composer.append_constant(tag)
    }

    fn add(&mut self, right: &Witness, left: &Witness) -> Witness {
        let constraint = Constraint::new().left(1).a(*left).right(1).b(*right);
        self.composer.gate_add(constraint)
    }
}

impl<'a> AnemoiHash<Witness> for GadgetPermutation<'a> {
    fn constant_addition(&mut self, state: &mut [Witness], round: usize) {
        for (i, witness) in state.iter_mut().enumerate() {
            let constraint = Constraint::new()
                .left(1)
                .a(*witness)
                .constant(ROUND_CONSTANTS[round][i]);
            *witness = self.composer.gate_add(constraint);
        }
    }

    fn linear_layer(&mut self, state: &mut [Witness]) {
        let mut result = [Composer::ZERO; WIDTH];

        for i in 0..L {
            for j in 0..L {
                let constraint = Constraint::new()
                    .left(1)
                    .a(result[i])
                    .right(MATRIX_X[i][j])
                    .b(state[j]);
                result[i] = self.composer.gate_add(constraint);

                let constraint = Constraint::new()
                    .left(1)
                    .a(result[i + L])
                    .right(MATRIX_Y[i][j])
                    .b(state[j + L]);
                result[i + L] = self.composer.gate_add(constraint);
            }
        }

        state.copy_from_slice(&result);
    }

    fn pht(&mut self, state: &mut [Witness]) {
        // Applies Y <- Y + X
        for i in L..WIDTH {
            let constraint = Constraint::new()
                .left(1)
                .a(state[i])
                .right(1)
                .b(state[i - L]);
            state[i] = self.composer.gate_add(constraint);
        }

        // Applies X <- X + Y
        for i in 0..L {
            let constraint = Constraint::new()
                .left(1)
                .a(state[i])
                .right(1)
                .b(state[i + L]);
            state[i] = self.composer.gate_add(constraint);
        }
    }

    fn sbox_layer(&mut self, state: &mut [Witness]) {
        for i in 0..L {
            // x_i = x_i - g * y_i^2 - g_inv
            let constraint = Constraint::new()
                .mult(-G)
                .a(state[i + L])
                .b(state[i + L])
                .fourth(1)
                .d(state[i])
                .constant(-(*G_INV));
            state[i] = self.composer.gate_add(constraint);

            // x_i^alpha_inv
            let scalar_x = self.composer[state[i]];
            let scalar_y = scalar_x.pow_vartime(&ALPHA_INV);
            let witness_y = self.composer.append_witness(scalar_y);

            // value^2
            let constraint =
                Constraint::new().mult(1).a(witness_y).b(witness_y);
            let power2 = self.composer.gate_mul(constraint);
            // value^4
            let constraint = Constraint::new().mult(1).a(power2).b(power2);
            let power4 = self.composer.gate_mul(constraint);

            let power = if ALPHA[0] == 7 {
                // value^6
                let constraint = Constraint::new().mult(1).a(power4).b(power2);
                self.composer.gate_mul(constraint)
            } else {
                power4
            };

            // value^5 or value^7
            let constraint = Constraint::new().mult(1).a(power).b(witness_y);
            state[i] = self.composer.gate_mul(constraint);

            // y_i = y_i - x_i^alpha_inv
            let constraint = Constraint::new()
                .left(1)
                .a(state[i + L])
                .right(-BlsScalar::from(1))
                .b(state[i]);
            state[i + L] = self.composer.gate_add(constraint);

            state[i] = witness_y;

            // x_i = x_i + g * y_i^2
            let constraint = Constraint::new()
                .mult(G)
                .a(state[i + L])
                .b(state[i + L])
                .fourth(1)
                .d(state[i]);
            state[i] = self.composer.gate_add(constraint);
        }
    }
}

#[cfg(feature = "encryption")]
impl dusk_safe::Encryption<Witness, WIDTH> for GadgetPermutation<'_> {
    fn subtract(&mut self, minuend: &Witness, subtrahend: &Witness) -> Witness {
        let constraint = Constraint::new()
            .left(1)
            .a(*minuend)
            .right(-BlsScalar::one())
            .b(*subtrahend);
        self.composer.gate_add(constraint)
    }

    fn is_equal(&mut self, lhs: &Witness, rhs: &Witness) -> bool {
        self.composer.assert_equal(*lhs, *rhs);
        // for the encryption to work we need to return true here, the proof
        // creation will fail at a later point if the above assertion isn't met
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::permutation::ScalarPermutation;

    use core::result::Result;
    use ff::Field;
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    #[derive(Default)]
    struct TestCircuit {
        i: [BlsScalar; WIDTH],
        o: [BlsScalar; WIDTH],
    }

    impl Circuit for TestCircuit {
        fn circuit(&self, composer: &mut Composer) -> Result<(), Error> {
            let zero = Composer::ZERO;

            let mut perm: [Witness; WIDTH] = [zero; WIDTH];

            let mut i_wit: [Witness; WIDTH] = [zero; WIDTH];
            self.i.iter().zip(i_wit.iter_mut()).for_each(|(i, w)| {
                *w = composer.append_witness(*i);
            });

            let mut o_wit: [Witness; WIDTH] = [zero; WIDTH];
            self.o.iter().zip(o_wit.iter_mut()).for_each(|(o, w)| {
                *w = composer.append_witness(*o);
            });

            // Apply gadget permutation.
            GadgetPermutation::new(composer).permute(&mut i_wit);

            // Copy the result of the permutation into the perm.
            perm.copy_from_slice(&i_wit);

            // Check that the Gadget perm results = BlsScalar perm results
            i_wit.iter().zip(o_wit.iter()).for_each(|(p, o)| {
                composer.assert_equal(*p, *o);
            });

            Ok(())
        }
    }

    /// Generate a random input and perform a permutation
    fn anemoi() -> ([BlsScalar; WIDTH], [BlsScalar; WIDTH]) {
        let mut input = [BlsScalar::zero(); WIDTH];

        let mut rng = StdRng::seed_from_u64(0xbeef);

        input
            .iter_mut()
            .for_each(|s| *s = BlsScalar::random(&mut rng));

        let mut output = [BlsScalar::zero(); WIDTH];

        output.copy_from_slice(&input);
        ScalarPermutation::new().permute(&mut output);

        (input, output)
    }

    /// Setup the test circuit prover and verifier
    fn setup() -> Result<(Prover, Verifier), Error> {
        const CAPACITY: usize = 1 << 12;

        let mut rng = StdRng::seed_from_u64(0xbeef);

        let pp = PublicParameters::setup(CAPACITY, &mut rng)?;
        let label = b"anemoi_gadget_tester";

        Compiler::compile::<TestCircuit>(&pp, label)
    }

    #[test]
    fn preimage() -> Result<(), Error> {
        let (prover, verifier) = setup()?;

        let (i, o) = anemoi();

        let circuit = TestCircuit { i, o };
        let mut rng = StdRng::seed_from_u64(0xbeef);

        // Proving
        let (proof, public_inputs) = prover.prove(&mut rng, &circuit)?;

        // Verifying
        verifier.verify(&proof, &public_inputs)?;

        Ok(())
    }

    #[test]
    fn preimage_constant() -> Result<(), Error> {
        let (prover, verifier) = setup()?;

        // Prepare input & output
        let i = [BlsScalar::from(5000u64); WIDTH];
        let mut o = [BlsScalar::from(5000u64); WIDTH];
        ScalarPermutation::new().permute(&mut o);

        let circuit = TestCircuit { i, o };
        let mut rng = StdRng::seed_from_u64(0xbeef);

        // Proving
        let (proof, public_inputs) = prover.prove(&mut rng, &circuit)?;

        // Verifying
        verifier.verify(&proof, &public_inputs)?;

        Ok(())
    }

    #[test]
    fn preimage_fails() -> Result<(), Error> {
        let (prover, _) = setup()?;

        // Generate [31, 0, 0, 0, 0] as real input to the perm but build the
        // proof with [31, 31, 31, 31, 31]. This should fail on verification
        // since the Proof contains incorrect statements.
        let x_scalar = BlsScalar::from(31u64);

        let mut i = [BlsScalar::zero(); WIDTH];
        i[1] = x_scalar;

        let mut o = [BlsScalar::from(31u64); WIDTH];
        ScalarPermutation::new().permute(&mut o);

        let circuit = TestCircuit { i, o };
        let mut rng = StdRng::seed_from_u64(0xbeef);

        // Proving should fail
        assert!(
            prover.prove(&mut rng, &circuit).is_err(),
            "proving should fail since the circuit is invalid"
        );

        Ok(())
    }
}

use dusk_bls12_381::BlsScalar;
use dusk_plonk::prelude::*;
use dusk_safe::Safe;

use super::GriffinPi;
use crate::griffin::{ALPHA, BETA, D, D_INV, MATRIX, ROUND_CONSTANTS, WIDTH};

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

impl<'a> GriffinPi<Witness> for GadgetPermutation<'a> {
    fn non_linear_layer(&mut self, state: &mut [Witness]) {
        let mut result = [Composer::ZERO; WIDTH];

        // i = 0
        // value^D_INV
        let scalar_x = self.composer[state[0]];
        let scalar_y = scalar_x.pow_vartime(&*D_INV);
        let witness_y = self.composer.append_witness(scalar_y);

        // witness_y^2
        let constraint = Constraint::new().mult(1).a(witness_y).b(witness_y);
        let power2 = self.composer.gate_mul(constraint);
        // witness_y^4
        let constraint = Constraint::new().mult(1).a(power2).b(power2);
        let power4 = self.composer.gate_mul(constraint);
        let power = if D[0] == 7 {
            // witness_y^6
            let constraint = Constraint::new().mult(1).a(power4).b(power2);
            self.composer.gate_mul(constraint)
        } else {
            power4
        };
        // witness_y^5 or witness_y^7
        let constraint = Constraint::new().mult(1).a(power).b(witness_y);
        result[0] = self.composer.gate_mul(constraint);
        result[0] = witness_y;

        // i = 1
        // value^D

        // value^2
        let constraint = Constraint::new().mult(1).a(state[1]).b(state[1]);
        let power2 = self.composer.gate_mul(constraint);
        // value^4
        let constraint = Constraint::new().mult(1).a(power2).b(power2);
        let power4 = self.composer.gate_mul(constraint);
        let power = if D[0] == 7 {
            // value^6
            let constraint = Constraint::new().mult(1).a(power4).b(power2);
            self.composer.gate_mul(constraint)
        } else {
            power4
        };
        let constraint = Constraint::new().mult(1).a(power).b(state[1]);
        result[1] = self.composer.gate_mul(constraint);

        // i >= 2
        for i in 2..WIDTH {
            let l = if i == 2 {
                // gamma_i = i -1 = 2 -1 = 1
                // l_2 = gamma_i * z_0 + z_1 + 0 = z_0 + z_1
                let constraint = Constraint::new()
                    .left(1)
                    .a(result[0])
                    .right(1)
                    .b(result[1]);
                self.composer.gate_add(constraint)
            } else {
                // gamma_i = i -1
                // l_i = gamma_i * z_0 + z_1 + state_(i-1)
                let constraint = Constraint::new()
                    .left(BlsScalar::from((i - 1) as u64))
                    .a(result[0])
                    .right(1)
                    .b(result[1])
                    .fourth(1)
                    .d(state[i - 1]);
                self.composer.gate_add(constraint)
            };

            // z_i = state_i * (l^2 + alpha_i * l + beta_i)
            let constraint = Constraint::new()
                .mult(1)
                .a(l)
                .b(l)
                .right(ALPHA[i])
                .constant(BETA[i]);
            let right_side = self.composer.gate_mul(constraint);
            let constraint =
                Constraint::new().mult(1).a(right_side).b(state[i]);
            result[i] = self.composer.gate_mul(constraint);
        }

        state.copy_from_slice(&result);
    }

    fn linear_layer(&mut self, state: &mut [Witness], _round: usize) {
        let mut result = [Composer::ZERO; WIDTH];

        for i in 0..WIDTH {
            for j in 0..WIDTH {
                let constraint = Constraint::new()
                    .left(1)
                    .a(result[i])
                    .right(MATRIX[i][j])
                    .b(state[j]);
                result[i] = self.composer.gate_add(constraint);
            }
        }

        state.copy_from_slice(&result);
    }

    fn constant_addition(&mut self, state: &mut [Witness], round: usize) {
        for i in 0..WIDTH {
            let constraint = Constraint::new()
                .left(1)
                .a(state[i])
                .constant(ROUND_CONSTANTS[round]);
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
            let mut i_wit: [Witness; WIDTH] = [Composer::ZERO; WIDTH];
            self.i.iter().zip(i_wit.iter_mut()).for_each(|(i, w)| {
                *w = composer.append_witness(*i);
            });

            let mut o_wit: [Witness; WIDTH] = [Composer::ZERO; WIDTH];
            self.o.iter().zip(o_wit.iter_mut()).for_each(|(o, w)| {
                *w = composer.append_witness(*o);
            });

            // Apply gadget permutation.
            GadgetPermutation::new(composer).permute(&mut i_wit);

            // Check that the Gadget perm results = BlsScalar perm results
            i_wit.iter().zip(o_wit.iter()).for_each(|(p, o)| {
                composer.assert_equal(*p, *o);
            });

            Ok(())
        }
    }

    /// Generate a random input and perform a permutation
    fn griffin() -> ([BlsScalar; WIDTH], [BlsScalar; WIDTH]) {
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
        let label = b"griffin_gadget_tester";

        Compiler::compile::<TestCircuit>(&pp, label)
    }

    #[test]
    fn preimage() -> Result<(), Error> {
        let (prover, verifier) = setup()?;

        let (i, o) = griffin();

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

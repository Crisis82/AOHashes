use dusk_bls12_381::BlsScalar;
use dusk_plonk::prelude::*;
use dusk_safe::Safe;

use super::ArionHash;
use crate::arion::{
    AFFINE_CONSTANTS, ALPHA1, ALPHA2, BETA, D1, D2, D2_INV, MATRIX, WIDTH,
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

impl<'a> ArionHash<Witness> for GadgetPermutation<'a> {
    fn gtds(&mut self, round: usize, state: &mut [Witness]) {
        let mut result = [Composer::ZERO; WIDTH];

        // High degree inverse S-Box is inverse of x^D2
        let scalar_x = self.composer[state[WIDTH - 1]];
        let scalar_y = scalar_x.pow_vartime(&*D2_INV);
        let witness_y = self.composer.append_witness(scalar_y);

        let constraint = if D2[0] == 125 {
            // y^2
            let constraint =
                Constraint::new().mult(1).a(witness_y).b(witness_y);
            let power2 = self.composer.gate_mul(constraint);
            // y^4
            let constraint = Constraint::new().mult(1).a(power2).b(power2);
            let power4 = self.composer.gate_mul(constraint);
            // y^5
            let constraint = Constraint::new().mult(1).a(power4).b(witness_y);
            let power5 = self.composer.gate_mul(constraint);
            // y^10
            let constraint = Constraint::new().mult(1).a(power5).b(power5);
            let power10 = self.composer.gate_mul(constraint);
            // y^20
            let constraint = Constraint::new().mult(1).a(power10).b(power10);
            let power20 = self.composer.gate_mul(constraint);
            // y^40
            let constraint = Constraint::new().mult(1).a(power20).b(power20);
            let power40 = self.composer.gate_mul(constraint);
            // y^80
            let constraint = Constraint::new().mult(1).a(power40).b(power40);
            let power80 = self.composer.gate_mul(constraint);
            // y^120
            let constraint = Constraint::new().mult(1).a(power80).b(power40);
            let power120 = self.composer.gate_mul(constraint);
            // y^125
            Constraint::new().mult(1).a(power120).b(power5)
        } else {
            // y^2
            let constraint =
                Constraint::new().mult(1).a(witness_y).b(witness_y);
            let power2 = self.composer.gate_mul(constraint);
            // y^4
            let constraint = Constraint::new().mult(1).a(power2).b(power2);
            let power4 = self.composer.gate_mul(constraint);
            // y^8
            let constraint = Constraint::new().mult(1).a(power4).b(power4);
            let power8 = self.composer.gate_mul(constraint);
            // y^16
            let constraint = Constraint::new().mult(1).a(power8).b(power8);
            let power16 = self.composer.gate_mul(constraint);
            // y^32
            let constraint = Constraint::new().mult(1).a(power16).b(power16);
            let power32 = self.composer.gate_mul(constraint);
            // y^64
            let constraint = Constraint::new().mult(1).a(power32).b(power32);
            let power64 = self.composer.gate_mul(constraint);
            // y^128
            let constraint = Constraint::new().mult(1).a(power64).b(power64);
            let power128 = self.composer.gate_mul(constraint);
            let power = if D2[0] == 161 {
                // y^160
                let constraint =
                    Constraint::new().mult(1).a(power128).b(power32);
                self.composer.gate_mul(constraint)
            } else if D2[0] == 193 {
                // y^192
                let constraint =
                    Constraint::new().mult(1).a(power128).b(power64);
                self.composer.gate_mul(constraint)
            } else {
                // y^256
                let constraint =
                    Constraint::new().mult(1).a(power128).b(power128);
                self.composer.gate_mul(constraint)
            };
            // y^161 or y^193 or y^256 depending on D2
            Constraint::new().mult(1).a(power).b(witness_y)
        };

        result[WIDTH - 1] = self.composer.gate_mul(constraint);

        let constraint = Constraint::new()
            .left(1)
            .a(result[WIDTH - 1])
            .right(1)
            .b(state[WIDTH - 1]);
        let mut s = self.composer.gate_add(constraint);

        // (WIDTH - 1)-th round to 2-nd component
        for i in (1..WIDTH - 1).rev() {
            // D1 S-Box
            // y^2
            let constraint = Constraint::new().mult(1).a(state[i]).b(state[i]);
            let power2 = self.composer.gate_mul(constraint);
            let power = if D1[0] == 5 {
                // y^4
                let constraint = Constraint::new().mult(1).a(power2).b(power2);
                self.composer.gate_mul(constraint)
            } else {
                power2
            };
            // y^3 or y^5
            let constraint = Constraint::new().mult(1).a(power).b(state[i]);
            result[i] = self.composer.gate_mul(constraint);

            // Evaluate g
            let constraint = Constraint::new()
                .mult(1)
                .a(s)
                .b(s)
                .right(ALPHA1[round][i])
                .constant(ALPHA2[round][i]);
            let g = self.composer.gate_add(constraint);

            // Evaluate h
            let constraint =
                Constraint::new().mult(1).a(s).b(s).right(BETA[round][i]);
            let h = self.composer.gate_add(constraint);

            // Evaluate x^d*g + h
            let constraint =
                Constraint::new().mult(1).a(result[i]).b(g).fourth(1).d(h);
            result[i] = self.composer.gate_add(constraint);

            let constraint = Constraint::new()
                .left(1)
                .a(s)
                .right(1)
                .b(state[i])
                .fourth(1)
                .d(result[i]);
            s = self.composer.gate_add(constraint);
        }

        // 1-st component

        // D1 S-Box
        // y^2
        let constraint = Constraint::new().mult(1).a(state[0]).b(state[0]);
        let power = self.composer.gate_mul(constraint);
        let power = if D1[0] == 5 {
            // y^4
            let constraint = Constraint::new().mult(1).a(power).b(power);
            self.composer.gate_mul(constraint)
        } else {
            power
        };
        // y^3 or y^5
        let constraint = Constraint::new().mult(1).a(power).b(state[0]);
        result[0] = self.composer.gate_mul(constraint);

        // Evaluate g
        let constraint = Constraint::new()
            .mult(1)
            .a(s)
            .b(s)
            .right(ALPHA1[round][0])
            .constant(ALPHA2[round][0]);
        let g = self.composer.gate_add(constraint);

        // Evaluate h
        let constraint =
            Constraint::new().mult(1).a(s).b(s).right(BETA[round][0]);
        let h = self.composer.gate_add(constraint);

        // Evaluate x^d*g + h
        let constraint =
            Constraint::new().mult(1).a(result[0]).b(g).fourth(1).d(h);
        result[0] = self.composer.gate_add(constraint);

        state.copy_from_slice(&result);
    }

    fn mul_matrix(&mut self, state: &mut [Witness]) {
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

    fn affine_layer(&mut self, round: usize, state: &mut [Witness]) {
        let mut result = [Composer::ZERO; WIDTH];
        result.copy_from_slice(&state);

        self.mul_matrix(&mut result);

        for i in 0..WIDTH {
            let constraint = Constraint::new()
                .left(1)
                .a(result[i])
                .constant(AFFINE_CONSTANTS[round][i]);
            result[i] = self.composer.gate_add(constraint);
        }

        state.copy_from_slice(&result);
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
    fn arion() -> ([BlsScalar; WIDTH], [BlsScalar; WIDTH]) {
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
        let label = b"arion_gadget_tester";

        Compiler::compile::<TestCircuit>(&pp, label)
    }

    #[test]
    fn preimage() -> Result<(), Error> {
        let (prover, verifier) = setup()?;

        let (i, o) = arion();

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

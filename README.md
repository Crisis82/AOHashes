# AOHashes

Pure Rust implementation of several Arimethization-Oriented primitives, providing two working modes:
- The **normal** hashing functionalities operating on the pairing Elliptic Curve `bls12-381`, re-implemented by the Dusk company and shipped via the crate [dusk-bls12_381](https://github.com/dusk-network/bls12_381).
- The **gadget** hashing functionalities that build a _Plonk_ circuit, based on it's Dusk implementation ([Dusk-Plonk](https://github.com/dusk-network/plonk)).

The primitives implemented are the following:
- [Anemoi](https://eprint.iacr.org/2022/840)
- [Arion](https://arxiv.org/abs/2303.04639)
- [GMiMC](https://eprint.iacr.org/2019/397)
- [Griffin](https://eprint.iacr.org/2022/403)
- [Poseidon](https://eprint.iacr.org/2019/458)
- [Rescue](https://eprint.iacr.org/2019/426)
- [Rescue-Prime](https://eprint.iacr.org/2020/1143)

Each primitive is linked to its original documentation source, plus for redundancy has been stored a copy of each paper under the [docs](./docs/) folder.

## Usage

To simplify the usage and exploiting the cargo functionalities, for each function has been implemented a cargo features that enables all the inherent methods. Plus are also privided two additional features, which are `encryption` and `zk`, to specify the working mode both for the _tests_ and the _benchmarks_.
The `encryption` feature enables hashing, encrypting and decrypting on bls12-381, while `zk` enables the zero-knowledge proof over a Plonk circuit.

### Simple run

To check the hashing functionalities of an primitive, the simplest command that can be run is the following:
```shell
cargo test --features={hash_name}
```

where the available `{hash_name}` are:
- **anemoi**: to specify the `Anemoi` hashing
- **arion**: to specify the `Arion` hashing
- **gmimc**: to specify the `GMiMC` hashing
- **griffin**: to specify the `Griffin` hashing
- **poseidon**: to specify the `Poseidon` hashing
- **rescue**: to specify the `Rescue` hashing
- **rescue_prime**: to specify the `Rescue-Prime` hashing

therefore, the following shell code
```shell
cargo test --features=gmimc
```
will execute `GMiMC`.

> [!CAUTION]  
> For simplicity, this library has been organized in such a way that only an hashing algorithms can be run at a time, thus if are provided more than 1 primitive features all together, it will run only the first, in alphabetical order, `{hash_name}` among the ones provided.

### Example

To run a simple hash example that consumes the same input in two different ways and compare its final digest, run:
```shell
cargo run --example hash_example --features={hash_name}
```

### Test

To execute the tests on the **normal** functionalities, run:
```shell
cargo test --features={hash_name},encryption
```
Instead, to test the **gadget** functionalities, use the command:
```shell
cargo test --features={hash_name},zk
```

Depending if we want to test the primitive on scalar or Plonk circuit, you can alternate between the two features, otherwise `encryption` and `zk`, can also be run altogheter with the following command:
```shell
cargo test --features={hash_name},zk,encryption
```

### Benchmarks

To perform some benchmarks, the commands and options are the same explained for the tests, except that we subsitute the cargo command `test` with `bench`, like that:
```shell
cargo bench --features={hash_name},encryption
```

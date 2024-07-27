# STARK101 Rust

This repository contains rust version of [Stark 101](https://starkware.co/stark-101/) by Starkware.

I have added comments here and there explaining the code, but I recommend going through the start101 videos, slides and github for better understanding.

## Structure

- [main.rs](./src/main.rs) - Code for executing prover and verifier.
- [utils.rs](./src/utils.rs) - All the useful functions are here.
- [field](./src/field/) - A minimal module for filed operations.
- [polynomial](./src/polynomial/) - A minimal module for polynomials in a field.
- [merkle tree](./src/merkle_tree/) - A wrapper around, rs_merkle crate.
- [channel](./src/channel/) - Simulates prover and verifier inteactions.
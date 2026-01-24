# Spartan-BN254

A Rust implementation of the [Spartan](https://eprint.iacr.org/2019/550) zkSNARK, adapted for the **BN254** curve for compatibility with Circom circuits.

## Features

- **BN254 Curve**: Compatible with Circom-generated R1CS files
- **Two PCS Modes**:
  - **Hyrax** (default): No trusted setup, O(√n) proof size
  - **KZG**: Requires trusted setup, O(1) proof size (smaller proofs)
- **Efficient Prover**: Optimized MSM and sumcheck implementations
- **Single-Threaded**: Deterministic benchmarking

## Quick Start

### Prerequisites

1. Rust 1.70+ installed
2. Circom circuit compiled to `.r1cs` format
3. Witness file in `.wtns` format

### Benchmark the Keyless Circuit

```bash
# Clone and build
cd spartan-bn254

# Run with Hyrax (default, no trusted setup)
RAYON_NUM_THREADS=1 cargo run --release --example keyless_benchmark

# Run with KZG (smaller proofs, requires trusted setup)
RAYON_NUM_THREADS=1 cargo run --release --features kzg --example keyless_benchmark
```

**Note**: The `RAYON_NUM_THREADS=1` ensures single-threaded execution for reproducible benchmarks.

## Expected Output

The benchmark produces a detailed prover time breakdown and proof size breakdown:

```
╔══════════════════════════════════════════════════════════════════════════╗
║   SPARTAN-BN254 KEYLESS BENCHMARK (HYRAX MODE)                           ║
╚══════════════════════════════════════════════════════════════════════════╝

  PROVER TIME BREAKDOWN (HYRAX)
  ──────────────────────────────
  R1CS SATISFIABILITY PROOF
    Sumcheck + Witness Opening:            3.45s  (  1.7%)
  
  R1CS EVAL PROOF (Lookup Argument)
    [a] EqPolynomial evaluation:           0.10s  (  0.0%)
    [b] Derefs computation:                0.14s  (  0.1%)
    [c] Derefs commitment (Hyrax MSM):   166.20s  ( 79.6%) ◀━━
    [d] Network construction:              4.07s  (  1.9%)
    [e] Network proof (Hyrax openings):   34.45s  ( 16.5%)
  ─────────────────────────────────────────────────
  TOTAL PROVE TIME:                      208.80s  (100.0%)

  SUMMARY (HYRAX)
  ───────────────
  Encode Time:           60.7s
  Prove Time:           208.8s
  Verify Time:          390.1ms
  Proof Size:           246.4 KB
  Trusted Setup:        No
```

## Comparison: Hyrax vs KZG

| Metric | Hyrax | KZG |
|--------|-------|-----|
| **Trusted Setup** | No | Yes |
| **Proof Size** | 246.4 KB | 117.6 KB |
| **Prove Time** | 208.8s | 264.6s |
| **Verify Time** | 390.1ms | 289.5ms |
| **Main Bottleneck** | Derefs Commitment (80%) | Network Proof (59%) |

## Project Structure

```
spartan-bn254/
├── src/
│   ├── lib.rs              # Library entry point
│   ├── snark.rs            # Top-level SNARK (prove/verify)
│   ├── r1cs.rs             # R1CS structures and commitment
│   ├── r1csproof.rs        # R1CS satisfiability proof
│   ├── hyrax.rs            # Hyrax PCS (dense multilinear polynomial)
│   ├── sparse_mlpoly_full.rs # Sparse matrix polynomials (lookup argument)
│   ├── kzg.rs              # KZG polynomial commitment scheme
│   ├── sumcheck.rs         # Sumcheck protocol
│   └── ...
├── examples/
│   └── keyless_benchmark.rs # Main benchmark (supports both modes)
├── Cargo.toml
└── README.md
```

## Configuration

### Feature Flags

| Flag | Description |
|------|-------------|
| `default` | Enables `parallel` and `serde` |
| `parallel` | Multi-threaded proving (uses Rayon) |
| `serde` | Serialization support for proofs |
| `kzg` | Use KZG instead of Hyrax for derefs commitment |

### Environment Variables

| Variable | Description |
|----------|-------------|
| `RAYON_NUM_THREADS=1` | Single-threaded mode for benchmarks |

## License

MIT License

## References

- [Spartan Paper](https://eprint.iacr.org/2019/550) by Srinath Setty
- [Original Spartan Implementation](https://github.com/microsoft/Spartan)
- [arkworks-spartan](https://github.com/arkworks-rs/spartan) - Arkworks port

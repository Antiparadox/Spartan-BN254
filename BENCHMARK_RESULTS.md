# Spartan-BN254 Benchmark Results

Benchmarks for the Aptos Keyless circuit on Apple M2 Max (single-threaded).

## Circuit: Keyless (Aptos)

| Parameter | Value |
|-----------|-------|
| Constraints | 1,040,083 |
| Variables | 1,016,724 |
| Public Inputs | 1 |
| Total NNZ | 7,132,133 |
| Padded Size | 2^20 × 2^20 |

---

## Comparison: Hyrax vs KZG

| Metric | Hyrax (Default) | KZG | Winner |
|--------|-----------------|-----|--------|
| **Trusted Setup** | No | Yes | Hyrax |
| **Encode Time** | 60.7s | 61.8s | ~same |
| **Prove Time** | 208.8s | 264.6s | **Hyrax** |
| **Verify Time** | 390.1ms | 289.5ms | **KZG** |
| **Proof Size** | 246.4 KB | 117.6 KB | **KZG** |

---

## Prover Time Breakdown

### Hyrax Mode

| Component | Time | % |
|-----------|------|---|
| R1CS Sat Proof (sumcheck) | 3.45s | 1.7% |
| Instance Evaluations | 0.36s | 0.2% |
| EqPolynomial evaluation | 0.10s | 0.0% |
| Derefs computation | 0.14s | 0.1% |
| **Derefs commitment (Hyrax MSM)** | **166.2s** | **79.6%** |
| Network construction | 4.07s | 1.9% |
| Network proof (Hyrax openings) | 34.5s | 16.5% |
| **TOTAL** | **208.8s** | 100% |

**Bottleneck**: Derefs commitment (Hyrax MSM) - 80% of proving time

### KZG Mode

| Component | Time | % |
|-----------|------|---|
| R1CS Sat Proof (sumcheck) | 3.55s | 1.3% |
| Instance Evaluations | 0.36s | 0.1% |
| EqPolynomial evaluation | 0.10s | 0.0% |
| Derefs computation | 0.20s | 0.1% |
| Derefs commitment (KZG) | 100.5s | 38.0% |
| Network construction | 5.03s | 1.9% |
| **Network proof (KZG openings)** | **154.8s** | **58.5%** |
| **TOTAL** | **264.6s** | 100% |

**Bottleneck**: Network proof (KZG quotient polynomials) - 59% of proving time

---

## Commitment Cost

The main commitment in Spartan is the **Derefs commitment** (looked-up eq-polynomial values).

| Parameter | Value |
|-----------|-------|
| Max NNZ per matrix | 3,151,183 |
| Padded NNZ | 2^22 = 4,194,304 |
| Derefs polynomials | 2 (row, col) |
| **Total field elements committed** | **8,388,608** (~8.4M) |

### Commitment Output Size

| PCS | Commitment Size | Formula |
|-----|-----------------|---------|
| **Hyrax** | ~130 KB | O(√n) = 2048 points × 2 polys × 32B |
| **KZG** | ~64 B | O(1) = 1 point × 2 polys × 32B |

---

## Proof Size Breakdown

### Top-Level Structure

| Component | Hyrax | KZG |
|-----------|-------|-----|
| R1CS Sat Proof | 45.9 KB | 45.9 KB |
| R1CS Eval Proof | 200.4 KB | 71.6 KB |
| inst_evals (3 scalars) | 96 B | 96 B |
| **TOTAL** | **246.4 KB** | **117.6 KB** |

### R1CS Sat Proof (Same for Both)

| Component | Size | Description |
|-----------|------|-------------|
| comm_vars (witness commitment) | ~2 KB | O(√n) points |
| sc_proof_phase1 (sumcheck) | ~20 KB | log(n) rounds × DotProductProof |
| claims_phase2 (4 points) | 128 B | 4 × GroupElement |
| pok_claims_phase2 | ~1 KB | KnowledgeProof + ProductProof |
| proof_eq_sc_phase1 | ~0.5 KB | EqualityProof |
| sc_proof_phase2 (sumcheck) | ~20 KB | log(n) rounds × DotProductProof |
| proof_eq_sc_phase2 | ~0.5 KB | EqualityProof |
| **Subtotal** | **~46 KB** | |

### R1CS Eval Proof (Lookup Argument) - THE KEY DIFFERENCE

| Component | Hyrax | KZG | Notes |
|-----------|-------|-----|-------|
| **DerefsCommitment** | **~130 KB** | **64 B** | Hyrax O(√n), KZG O(1) |
| ProductLayerProof | ~35 KB | ~35 KB | Scalars + batched circuit proofs |
| HashLayerProof.eval_* | ~3 KB | ~3 KB | Evaluation scalars |
| HashLayerProof.proof_ops | ~15 KB | ~15 KB | Hyrax log-sized opening |
| HashLayerProof.proof_mem | ~15 KB | ~15 KB | Hyrax log-sized opening |
| **HashLayerProof.proof_derefs** | **~2 KB** | **~64 B** | Hyrax log(n), KZG O(1) |
| **Subtotal** | **~200 KB** | **~72 KB** | |

### Why Hyrax DerefsCommitment is Large

```
Derefs = 2 polynomials (row_derefs, col_derefs)
Each polynomial has NNZ_padded = 2^22 = 4,194,304 evaluations

Hyrax commits to √n rows:
  √(4,194,304) = 2,048 points per polynomial
  2 polynomials × 2,048 points × 32 bytes = 131,072 bytes ≈ 128 KB

KZG commits to single point:
  2 polynomials × 1 point × 32 bytes = 64 bytes
```

---

## How to Run

```bash
# Hyrax mode (no trusted setup)
RAYON_NUM_THREADS=1 cargo run --release --example keyless_benchmark

# KZG mode (smaller proofs)
RAYON_NUM_THREADS=1 cargo run --release --features kzg --example keyless_benchmark
```

---

## Summary

| Use Case | Recommended Mode |
|----------|------------------|
| **Fastest proving** | Hyrax (209s vs 265s) |
| **Smallest proofs** | KZG (118 KB vs 246 KB) |
| **No trusted setup** | Hyrax |
| **Fastest verification** | KZG (290ms vs 390ms) |

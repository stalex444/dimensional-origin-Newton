# The Dimensional Origin of Newton's Constant

Companion repository for the Gravity Research Foundation essay (2026)

**Stephanie Alexander** · Double alumna, University of Chicago

## Overview

This repository contains the computational materials supporting:

> S. Alexander, "The Dimensional Origin of Newton's Constant," submitted to the Gravity Research Foundation Essay Competition (2026).

The essay shows that Newton's constant emerges from the polynomial family xⁿ = x + 1, evaluated at the dimensions of spacetime. The 3D root ρ = 1.32472 (x³ = x + 1) forces convergence; the 4D root Q = 1.22074 (x⁴ = x + 1) does not. A Lagrangian built from these roots with zero free parameters predicts gauge couplings, mass hierarchies, and mixing angles. A single correction factor — Q²/(2Q − 1), the quantitative completion of Ehrenfest's 1917 result — yields Newton's constant to 0.003%.

## Repository Contents

| File | Reproduces | Description |
|------|-----------|-------------|
| `table2_predictions.py` | Table 2 | Five predictions from the common ruler ρQ: α⁻¹, sin²θ_W, α_s, Y_p, α/α_G |
| `gravity_correction.py` | Table 4 | The gravity correction: tree-level → corrected (3.27% → 0.003%), Jacobson's η, the −23 discriminant bridge |
| `expanded_exclusion_test.py` | Figure 1 | Exclusion of 3,828 algebraic pairs against five simultaneous observables (α⁻¹, sin²θ_W, α_s, Y_p, m_τ/m_e), with parasite analysis |

Each script is self-contained, requires only NumPy, and runs in under one second.

## Quick Start

```
pip install numpy
python table2_predictions.py
python gravity_correction.py
python expanded_exclusion_test.py
```

Or run everything in the browser: [Open in Colab](notebooks/)

## Key Results

**Table 2 — Five predictions, zero free parameters:**

| Quantity | Formula | Predicted | Measured | Error |
|----------|---------|-----------|----------|-------|
| α⁻¹ | (ρQ)¹⁵/π² | 137.063 | 137.036 | 0.020% |
| sin²θ_W | λ₄/ψ³ | 0.2311 | 0.2312 | 0.057% |
| α_s(M_Z) | λ₄³/(4λ₃³ψ²) | 0.1182 | 0.1180 | 0.161% |
| Y_p | λ₃ = 1−1/ρ | 0.2451 | 0.2449 | 0.091% |
| α/α_G | (ρQ)²⁰⁹/π² | 4.31×10⁴² | 4.17×10⁴² | 3.27% |

**Table 4 — The gravity correction:**

| Quantity | Tree level | Corrected | Measured | Improvement |
|----------|-----------|-----------|----------|-------------|
| G (×10⁻¹¹) | 6.456 | 6.674 | 6.674 | 3.27% → 0.003% |
| M_P/m_e (×10²²) | 2.429 | 2.389 | 2.389 | 1.68% → 0.001% |

**Figure 1 — Exclusion test:** of 3,828 algebraic pairs tested against five simultaneous constraints (α⁻¹, sin²θ_W, α_s, Y_p, and m_τ/m_e), exactly one passes all five: the Pisot boundary pair (ρ, Q). Combined error: 0.41%. Every near-miss borrows ρ from PDT and pairs it with a numerical approximation to Q. No pair not containing ρ passes even 3/5.

## License

MIT
**Figure 1** — Exclusion test: of 3,828 algebraic pairs, exactly one passes all five constraints. Every near-miss borrows ρ from PDT and pairs it with a numerical approximation to Q. No pair not containing ρ falls below 25%.

## License

MIT

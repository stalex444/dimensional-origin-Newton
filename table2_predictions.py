"""
TABLE 2 REPRODUCER
==================
Reproduces Table 2 from:
  S. Alexander, "The Dimensional Origin of Newton's Constant" (2026)
  Gravity Research Foundation Essay.

Five predictions from zero free parameters.

The two roots:
  ρ = 1.32472 — real root of x³ = x + 1 (3D inflation constant)
  Q = 1.22074 — real root of x⁴ = x + 1 (4D inflation constant)

Derived quantities:
  λ₃ = 1 − 1/ρ   (3D self-coupling)
  λ₄ = 1 − 1/Q   (4D self-coupling)
  ψ  = Q/ρ        (dimensional ratio)
  w  = ρQ          (common ruler)

All experimental values: CODATA 2022 / PDG 2024 / Schramm & Turner (1998).
"""

import numpy as np

# ============================================================
# ROOTS — computed from the defining polynomials
# ============================================================

# x³ − x − 1 = 0
rho = max(r.real for r in np.roots([1, 0, -1, -1])
          if abs(r.imag) < 1e-10 and r.real > 1)

# x⁴ − x − 1 = 0
Q = max(r.real for r in np.roots([1, 0, 0, -1, -1])
        if abs(r.imag) < 1e-10 and r.real > 1)

# Derived quantities
w = rho * Q               # common ruler
lam3 = 1 - 1 / rho        # λ₃
lam4 = 1 - 1 / Q          # λ₄
psi = Q / rho              # ψ

print("=" * 72)
print("TABLE 2 REPRODUCER")
print("The Dimensional Origin of Newton's Constant (Alexander, 2026)")
print("=" * 72)

print(f"""
  Roots:
    ρ  = {rho:.6f}    (real root of x³ = x + 1)
    Q  = {Q:.6f}    (real root of x⁴ = x + 1)

  Derived:
    ρQ = {w:.5f}     (common ruler)
    λ₃ = 1 − 1/ρ = {lam3:.4f}
    λ₄ = 1 − 1/Q = {lam4:.4f}
    ψ  = Q/ρ     = {psi:.4f}
    ψ² = (Q/ρ)²  = {psi**2:.4f}
""")

# ============================================================
# FIVE PREDICTIONS
# ============================================================

predictions = []

# 1. Fine structure constant inverse
#    α⁻¹ = (ρQ)¹⁵ / π²
#    Exponent 15 = dim(SO(4,2)), the conformal group
pred_alpha = w ** 15 / np.pi ** 2
meas_alpha = 137.035999
err_alpha = abs(pred_alpha - meas_alpha) / meas_alpha * 100
predictions.append(("α⁻¹ (fine structure)", "(ρQ)¹⁵/π²",
                     pred_alpha, meas_alpha, err_alpha, "CODATA 2022"))

# 2. Weak mixing angle
#    sin²θ_W = λ₄ / ψ³
pred_sw = lam4 / psi ** 3
meas_sw = 0.23121
err_sw = abs(pred_sw - meas_sw) / meas_sw * 100
predictions.append(("sin²θ_W (weak mixing)", "λ₄/ψ³",
                     pred_sw, meas_sw, err_sw, "PDG 2024"))

# 3. Strong coupling at M_Z
#    α_s(M_Z) = λ₄³ / (4λ₃³ψ²)
pred_as = lam4 ** 3 / (4 * lam3 ** 3 * psi ** 2)
meas_as = 0.1180
err_as = abs(pred_as - meas_as) / meas_as * 100
predictions.append(("α_s(M_Z) (strong coupling)", "λ₄³/(4λ₃³ψ²)",
                     pred_as, meas_as, err_as, "PDG 2024"))

# 4. Primordial helium fraction
#    Y_p = λ₃ = 1 − 1/ρ
pred_yp = lam3
meas_yp = 0.2449
err_yp = abs(pred_yp - meas_yp) / meas_yp * 100
predictions.append(("Y_p (primordial helium)", "λ₃ = 1−1/ρ",
                     pred_yp, meas_yp, err_yp, "Schramm & Turner 1998"))

# 5. Gauge hierarchy (gravitational coupling ratio)
#    α/α_G = (ρQ)²⁰⁹ / π²
#    Exponent 209 = 15² − 15 − 1 (polynomial self-map at x = 15)
pred_hier = w ** 209 / np.pi ** 2
meas_hier = 4.17e42
err_hier = abs(pred_hier - meas_hier) / meas_hier * 100
predictions.append(("α/α_G (gauge hierarchy)", "(ρQ)²⁰⁹/π²",
                     pred_hier, meas_hier, err_hier, "CODATA 2022"))

# ============================================================
# DISPLAY
# ============================================================

print("  TABLE 2: Predictions from the common ruler ρQ")
print("  " + "-" * 68)
print(f"  {'Quantity':<30s} {'Formula':<16s} {'Predicted':>12s}"
      f" {'Measured':>12s} {'Error':>8s}")
print("  " + "-" * 68)

for name, formula, pred, meas, err, source in predictions:
    if pred > 1e6:
        pred_str = f"{pred:.2e}"
        meas_str = f"{meas:.2e}"
    else:
        # Match significant figures to the measured value
        if meas > 100:
            pred_str = f"{pred:.3f}"
            meas_str = f"{meas:.6f}"
        else:
            pred_str = f"{pred:.4f}"
            meas_str = f"{meas:.4f}"
    print(f"  {name:<30s} {formula:<16s} {pred_str:>12s}"
          f" {meas_str:>12s} {err:>7.3f}%")

print("  " + "-" * 68)

total_err = sum(p[4] for p in predictions)
mean_err = total_err / len(predictions)

print(f"""
  Combined error (sum):  {total_err:.3f}%
  Mean error:            {mean_err:.3f}%
  Free parameters:       0

  Sources:
    α⁻¹:     CODATA 2022 (NIST)
    sin²θ_W:  PDG 2024 (MS-bar, M_Z)
    α_s:      PDG 2024 (M_Z)
    Y_p:      Schramm & Turner, Rev. Mod. Phys. 70, 303 (1998)
    α/α_G:    Derived from CODATA 2022 (G, m_e, ħ, c)

  Note: The exponent 15 = dim(SO(4,2)), the conformal group of 3+1
  spacetime. The exponent 209 = 15² − 15 − 1, the polynomial self-map
  x³ = x + 1 evaluated at x = 15. Both are determined by group theory,
  not by fitting.
""")

"""
GRAVITY CORRECTION REPRODUCER
==============================
Reproduces Tables 4 and the Jacobson η computation from:
  S. Alexander, "The Dimensional Origin of Newton's Constant" (2026)
  Gravity Research Foundation Essay.

Key result: Forces propagating THROUGH spacetime see ρQ as a product.
Gravity IS spacetime and must resolve the product into its roots.

The correction factor Q²/(2Q − 1) emerges from the 4D self-coupling
acting on itself — the algebraic expression of Ehrenfest's degeneracy
(gravitational and centrifugal potentials share the same exponent in
4 spatial dimensions).

One correction, already in the published Lagrangian, improves
every gravitational prediction by a factor of ~1000.

All experimental values: CODATA 2022.
"""

import numpy as np

# ============================================================
# ROOTS
# ============================================================

# x³ − x − 1 = 0
rho = max(r.real for r in np.roots([1, 0, -1, -1])
          if abs(r.imag) < 1e-10 and r.real > 1)

# x⁴ − x − 1 = 0
Q = max(r.real for r in np.roots([1, 0, 0, -1, -1])
        if abs(r.imag) < 1e-10 and r.real > 1)

# Derived
w = rho * Q
lam3 = 1 - 1 / rho
lam4 = 1 - 1 / Q
psi = Q / rho

# ============================================================
# PHYSICAL CONSTANTS (CODATA 2022)
# ============================================================

G_measured = 6.67430e-11      # m³ kg⁻¹ s⁻²
m_e = 9.1093837015e-31        # kg (electron mass)
hbar = 1.054571817e-34        # J·s
c = 2.99792458e8              # m/s

# Derived measured quantities
# Planck mass: M_P = √(ħc/G)
M_P = np.sqrt(hbar * c / G_measured)
MP_me_measured = M_P / m_e
alpha_G_measured = (m_e / M_P) ** 2

print("=" * 72)
print("GRAVITY CORRECTION REPRODUCER")
print("The Dimensional Origin of Newton's Constant (Alexander, 2026)")
print("=" * 72)

# ============================================================
# THE CORRECTION FACTOR
# ============================================================

correction = (2 * Q - 1) / Q ** 2
lam4_sq = lam4 ** 2

print(f"""
  Roots:
    ρ  = {rho:.6f}
    Q  = {Q:.6f}
    ρQ = {w:.5f}

  The gravitational correction:
    λ₄² = (1 − 1/Q)²     = {lam4_sq:.5f}   ({lam4_sq * 100:.2f}%)
    (2Q−1)/Q²             = {correction:.5f}   ({correction * 100:.2f}%)
    1 − λ₄²               = {1 - lam4_sq:.5f}
    Identity: (2Q−1)/Q² ≡ 1 − λ₄²  ✓  (difference: {abs(correction - (1 - lam4_sq)):.1e})

  Physical meaning:
    λ₄² = {lam4_sq * 100:.2f}% of the 4D sector screens against itself
    (Ehrenfest's degeneracy: same-exponent cancellation in 4D)
    {correction * 100:.2f}% survives to contribute to gravitational binding
""")

# ============================================================
# TABLE 4: Tree level vs. corrected
# ============================================================

# The derivation:
#   α⁻¹ = (ρQ)¹⁵/π²           → α = π²/(ρQ)¹⁵
#   α/α_G = (ρQ)²⁰⁹/π²        → α_G = α × π²/(ρQ)²⁰⁹
#                                      = π⁴/(ρQ)²²⁴       [tree level]
#
# Gravity resolves ρ and Q separately. The correction:
#   α_G = π⁴ Q² / [(ρQ)²²⁴ (2Q − 1)]                     [corrected]
#
# Converting to G: α_G = G m_e²/(ħc), so G = α_G × ħc/m_e²

hbar_c = hbar * c  # ħc in SI

# Tree level: α_G = π⁴/(ρQ)²²⁴
alpha_G_tree = np.pi ** 4 / w ** 224

# Corrected: α_G = π⁴ Q² / [(ρQ)²²⁴ (2Q − 1)]
alpha_G_corr = np.pi ** 4 * Q ** 2 / (w ** 224 * (2 * Q - 1))

# Convert to G
G_tree = alpha_G_tree * hbar_c / m_e ** 2
G_corr = alpha_G_corr * hbar_c / m_e ** 2

# Convert to M_P/m_e (since α_G = (m_e/M_P)², M_P/m_e = 1/√α_G)
MP_me_tree = 1.0 / np.sqrt(alpha_G_tree)
MP_me_corr = 1.0 / np.sqrt(alpha_G_corr)

# Errors
G_tree_err = abs(G_tree - G_measured) / G_measured * 100
G_corr_err = abs(G_corr - G_measured) / G_measured * 100
MP_tree_err = abs(MP_me_tree - MP_me_measured) / MP_me_measured * 100
MP_corr_err = abs(MP_me_corr - MP_me_measured) / MP_me_measured * 100
aG_tree_err = abs(alpha_G_tree - alpha_G_measured) / alpha_G_measured * 100
aG_corr_err = abs(alpha_G_corr - alpha_G_measured) / alpha_G_measured * 100

print("  TABLE 4: Gravitational predictions")
print("  " + "=" * 66)
print(f"  {'Quantity':<20s} {'Tree level':>14s} {'Corrected':>14s}"
      f" {'Measured':>14s} {'Tree err':>9s} {'Corr err':>9s}")
print("  " + "-" * 66)

print(f"  {'G (×10⁻¹¹)':<20s} {G_tree * 1e11:14.3f} {G_corr * 1e11:14.3f}"
      f" {G_measured * 1e11:14.3f} {G_tree_err:8.2f}% {G_corr_err:8.3f}%")

print(f"  {'M_P/m_e (×10²²)':<20s} {MP_me_tree / 1e22:14.3f} {MP_me_corr / 1e22:14.3f}"
      f" {MP_me_measured / 1e22:14.3f} {MP_tree_err:8.2f}% {MP_corr_err:8.3f}%")

print(f"  {'α_G (×10⁻⁴⁵)':<20s} {alpha_G_tree * 1e45:14.3f} {alpha_G_corr * 1e45:14.3f}"
      f" {alpha_G_measured * 1e45:14.3f} {aG_tree_err:8.2f}% {aG_corr_err:8.3f}%")

print("  " + "-" * 66)
print(f"\n  Improvement factor: ~{G_tree_err / G_corr_err:.0f}×")

# Verify the correction is exactly Q²/(2Q-1)
ratio = alpha_G_corr / alpha_G_tree
expected_ratio = Q ** 2 / (2 * Q - 1)
print(f"\n  Correction ratio α_G(corr)/α_G(tree) = {ratio:.6f}")
print(f"  Q²/(2Q−1)                            = {expected_ratio:.6f}")
print(f"  Match: {'✓' if abs(ratio - expected_ratio) < 1e-10 else '✗'}")

# ============================================================
# JACOBSON'S η
# ============================================================

print(f"""
  {"=" * 66}
  JACOBSON'S ENTROPY DENSITY η
  {"=" * 66}

  Jacobson (1995) derived Einstein's equation as an equation of state,
  with Newton's constant set by an entropy density η that he left
  undetermined. The present framework provides:

    α_G = π⁴ Q² / [(ρQ)²²⁴ (2Q − 1)]

  Every factor is determined by the polynomial family and spacetime
  dimension:
    π⁴ = (π²)²    — squared projective volume (mass ratio → coupling)
    224 = 15 + 209 — electromagnetic + hierarchy staircase steps
    Q²/(2Q−1)      — gravitational screening (Ehrenfest correction)

  The resulting value of G:
    G_predicted = {G_corr:.5e} m³ kg⁻¹ s⁻²
    G_measured  = {G_measured:.5e} m³ kg⁻¹ s⁻²
    Error:        {G_corr_err:.3f}%
""")

# ============================================================
# THE −23 DISCRIMINANT BRIDGE
# ============================================================

print(f"  {'=' * 66}")
print(f"  THE −23 DISCRIMINANT BRIDGE")
print(f"  {'=' * 66}")
print(f"""
  Two independent computations in two different number fields:

  1. Norm of (2Q − 1) in ℚ(Q):
     The minimal polynomial of Q is f(x) = x⁴ − x − 1.
     N(2Q − 1) = 2⁴ × f(1/2) = 16 × [(1/2)⁴ − 1/2 − 1]
               = 16 × [1/16 − 3/2] = 16 × (−23/16) = −23""")

f_half = (0.5) ** 4 - 0.5 - 1
norm_2Q1 = 16 * f_half
print(f"\n     Computed: N(2Q − 1) = {norm_2Q1:.0f}")

# Discriminant of x³ − x − 1
disc = -4 * (-1) ** 3 - 27 * (-1) ** 2
print(f"""
  2. Discriminant of x³ − x − 1 (the polynomial defining ρ):
     Δ = −4(−1)³ − 27(−1)² = 4 − 27 = −23

     Computed: Δ = {disc}

  Both equal −23. The gravitational correction and the 3D polynomial
  share the same ramification prime. This is the algebraic bridge
  through which 3D and 4D communicate in gravitational physics.
""")

# ============================================================
# THE DEFICIT
# ============================================================

phi = (1 + np.sqrt(5)) / 2

print(f"  {'=' * 66}")
print(f"  THE DEFICIT: φ vs ρQ")
print(f"  {'=' * 66}")
print(f"""
  φ  = {phi:.5f}   (golden ratio — exact coherence in 2D)
  ρQ = {w:.5f}   (closest approach in 3+1 dimensions)

  Deficit: φ − ρQ = {phi - w:.5f}  ({(phi - w) / phi * 100:.3f}%)

  Unit-norm identity: N(ρQ) = −1, proving the relationship
  between ρQ and φ is algebraically exact, not approximate.

  Newton's constant is this deficit — the algebraic cost
  of living in three dimensions.
""")

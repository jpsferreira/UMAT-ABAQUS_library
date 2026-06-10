# UMAT verification suite (standalone, no ABAQUS)

Numerical-correctness tests for the modular UMAT. They concatenate `src/` into a
`umat.f90` (via `generate.py`, the single source of truth for PROPS/NSTATEV
layout), compile a small Fortran driver with `gfortran`, run it, and assert in
Python. No ABAQUS, no third-party deps (numpy not required).

```bash
.venv/bin/python tests/run_all.py        # full suite
.venv/bin/python tests/tangent_check.py  # T1 only (battery)
.venv/bin/python tests/tangent_check.py neo_hooke   # one model, per-state
```

## What each module does

| Module | What it verifies |
|--------|------------------|
| `tangent_check.py` (T1) | Analytic `DDSDDE` vs a central-difference tangent built from the **Kirchhoff** stress `τ = J·σ` under the standard Miehe symmetric `F`-perturbation. This reproduces ABAQUS's Jaumann material Jacobian `C_ijkl = (1/J)∂(Jσ_ij)/∂ε_kl`. |
| `additive_check.py` (T2) | Volumetric + isotropic + anisotropic in isolation, and `σ(all) == σ(iso)+σ(aniso)−σ(vol)` (the additive-decomposition invariant). Also: vol stress spherical, iso/aniso vanish at `F=I`. |
| `oracle_check.py` (T3) | Absolute closed-form oracles: volumetric `p(J)=K/2(J−1/J)`, Neo-Hooke isochoric uniaxial, Ogden(α=2) ≡ Neo-Hooke(C10=μ/2), Ogden isochoric-stress J-independence (neig==3), and fiber energy vanishing in compression. |
| `evolution_check.py` (T4) | Multi-step time evolution (carries STATEV): viscoelastic stress relaxes to the elastic equilibrium under hold; damage is monotone non-decreasing and frozen on unload. |
| `doc_consistency.py` (T5) | Pure-Python: parses the `ip = ip + N` advances in `umat_builder.f90` (ground truth) and asserts the `PROPS_REFERENCE.md` table counts and `generate.py` examples + NSTATEV all agree. No compile. |
| `network_check.py` (T6) | Filament-network battery: smoke (all 6 types finite/non-NaN, symmetric DDSDDE), RW quadrature loading (non-zero with `sphere_int<N>c.inp`, silent-zero without → E6), B3 prefdir sensitivity, and **gating** network tangent-consistency (all 6 types < 1e-4 after C5). |
| `recombination_check.py` (T7) | Combined-feature recombination: AI-fiber must survive viscosity (D1), and a pure network + viscosity must not be double-counted (D2). Stress oracles (FD can't see these — stress and tangent are consistently wrong). |
| `robustness_check.py` (T8) | Graceful degradation: an inverted element (det F < 0) returns `pnewdt < 1` (cutback) with finite stress, not NaN (E3); normal F leaves `pnewdt = 1`. |
| `contractile_tangent_check.py` (B6) | Contractile consistent tangent via **evolve-then-perturb**: evolve to active contraction, save STATEV, then compare analytic DDSDDE vs Kirchhoff-FD from that state. Catches the chi/raw-stretch/density bug the single-step checks can't (contractile gives zero stress at the first increment). |

`harness.py` holds the shared build/run helpers and deformation-gradient
generators.

## Calibration note (why the FD differences Kirchhoff, not Cauchy)

`DDSDDE_ijkl = ∂σ_ij/∂ε_kl + σ_ij·δ_kl`. Differencing Cauchy `σ` alone omits the
`σ_ij·δ_kl` term and overstates the error by ~1% of the stress. Differencing
`τ = J·σ` (and dividing by `J_base`) injects that term automatically, because the
normal perturbations change `J`. With this, all verified-correct models match to
~1e-10 (central differences, ε=1e-6). `det` is read from `statev(1)`.

## Known-bug markers (see `tasks/todo.md`)

Some checks are expected to FAIL until a fix lands; they print inline but do not
gate the suite:
- **B1** — `humphrey_hgo` (κ=0.226) tangent fails ~10% (dispersion coupling
  dropped). The κ=0 control passes at ~1e-10, pinning the bug to dispersion.
- **C1** — active-damage tangent fails ~8% (missing `dmg_diff·S⊗S` softening term;
  stress is correct, only the tangent is off, and only while damage grows).
- **B4** — fiber models are non-smooth at `I4=1`; a tangent check exactly on that
  kink (simple shear with fiber along the shear direction) is ill-posed, so
  fiber tangents are tested with the fiber engaged (I4>1).

When B1/C1 are fixed, flip their `tag` from `suspect` to `good` so the suite
gates on them.

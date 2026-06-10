# UMAT Builder — PROPS Layout Reference

## Overview

The UMAT builder uses a single `SUBROUTINE UMAT` that assembles material
behaviour from composable building blocks. The user selects which components
to combine through the ABAQUS `*USER MATERIAL` property array (PROPS).

## PROPS Header (positions 1–7)

| Position | Name          | Description                                          |
|----------|---------------|------------------------------------------------------|
| 1        | KBULK         | Bulk modulus (always required)                       |
| 2        | ISO_TYPE      | Isotropic model ID (see table below)                 |
| 3        | ANISO_TYPE    | Anisotropic model ID (see table below)               |
| 4        | N_FIBER_FAM   | Number of fiber families (0, 1, or 2)                |
| 5        | NETWORK_TYPE  | Network model ID (see table below)                   |
| 6        | DAMAGE_TYPE   | Damage model ID (see table below)                    |
| 7        | N_VISCO       | Number of Maxwell viscoelastic branches (0 = none)   |

## Model Type IDs

### Isotropic (ISO_TYPE)

| ID | Model          | Parameters                  | Count |
|----|----------------|-----------------------------|-------|
| 0  | None           | —                           | 0     |
| 1  | Neo-Hookean    | C10                         | 1     |
| 2  | Mooney-Rivlin  | C10, C01                    | 2     |
| 3  | Ogden          | N, mu_1, α_1, ..., mu_N, α_N | 1+2N |
| 4  | Humphrey (exp) | C10, C01                    | 2     |

### Anisotropic (ANISO_TYPE) — per fiber family

| ID | Model                       | Parameters per family                          | Count |
|----|-----------------------------|------------------------------------------------|-------|
| 0  | None                        | —                                              | 0     |
| 1  | HGO (dispersed)             | K1, K2, κ, fiber_x, fiber_y, fiber_z           | 6     |
| 2  | Humphrey fiber              | K1, K2, fiber_x, fiber_y, fiber_z              | 5     |
| 3  | HGO — angular integration   | K1, K2, b_disp, factor, fiber_x, fiber_y, fiber_z | 7  |
| 4  | Humphrey — angular integ.   | K1, K2, b_disp, factor, fiber_x, fiber_y, fiber_z | 7  |
| 5  | Humphrey + activation       | K1, K2, fiber_x, fiber_y, fiber_z, T0M         | 6     |

- IDs 3/4 integrate the fiber response over an icosahedral sphere: `factor` is the
  integer refinement level and `b_disp` the von Mises fiber-dispersion concentration.
- ID 5 adds an active (muscle) stress driven by the `PREDEF` field; `T0M` is the
  maximum active fiber stress.

### Network (NETWORK_TYPE)

All network models share an eight-element filament-property block
`filprops = [L, R0F, μ0, β, B0, λ0, R0C, ETAC]`:

- `L`   — filament contour length
- `R0F` — filament end-to-end distance
- `μ0`  — bending stiffness
- `β`   — power-law exponent
- `B0`  — thermal persistence parameter
- `λ0`  — reference length density
- `R0C` — crosslinker end-to-end distance
- `ETAC` — crosslinker fraction (linkers active only when `0 < ETAC < 1`; set
  `R0C = ETAC = 0` to disable). **Both values must always be present in PROPS**
  because the reader always consumes them.

`factor` = icosahedron refinement level (integer); factor=6 gives ~720 integration directions.
`pdir`  = preferred direction (reference config) for the orientation-density function.

**Quadrature-weight integration** (requires `sphere_intXXc.inp` loaded via UEXTERNALDB):

| ID | Model       | Parameters                                                  | Count |
|----|-------------|-------------------------------------------------------------|-------|
| 0  | None        | —                                                           | 0     |
| 1  | Affine      | PHI, N, B_orient, EFI, pdir_x, pdir_y, pdir_z, filprops(8)  | 15    |
| 2  | Non-affine  | PHI, N, B_orient, EFI, PP, pdir_x, pdir_y, pdir_z, filprops(8) | 16  |

**Angular integration** (self-contained, no external files needed):

| ID | Model          | Parameters                                                              | Count |
|----|----------------|-------------------------------------------------------------------------|-------|
| 3  | Mixed-AI       | PHI, N_naff, PP, N_aff, B_orient, EFI, factor, pdir_x, pdir_y, pdir_z, filprops(8)      | 18    |
| 4  | Contractile-AI | PHI, N, B_orient, EFI, FRIC, FFMAX, factor, pdir_x, pdir_y, pdir_z, filprops(8), KCH(7) | 25    |
| 5  | Affine-AI      | PHI, N, B_orient, EFI, factor, pdir_x, pdir_y, pdir_z, filprops(8)                      | 16    |
| 6  | Non-affine-AI  | PHI, N, PP, factor, filprops(8)                                                         | 12    |

`filprops(8)` denotes the eight filament properties listed above, packed in order.
Counts are the exact number of values consumed by `umat_builder.f90`
(`network_contribution`); they match `generate.py --list`.

### Damage (DAMAGE_TYPE)

| ID | Model   | Parameters        | Count |
|----|---------|-------------------|-------|
| 0  | None    | —                 | 0     |
| 1  | Sigmoid | β_damage, ψ_half | 2     |

### Viscosity — per branch

Each Maxwell branch requires: τ (relaxation time), θ (viscous fraction) = 2 params.

## Parameter Packing (PROPS 8 onwards)

Parameters are packed sequentially after the 7-element header:

```
PROPS(8:) = [iso_params] [aniso_params × N_FIBER_FAM] [network_params] [damage_params] [visco_params × N_VISCO]
```

## Examples

### Example 1: Neo-Hookean (simplest case)

```
*USER MATERIAL, CONSTANTS=8
** KBULK, ISO_TYPE, ANISO, NFIB, NET, DMG, NVISCO, C10
  1000.0,  1,  0,  0,  0,  0,  0,  10.0
```

NSTATEV = 1

### Example 2: Mooney-Rivlin

```
*USER MATERIAL, CONSTANTS=9
** KBULK, ISO, ANISO, NFIB, NET, DMG, NVISCO, C10, C01
  1000.0,  2,  0,  0,  0,  0,  0,  6.3,  0.012
```

### Example 3: HGO anisotropic (Neo-Hooke + 1 fiber family with dispersion)

```
*USER MATERIAL, CONSTANTS=14
** KBULK, ISO, ANISO, NFIB, NET, DMG, NVISCO
  1000.0,  1,  1,  1,  0,  0,  0,
** C10 (NH)
  10.0,
** K1, K2, kappa, fiber_x, fiber_y, fiber_z
  100.0,  10.0,  0.226,  1.0,  0.0,  0.0
```

### Example 4: Humphrey matrix + 2 HGO fiber families

```
*USER MATERIAL, CONSTANTS=21
** KBULK, ISO, ANISO, NFIB, NET, DMG, NVISCO
  500.0,  4,  1,  2,  0,  0,  0,
** C10, C01 (Humphrey iso)
  2.0,  1.5,
** Family 1: K1, K2, kappa, fiber_x, fiber_y, fiber_z
  50.0,  5.0,  0.1,  0.707,  0.707,  0.0,
** Family 2: K1, K2, kappa, fiber_x, fiber_y, fiber_z
  50.0,  5.0,  0.1,  0.707,  -0.707,  0.0
```

### Example 5: Neo-Hookean with sigmoid damage

```
*USER MATERIAL, CONSTANTS=10
** KBULK, ISO, ANISO, NFIB, NET, DMG, NVISCO
  1000.0,  1,  0,  0,  0,  1,  0,
** C10
  10.0,
** beta_damage, psi_half
  5.0,  50.0
```

NSTATEV = 3 (det + damage + max_sef)

### Example 6: Mooney-Rivlin with 2 Maxwell branches

```
*USER MATERIAL, CONSTANTS=13
** KBULK, ISO, ANISO, NFIB, NET, DMG, NVISCO
  1000.0,  2,  0,  0,  0,  0,  2,
** C10, C01
  6.3,  0.012,
** tau_1, theta_1, tau_2, theta_2
  0.1,  0.3,  1.0,  0.2
```

NSTATEV = 1 + 9×2 = 19

### Example 7: Full composite — Humphrey + HGO fiber + damage + viscosity

```
*USER MATERIAL, CONSTANTS=23
** KBULK, ISO, ANISO, NFIB, NET, DMG, NVISCO
  500.0,  4,  1,  1,  0,  1,  1,
** C10, C01 (Humphrey iso)
  2.0,  1.5,
** K1, K2, kappa, fiber_x, fiber_y, fiber_z (HGO fiber)
  100.0,  10.0,  0.226,  1.0,  0.0,  0.0,
** beta_damage, psi_half
  5.0,  50.0,
** tau_1, theta_1
  0.5,  0.25
```

NSTATEV = 3 + 9 = 12

### Example 8: Affine network with angular integration (no external files)

```
*USER MATERIAL, CONSTANTS=23
** KBULK, ISO, ANISO, NFIB, NET, DMG, NVISCO
  500.0,  0,  0,  0,  5,  0,  0,
** PHI, N, B_orient, EFI, factor, pdir_x, pdir_y, pdir_z
  0.5,  1.0e6,  2.0,  1.0,  6,  1.0,  0.0,  0.0,
** L, R0F, mu0, beta, B0, lambda0, R0C, ETAC
  1.0,  0.1,  0.01,  2.0,  0.1,  1.0,  0.0,  0.0
```

### Example 9: Non-affine network with angular integration

```
*USER MATERIAL, CONSTANTS=19
** KBULK, ISO, ANISO, NFIB, NET, DMG, NVISCO
  500.0,  0,  0,  0,  6,  0,  0,
** PHI, N, PP, factor
  0.5,  1.0e6,  2.0,  6,
** L, R0F, mu0, beta, B0, lambda0, R0C, ETAC
  1.0,  0.1,  0.01,  2.0,  0.1,  1.0,  0.0,  0.0
```

### Example 10: Humphrey matrix + HGO fiber (angular integration, ANISO_TYPE=3)

```
*USER MATERIAL, CONSTANTS=16
** KBULK, ISO, ANISO, NFIB, NET, DMG, NVISCO
  500.0,  4,  3,  1,  0,  0,  0,
** C10, C01 (Humphrey iso)
  2.0,  1.5,
** K1, K2, b_disp, factor, fiber_x, fiber_y, fiber_z
  100.0,  10.0,  5.0,  6,  1.0,  0.0,  0.0
```

### Example 11: Humphrey matrix + Humphrey fiber with activation (ANISO_TYPE=5)

```
*USER MATERIAL, CONSTANTS=15
** KBULK, ISO, ANISO, NFIB, NET, DMG, NVISCO
  500.0,  4,  5,  1,  0,  0,  0,
** C10, C01 (Humphrey iso)
  2.0,  1.5,
** K1, K2, fiber_x, fiber_y, fiber_z, T0M  (T0M = max active stress; ACT from PREDEF)
  100.0,  10.0,  1.0,  0.0,  0.0,  50.0
```

### Example 12: Mixed network — affine + non-affine (angular integration, NETWORK_TYPE=3)

```
*USER MATERIAL, CONSTANTS=25
** KBULK, ISO, ANISO, NFIB, NET, DMG, NVISCO
  500.0,  0,  0,  0,  3,  0,  0,
** PHI, N_naff, PP, N_aff, B_orient, EFI, factor, pdir_x, pdir_y, pdir_z
  0.5,  5.0e5,  2.0,  5.0e5,  2.0,  1.0,  6,  1.0,  0.0,  0.0,
** L, R0F, mu0, beta, B0, lambda0, R0C, ETAC
  1.0,  0.1,  0.01,  2.0,  0.1,  1.0,  0.0,  0.0
```

### Example 13: Contractile network (angular integration, NETWORK_TYPE=4)

```
*USER MATERIAL, CONSTANTS=32
** KBULK, ISO, ANISO, NFIB, NET, DMG, NVISCO
  500.0,  0,  0,  0,  4,  0,  0,
** PHI, N, B_orient, EFI, FRIC, FFMAX, factor, pdir_x, pdir_y, pdir_z
  0.5,  0.2,  2.0,  1.0,  11.0,  11.0,  6,  1.0,  0.0,  0.0,
** L, R0F, mu0, beta, B0, lambda0, R0C, ETAC
  0.988,  0.804,  38600.0,  0.438,  0.065,  1.007,  0.014,  0.667,
** KCH(1..7) chemical rate constants
  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0
```

NSTATEV = 1 + 4 + 20·factor² = 725 (factor = 6)

## State Variables (NSTATEV)

| Range            | Content                                                        |
|------------------|----------------------------------------------------------------|
| 1                | Jacobian determinant J                                         |
| 2                | Damage variable d (if damage active)                           |
| 3                | Max historical SEF (if damage active)                          |
| ... : +9×V       | Hidden stress tensors (V = N_VISCO branches)                   |
| next 4           | Contractile chemical fractions FRAC(4) (if NETWORK_TYPE == 4)  |
| next nwp         | Contractile sliding displacements RU0 (if NETWORK_TYPE == 4)   |

Formula: `NSTATEV = 1 + (2 if damage) + (9 × N_VISCO) + (4 + nwp if NETWORK_TYPE == 4)`,
where `nwp = 20 × factor²` (icosahedron integration points; factor=6 → 720).

The blocks are laid out in this order: det, damage(2), viscous(9×V), contractile(4+nwp).

## Module Architecture

```
mod_constants       ← Named constants and precision
  ↓
mod_tensor          ← Tensor algebra (contractions, push/pull, eigenvalues)
  ↓
mod_kinematics      ← Deformation measures (F, C, B, invariants, stretch)
  ↓
mod_continuum       ← Stress/stiffness framework (vol/iso split, Voigt)
  ↓
mod_hyperelastic    ← Isotropic SEFs (NH, MR, Ogden, Humphrey)
mod_anisotropic     ← Fiber models (HGO, Humphrey-fiber)
mod_icosahedron     ← Icosahedron geometry for angular integration
mod_network         ← Filament networks (affine, non-affine, AI variants)
mod_damage          ← Damage evolution (sigmoid)
mod_viscosity       ← Viscoelasticity (generalized Maxwell)
  ↓
umat_builder        ← Single UMAT: reads PROPS config, dispatches to modules
```

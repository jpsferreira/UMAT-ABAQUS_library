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

| ID | Model            | Parameters per family                  | Count |
|----|------------------|----------------------------------------|-------|
| 0  | None             | —                                      | 0     |
| 1  | HGO (dispersed)  | K1, K2, κ, fiber_x, fiber_y, fiber_z  | 6     |
| 2  | Humphrey fiber   | K1, K2, fiber_x, fiber_y, fiber_z     | 5     |

### Network (NETWORK_TYPE)

**Quadrature-weight integration** (requires `sphere_intXXc.inp` loaded via UEXTERNALDB):

| ID | Model       | Parameters                                               | Count |
|----|-------------|----------------------------------------------------------|-------|
| 0  | None        | —                                                        | 0     |
| 1  | Affine      | PHI, N, B_orient, EFI, pdir_x, pdir_y, pdir_z, L, R0, μ0, β, B0, λ0 | 13 |
| 2  | Non-affine  | PHI, N, B_orient, EFI, PP, L, R0, μ0, β, B0, λ0        | 11    |

**Angular integration** (self-contained, no external files needed):

| ID | Model            | Parameters                                                            | Count |
|----|------------------|-----------------------------------------------------------------------|-------|
| 5  | Affine-AI        | PHI, N, B_orient, EFI, factor, pdir_x, pdir_y, pdir_z, L, R0, μ0, β, B0, λ0 | 14 |
| 6  | Non-affine-AI    | PHI, N, PP, factor, L, R0, μ0, β, B0, λ0                            | 10    |

- `factor` = icosahedron refinement level (integer). factor=6 gives ~720 integration points.
- `pdir` = preferred direction (reference config) for orientation density in affine models.

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
*USER MATERIAL, CONSTANTS=21
** KBULK, ISO, ANISO, NFIB, NET, DMG, NVISCO
  500.0,  0,  0,  0,  5,  0,  0,
** PHI, N, B_orient, EFI, factor, pdir_x, pdir_y, pdir_z
  0.5,  1.0e6,  2.0,  1.0,  6,  1.0,  0.0,  0.0,
** L, R0, mu0, beta, B0, lambda0
  1.0,  0.1,  0.01,  2.0,  0.1,  1.0
```

### Example 9: Non-affine network with angular integration

```
*USER MATERIAL, CONSTANTS=17
** KBULK, ISO, ANISO, NFIB, NET, DMG, NVISCO
  500.0,  0,  0,  0,  6,  0,  0,
** PHI, N, PP, factor
  0.5,  1.0e6,  2.0,  6,
** L, R0, mu0, beta, B0, lambda0
  1.0,  0.1,  0.01,  2.0,  0.1,  1.0
```

## State Variables (NSTATEV)

| Range          | Content                                      |
|----------------|----------------------------------------------|
| 1              | Jacobian determinant J                       |
| 2              | Damage variable d (if damage active)         |
| 3              | Max historical SEF (if damage active)        |
| 4 : 3+9×V     | Hidden stress tensors (V = number of branches)|

Formula: `NSTATEV = 1 + (2 if damage) + (9 × N_VISCO)`

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

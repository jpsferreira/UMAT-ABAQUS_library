# UMAT Generator Configuration Guide

## Quick Start

```bash
# Install
uv sync

# Generate from built-in example
uv run umat-generate --example neo_hooke
cd neo_hooke && make run

# Generate from custom JSON
uv run umat-generate my_material.json
cd my_material && make run

# List all model types and examples
uv run umat-generate --list
```

## JSON Config Schema

```json
{
  "name": "my_material",
  "kbulk": 1000.0,
  "iso_type": 1,
  "iso_params": [10.0],
  "aniso_type": 0,
  "n_fiber_fam": 0,
  "aniso_params": [],
  "network_type": 0,
  "network_params": [],
  "damage_type": 0,
  "damage_params": [],
  "n_visco": 0,
  "visco_params": [],
  "test": {
    "stretch_max": 1.5,
    "gamma_max": 0.6,
    "nsteps": 400,
    "dtime": 0.01
  }
}
```

### Required Fields

| Field | Type | Description |
|-------|------|-------------|
| `name` | string | Output directory name |
| `kbulk` | float | Bulk modulus (volumetric penalty, typically 100-10000) |
| `iso_type` | int | Isotropic model type (see below) |
| `iso_params` | list | Parameters for the isotropic model |
| `aniso_type` | int | Anisotropic model type (see below) |
| `n_fiber_fam` | int | Number of fiber families (0, 1, or 2) |
| `aniso_params` | list | Parameters per fiber family |
| `network_type` | int | Filament network model type |
| `network_params` | list | Network model parameters |
| `damage_type` | int | Damage model type |
| `damage_params` | list | Damage parameters |
| `n_visco` | int | Number of Maxwell branches (0-3) |
| `visco_params` | list | Viscous params: [tau1, theta1, tau2, theta2, ...] |

### Optional Fields

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `test.stretch_max` | float | 1.5 | Max stretch for uniaxial/biaxial tests |
| `test.gamma_max` | float | 0.6 | Max shear strain |
| `test.nsteps` | int | 400 | Number of loading steps |
| `test.dtime` | float | 0.01 | Time increment (important for visco models) |

## Model Types

### Isotropic (`iso_type`)

| ID | Model | Parameters | Typical Values |
|----|-------|-----------|----------------|
| 0 | None | — | — |
| 1 | Neo-Hookean | `[C10]` | C10 = 1-100 kPa |
| 2 | Mooney-Rivlin | `[C10, C01]` | C10 = 1-100, C01 = 0.01-10 |
| 3 | Ogden (N-term) | `[N, mu1, alpha1, ..., muN, alphaN]` | N=1-3, mu~1, alpha~2-5 |
| 4 | Humphrey exponential | `[C10, C01]` | C10 = 0.1-10, C01 = 0.1-10 |

### Anisotropic (`aniso_type`)

Each fiber family gets the same parameter set. Set `n_fiber_fam` to the number of families.

| ID | Model | Parameters per family | Notes |
|----|-------|----------------------|-------|
| 0 | None | — | — |
| 1 | HGO (dispersed) | `[K1, K2, kappa, fx, fy, fz]` | kappa in [0, 1/3], (fx,fy,fz) = unit vector |
| 2 | Humphrey fiber | `[K1, K2, fx, fy, fz]` | (fx,fy,fz) = fiber direction unit vector |
| 3 | HGO (AI discrete) | `[K1, K2, bdisp, factor, fx, fy, fz]` | Icosahedron integration, bdisp = dispersion |
| 4 | Humphrey (AI discrete) | `[K1, K2, bdisp, factor, fx, fy, fz]` | Icosahedron integration, stretch-based |
| 5 | Humphrey + activation | `[K1, K2, fx, fy, fz, T0M]` | T0M = max active stress, ACT from PREDEF |

### Network (`network_type`)

| ID | Model | Parameters |
|----|-------|-----------|
| 0 | None | — |
| 1 | Affine (RW) | `[PHI, N, B_orient, EFI, pdir(3), L, R0F, mu0, beta, B0, lambda0, R0C, ETAC]` |
| 2 | Non-affine (RW) | `[PHI, N, B_orient, EFI, PP, L, R0F, mu0, beta, B0, lambda0, R0C, ETAC]` |
| 3 | Mixed (AI) | `[PHI, N_naff, PP, N_aff, B_orient, EFI, factor, pdir(3), L, R0F, mu0, beta, B0, lambda0, R0C, ETAC]` |
| 4 | Contractile (AI) | `[PHI, N, B_orient, EFI, FRIC, FFMAX, factor, pdir(3), L, R0F, mu0, beta, B0, lambda0, R0C, ETAC, KCH(7)]` |
| 5 | Affine (AI) | `[PHI, N, B_orient, EFI, factor, pdir(3), L, R0F, mu0, beta, B0, lambda0, R0C, ETAC]` |
| 6 | Non-affine (AI) | `[PHI, N, PP, factor, L, R0F, mu0, beta, B0, lambda0, R0C, ETAC]` |

- **RW** = Random Walk (requires sphere quadrature files at runtime)
- **AI** = Angular Integration (icosahedron, self-contained)
- **factor** = icosahedron refinement level (1-10, higher = more integration points)
- **PP** = non-affinity exponent

#### Network Parameter Meanings

| Param | Meaning | Typical Range |
|-------|---------|--------------|
| PHI | Volume fraction of network | 0-1 |
| N | Number of filaments per chain | 1e4-1e8 |
| B_orient | Orientation bias | 0-5 |
| EFI | End-to-end factor | 0.5-2.0 |
| L | Filament contour length | 0.1-10 |
| R0 | End-to-end distance | 0.01-1 |
| mu0 | Bending stiffness | 0.001-1 |
| beta | Stretching exponent | 1-5 |
| B0 | Initial stiffness | 0.01-1 |
| lambda0 | Reference stretch | 0.5-2 |

### Damage (`damage_type`)

| ID | Model | Parameters |
|----|-------|-----------|
| 0 | None | — |
| 1 | Sigmoid | `[beta_d, psi_half]` — beta_d controls steepness, psi_half is strain energy at 50% damage |

### Viscosity (`n_visco`)

Set `n_visco` = 1, 2, or 3 for Maxwell branches. Parameters: `[tau1, theta1, tau2, theta2, ...]`

| Param | Meaning |
|-------|---------|
| tau | Relaxation time (seconds) |
| theta | Viscous weight fraction (0-1) |

## Combining Features

Features are composable. Any combination of iso + aniso + network + damage + visco is valid:

```json
{
  "name": "humphrey_hgo_visco",
  "kbulk": 500.0,
  "iso_type": 4, "iso_params": [2.0, 1.5],
  "aniso_type": 1, "n_fiber_fam": 1,
  "aniso_params": [100.0, 10.0, 0.226, 1.0, 0.0, 0.0],
  "network_type": 0, "network_params": [],
  "damage_type": 0, "damage_params": [],
  "n_visco": 1, "visco_params": [0.5, 0.25]
}
```

## Built-in Examples

Run `uv run umat-generate --list` for the full list. Current examples:

| Example | Description |
|---------|-------------|
| `neo_hooke` | Simplest isotropic hyperelastic |
| `mooney_rivlin` | Two-parameter isotropic |
| `ogden_3term` | Three-term Ogden |
| `humphrey_hgo` | Humphrey matrix + HGO fiber family |
| `humphrey_fiber` | Humphrey matrix + Humphrey fiber family |
| `neo_hooke_damage` | Neo-Hookean + sigmoid damage |
| `humphrey_hgo_damage` | HGO model + damage |
| `humphrey_fiber_damage` | Humphrey fiber + damage |
| `neo_hooke_visco` | Neo-Hookean + 1 Maxwell branch |
| `mooney_rivlin_visco` | Mooney-Rivlin + viscosity |
| `ogden_visco` | Ogden + viscosity |
| `humphrey_hgo_visco` | HGO + viscosity |
| `humphrey_fiber_visco` | Humphrey fiber + viscosity |
| `affine_network` | Affine filament network (AI) |
| `nonaffine_network` | Non-affine filament network (AI) |

## Generated Output

Each generated directory contains:

```
my_material/
  umat.f90          # Self-contained UMAT source
  test_umat.f90     # Standalone test (uniaxial, biaxial, shear, simple shear)
  Makefile          # make run → compile and test
  aba_param.inc     # ABAQUS include stub
  config.json       # Config snapshot for regeneration
  abaqus/
    cube.inp        # C3D8 single-element mesh
    sec.inp         # *User Material card with PROPS
    bcs_uni.inp     # Uniaxial BCs
    bcs_bi.inp      # Biaxial BCs
    bcs_sh.inp      # Shear BCs
    run.sh          # abaqus job=cube user=../umat.f90
```

## Running in ABAQUS

```bash
cd my_material/abaqus
./run.sh                    # Uniaxial (default)
./run.sh bcs_bi.inp         # Biaxial
```

## State Variables

`NSTATEV = 1 + (2 if damage) + (9 * n_visco) + (4 + nwp if contractile)`

Where `nwp = 20 * factor²` for contractile networks (icosahedron integration points).

| Index | Content | Condition |
|-------|---------|-----------|
| 1 | det(F) | Always |
| 2 | Damage variable d | If damage_type > 0 |
| 3 | Max strain energy | If damage_type > 0 |
| 4+ | Hidden stress tensors (3x3 per branch) | If n_visco > 0 |
| next 4 | Chemical fractions (4 species) | If network_type == 4 |
| next nwp | Sliding displacements per direction | If network_type == 4 |

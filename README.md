# UMAT-ABAQUS Library

A composable UMAT framework for finite-strain constitutive modelling of soft biological tissues and filamentous networks in ABAQUS.

## What this repository does

This library provides a **single UMAT subroutine** that assembles material behaviour at runtime from modular building blocks. Instead of maintaining dozens of nearly-identical Fortran files (one per material law), users combine isotropic, anisotropic, network, damage, and viscous contributions through the ABAQUS `*USER MATERIAL` property array (PROPS).

### Available building blocks

| Category | Models |
|----------|--------|
| **Isotropic** | Neo-Hookean, Mooney-Rivlin, Ogden (N-term), Humphrey exponential |
| **Anisotropic** | HGO (analytical dispersion or AI discrete), Humphrey fiber, Humphrey fiber with muscle activation |
| **Network** | Affine, Non-affine, Mixed, Contractile (RW or AI integration), with optional crosslinker coupling |
| **Damage** | Sigmoid evolution |
| **Viscosity** | Generalized Maxwell (up to 3 branches) |

Any combination of these can be activated in a single analysis by setting the appropriate flags in the PROPS header. See `src/PROPS_REFERENCE.md` for the full parameter layout and ABAQUS input file examples.

### Legacy material laws

The original single-file material laws are preserved in:
- `soft_tissues/` — hyperelastic and viscoelastic soft tissue models
- `biofilaments/` — filament network models (affine, non-affine, mixed, contractile)

These remain functional but are superseded by the modular `src/` implementation.

## Prerequisites

- **Python 3.8+** with [uv](https://docs.astral.sh/uv/) (for the generator)
- **Fortran compiler** — gfortran 7+ or Intel Fortran 2017+
- **GNU Make**
- **ABAQUS** 2017+ (only for ABAQUS testing, not for standalone)

```bash
# Set up Python environment
uv sync

# Install gfortran on Debian/Ubuntu
sudo apt install gcc gfortran make
```

See [docs/config_guide.md](docs/config_guide.md) for the full JSON config schema and parameter reference.

## Repository structure

```
generate.py                 Material law generator (main entry point)
validate.py                 Legacy vs modular comparison test suite
docs/config_guide.md        Full JSON config schema and parameter reference
src/                        Modular UMAT source modules
  mod_constants.f90           Precision parameters and model type IDs
  mod_tensor.f90              Tensor algebra (contractions, push/pull, spectral)
  mod_kinematics.f90          Deformation measures (F, C, B, invariants, stretches)
  mod_continuum.f90           Stress/stiffness framework (vol/iso split, Voigt, Jaumann)
  mod_hyperelastic.f90        Isotropic strain energy functions
  mod_anisotropic.f90         Fiber-reinforced models (HGO, Humphrey, AI discrete, activation)
  mod_icosahedron.f90         Icosahedron geometry for angular integration
  mod_network.f90             Filament network assembly (affine, non-affine, mixed, contractile)
  mod_damage.f90              Damage evolution
  mod_viscosity.f90           Generalized Maxwell viscoelasticity
  umat_builder.f90            Single UMAT entry point + network dispatch
  uexternaldb.f90             UEXTERNALDB for loading quadrature data (RW networks)
  aba_param.inc               ABAQUS implicit typing include
  PROPS_REFERENCE.md          Full PROPS layout documentation with examples

soft_tissues/               Legacy single-file soft tissue UMATs
biofilaments/               Legacy single-file filament network UMATs
```

## Generating a material law

The `generate.py` script creates a self-contained directory with everything needed to test a material law both standalone and in ABAQUS.

### From a built-in example

```bash
uv run umat-generate --example neo_hooke
```

Run `uv run umat-generate --list` for the full list. Current examples:

| Category | Examples |
|----------|----------|
| Isotropic | `neo_hooke`, `mooney_rivlin`, `ogden_3term` |
| Anisotropic | `humphrey_hgo`, `humphrey_fiber`, `humphrey_hgo_ai`, `humphrey_muscle` |
| Damage | `neo_hooke_damage`, `humphrey_hgo_damage`, `humphrey_fiber_damage` |
| Viscosity | `neo_hooke_visco`, `mooney_rivlin_visco`, `ogden_visco`, `humphrey_hgo_visco`, `humphrey_fiber_visco` |
| Networks | `affine_network`, `nonaffine_network`, `affine_network_linkers`, `mixed_network`, `contractile_network` |

### From a JSON configuration

```bash
uv run umat-generate my_material.json
```

Example `my_material.json`:

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

### List available model types and parameters

```bash
uv run umat-generate --list
```

### What gets generated

```
my_material/
  umat.f90            Concatenated UMAT source (for ABAQUS user= submission)
  test_umat.f90       Standalone test driver (uniaxial, biaxial, shear, simple shear)
  Makefile            Build & run the standalone test
  aba_param.inc       ABAQUS include file
  config.json         Saved configuration (for regeneration)
  abaqus/
    cube.inp          Single C3D8 element mesh
    sec.inp           *User Material card with PROPS values
    bcs_uni.inp       Uniaxial boundary conditions
    bcs_bi.inp        Biaxial boundary conditions
    bcs_sh.inp        Shear boundary conditions
    run.sh            ABAQUS submission script
```

### Running the standalone test

```bash
cd my_material
make run          # Compiles umat.f90 + test_umat.f90, runs all load cases
```

Output files are written to `my_material/results/` (uniaxial.dat, biaxial.dat, shear.dat, simple_shear.dat).

### Running in ABAQUS

```bash
cd my_material/abaqus
./run.sh                    # Uniaxial (default)
./run.sh bcs_bi.inp         # Biaxial
./run.sh bcs_sh.inp         # Shear
```

See `src/PROPS_REFERENCE.md` for the full parameter reference.

### Network models — integration scheme choice

Network models offer two integration approaches:

| NETWORK_TYPE | Model | Scheme | External files? |
|:---:|--------|--------|:---:|
| 1 | Affine | Quadrature weights (RW) | Yes (`sphere_intXXc.inp`) |
| 2 | Non-affine | Quadrature weights (RW) | Yes |
| 3 | Mixed (affine + non-affine) | Icosahedron (AI) | No |
| 4 | Contractile (active contraction) | Icosahedron (AI) | No |
| 5 | Affine | Icosahedron (AI) | No |
| 6 | Non-affine | Icosahedron (AI) | No |

The RW approach loads quadrature directions from `sphere_intXXc.inp` at analysis start via the `UEXTERNALDB` subroutine. The AI approach generates integration directions dynamically from an icosahedron subdivision, controlled by the `factor` parameter in PROPS.

All network types support optional crosslinker coupling via R0C and ETAC parameters (set ETAC=0 to disable).

## State variables (DEPVAR)

| Slot | Content | Condition |
|------|---------|-----------|
| 1 | Jacobian determinant J | Always |
| 2 | Damage variable d | If damage_type > 0 |
| 3 | Max strain energy | If damage_type > 0 |
| next 9 per branch | Hidden stress tensors (3x3) | If n_visco > 0 |
| next 4 | Chemical fractions (4 species) | If network_type == 4 |
| next nwp | Sliding displacements per direction | If network_type == 4 |

Formula: `NSTATEV = 1 + (2 if damage) + (9 * N_VISCO) + (4 + nwp if contractile)`

Where `nwp = 20 * factor^2` for contractile networks.

## License

GNU GPL 3.0

## References

J P S Ferreira et al. "On the mechanical response of the actomyosin cortex during cell indentations". In: *Biomechanics and Modeling in Mechanobiology* (2020). doi: [10.1007/s10237-020-01324-5](https://link.springer.com/article/10.1007%2Fs10237-020-01324-5)

J P S Ferreira et al. "Altered mechanics of vaginal smooth muscle cells due to the lysyl oxidase-like1 knockout". In: *Acta Biomaterialia* (2020). doi: [10.1016/j.actbio.2020.03.046](https://doi.org/10.1016/j.actbio.2020.03.046)

J A Lopez-Campos et al. "Characterization of hyperelastic and damage behavior of tendons". In: *Computer Methods in Biomechanics and Biomedical Engineering* (2020). doi: [10.1080/10255842.2019.1710742](https://doi.org/10.1080/10255842.2019.1710742)

M C P Vila Pouca et al. "Viscous effects in pelvic floor muscles during childbirth: A numerical study". In: *International Journal for Numerical Methods in Biomedical Engineering* (2018). doi: [10.1002/cnm.2927](https://doi.org/10.1002/cnm.2927)

J P S Ferreira, M P L Parente, and R M Natal Jorge. "Continuum mechanical model for cross-linked actin networks with contractile bundles". In: *Journal of the Mechanics and Physics of Solids* (2017). doi: [10.1016/j.jmps.2017.09.009](https://doi.org/10.1016/j.jmps.2017.09.009)

J P S Ferreira et al. "A general framework for the numerical implementation of anisotropic hyperelastic material models including non-local damage". In: *Biomechanics and Modeling in Mechanobiology* (2017). doi: [10.1007/s10237-017-0875-9](https://link.springer.com/article/10.1007/s10237-017-0875-9)

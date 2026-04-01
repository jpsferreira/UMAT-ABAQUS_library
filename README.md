# UMAT-ABAQUS Library

A composable UMAT framework for finite-strain constitutive modelling of soft biological tissues and filamentous networks in ABAQUS.

## What this repository does

This library provides a **single UMAT subroutine** that assembles material behaviour at runtime from modular building blocks. Instead of maintaining dozens of nearly-identical Fortran files (one per material law), users combine isotropic, anisotropic, network, damage, and viscous contributions through the ABAQUS `*USER MATERIAL` property array (PROPS).

### Available building blocks

| Category | Models |
|----------|--------|
| **Isotropic** | Neo-Hookean, Mooney-Rivlin, Ogden (N-term), Humphrey exponential |
| **Anisotropic** | HGO with dispersion, Humphrey fiber (per fiber family, up to 2 families) |
| **Network** | Affine, Non-affine (quadrature weights or icosahedron angular integration) |
| **Damage** | Sigmoid evolution |
| **Viscosity** | Generalized Maxwell (up to 3 branches) |

Any combination of these can be activated in a single analysis by setting the appropriate flags in the PROPS header. See `src/PROPS_REFERENCE.md` for the full parameter layout and ABAQUS input file examples.

### Legacy material laws

The original single-file material laws are preserved in:
- `soft_tissues/` — hyperelastic and viscoelastic soft tissue models
- `biofilaments/` — filament network models (affine, non-affine, mixed, contractile)

These remain functional but are superseded by the modular `src/` implementation.

## Prerequisites

- **ABAQUS** 2017 or later (for UMAT/UEXTERNALDB support)
- **Fortran compiler** — one of:
  - Intel Fortran (ifort/ifx) 2017+
  - GNU Fortran (gfortran) 7+
- **GNU Make** (for building outside ABAQUS)

```bash
# Install gfortran on Debian/Ubuntu
sudo apt install gcc gfortran make
```

## Repository structure

```
generate.py                 Material law generator (main entry point)
src/                        Modular UMAT source modules
  mod_constants.f90           Precision parameters and model type IDs
  mod_tensor.f90              Tensor algebra (contractions, push/pull, spectral)
  mod_kinematics.f90          Deformation measures (F, C, B, invariants, stretches)
  mod_continuum.f90           Stress/stiffness framework (vol/iso split, Voigt, Jaumann)
  mod_hyperelastic.f90        Isotropic strain energy functions
  mod_anisotropic.f90         Fiber-reinforced models (HGO, Humphrey fiber)
  mod_icosahedron.f90         Icosahedron geometry for angular integration
  mod_network.f90             Filament network assembly (RW and AI variants)
  mod_damage.f90              Damage evolution
  mod_viscosity.f90           Generalized Maxwell viscoelasticity
  umat_builder.f90            Single UMAT entry point + network dispatch
  uexternaldb.f90             UEXTERNALDB for loading quadrature data (RW networks)
  aba_param.inc               ABAQUS implicit typing include
  PROPS_REFERENCE.md          Full PROPS layout documentation with examples
  Makefile                    Build rules for module compilation

soft_tissues/               Legacy single-file soft tissue UMATs
biofilaments/               Legacy single-file filament network UMATs
```

## Generating a material law

The `generate.py` script creates a self-contained directory with everything needed to test a material law both standalone and in ABAQUS.

### From a built-in example

```bash
python generate.py --example neo_hooke
```

Available examples: `neo_hooke`, `mooney_rivlin`, `humphrey_hgo`, `ogden_3term`, `neo_hooke_damage`, `neo_hooke_visco`, `affine_network`.

### From a JSON configuration

```bash
python generate.py my_material.json
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
python generate.py --list
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

| NETWORK_TYPE | Scheme | Requires external files? |
|:---:|--------|:---:|
| 1, 2 | Pre-loaded quadrature weights (RW) | Yes (`sphere_intXXc.inp`) |
| 5, 6 | Icosahedron angular integration (AI) | No |

The RW approach loads quadrature directions from `sphere_intXXc.inp` at analysis start via the `UEXTERNALDB` subroutine. The AI approach generates integration directions dynamically from an icosahedron subdivision, controlled by the `factor` parameter in PROPS.

## State variables (DEPVAR)

| Slot | Content |
|------|---------|
| 1 | Jacobian determinant J |
| 2 | Damage variable d (if damage active) |
| 3 | Max historical strain energy (if damage active) |
| 4+ | Hidden stress tensors for Maxwell branches (9 per branch) |

Formula: `NSTATEV = 1 + (2 if damage) + (9 * N_VISCO)`

## License

GNU GPL 3.0

## References

J P S Ferreira et al. "On the mechanical response of the actomyosin cortex during cell indentations". In: *Biomechanics and Modeling in Mechanobiology* (2020). doi: [10.1007/s10237-020-01324-5](https://link.springer.com/article/10.1007%2Fs10237-020-01324-5)

J P S Ferreira et al. "Altered mechanics of vaginal smooth muscle cells due to the lysyl oxidase-like1 knockout". In: *Acta Biomaterialia* (2020). doi: [10.1016/j.actbio.2020.03.046](https://doi.org/10.1016/j.actbio.2020.03.046)

J A Lopez-Campos et al. "Characterization of hyperelastic and damage behavior of tendons". In: *Computer Methods in Biomechanics and Biomedical Engineering* (2020). doi: [10.1080/10255842.2019.1710742](https://doi.org/10.1080/10255842.2019.1710742)

M C P Vila Pouca et al. "Viscous effects in pelvic floor muscles during childbirth: A numerical study". In: *International Journal for Numerical Methods in Biomedical Engineering* (2018). doi: [10.1002/cnm.2927](https://doi.org/10.1002/cnm.2927)

J P S Ferreira, M P L Parente, and R M Natal Jorge. "Continuum mechanical model for cross-linked actin networks with contractile bundles". In: *Journal of the Mechanics and Physics of Solids* (2017). doi: [10.1016/j.jmps.2017.09.009](https://doi.org/10.1016/j.jmps.2017.09.009)

J P S Ferreira et al. "A general framework for the numerical implementation of anisotropic hyperelastic material models including non-local damage". In: *Biomechanics and Modeling in Mechanobiology* (2017). doi: [10.1007/s10237-017-0875-9](https://link.springer.com/article/10.1007/s10237-017-0875-9)

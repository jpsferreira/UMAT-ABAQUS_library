#!/usr/bin/env python3
"""Generate a self-contained material law directory from a JSON configuration.

Usage:
    python generate.py config.json              # Generate from config
    python generate.py --example neo_hooke      # Generate example config + material
    python generate.py --list                   # List available model types
"""
import argparse
import json
import os
import sys
import stat
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
SRC_DIR = SCRIPT_DIR / "src"

# Source files in concatenation order
SOURCE_FILES = [
    "mod_constants.f90",
    "mod_tensor.f90",
    "mod_kinematics.f90",
    "mod_continuum.f90",
    "mod_hyperelastic.f90",
    "mod_icosahedron.f90",
    "mod_anisotropic.f90",
    "mod_network.f90",
    "mod_damage.f90",
    "mod_viscosity.f90",
    "umat_builder.f90",
    "uexternaldb.f90",
]

# Element-layer sources appended after SOURCE_FILES when emitting uel.f90
UEL_SOURCE_FILES = [
    "element/mod_uel_config.f90",
    "element/mod_uel_shape.f90",
    "element/mod_uel_element.f90",
    "element/uel_entry.f90",
]

# Defaults for the optional "element" config block (UEL emission)
DEFAULT_ELEMENT = {
    "type": "u3d8",      # 8-node brick (only type currently)
    "nint": 8,           # volume integration points (1 or 8)
    "fbar": True,        # F-bar locking treatment (active on the 8-pt brick)
    "num_elem": 1,       # UEL elements in the real mesh (sizes globalSdv)
    "elem_offset": 1000, # dummy-mesh element-number offset (must match deck)
}

# --- Example configurations ---------------------------------------------------

EXAMPLES = {
    "neo_hooke": {
        "name": "neo_hooke",
        "kbulk": 1000.0,
        "iso_type": 1, "iso_params": [10.0],
        "aniso_type": 0, "n_fiber_fam": 0, "aniso_params": [],
        "network_type": 0, "network_params": [],
        "damage_type": 0, "damage_params": [],
        "n_visco": 0, "visco_params": [],
        "test": {"stretch_max": 1.5, "gamma_max": 0.6, "nsteps": 400, "dtime": 0.01},
    },
    "mooney_rivlin": {
        "name": "mooney_rivlin",
        "kbulk": 1000.0,
        "iso_type": 2, "iso_params": [6.3, 0.012],
        "aniso_type": 0, "n_fiber_fam": 0, "aniso_params": [],
        "network_type": 0, "network_params": [],
        "damage_type": 0, "damage_params": [],
        "n_visco": 0, "visco_params": [],
        "test": {"stretch_max": 1.5, "gamma_max": 0.6, "nsteps": 400, "dtime": 0.01},
    },
    "humphrey_hgo": {
        "name": "humphrey_hgo",
        "kbulk": 500.0,
        "iso_type": 4, "iso_params": [2.0, 1.5],
        "aniso_type": 1, "n_fiber_fam": 1,
        "aniso_params": [100.0, 10.0, 0.226, 1.0, 0.0, 0.0],
        "network_type": 0, "network_params": [],
        "damage_type": 0, "damage_params": [],
        "n_visco": 0, "visco_params": [],
        "test": {"stretch_max": 1.3, "gamma_max": 0.4, "nsteps": 400, "dtime": 0.01},
    },
    "ogden_3term": {
        "name": "ogden_3term",
        "kbulk": 1000.0,
        "iso_type": 3,
        "iso_params": [3, 1.3, 5.0, 0.5, -2.0, 0.012, 2.0],
        "aniso_type": 0, "n_fiber_fam": 0, "aniso_params": [],
        "network_type": 0, "network_params": [],
        "damage_type": 0, "damage_params": [],
        "n_visco": 0, "visco_params": [],
        "test": {"stretch_max": 1.5, "gamma_max": 0.6, "nsteps": 400, "dtime": 0.01},
    },
    "neo_hooke_damage": {
        "name": "neo_hooke_damage",
        "kbulk": 1000.0,
        "iso_type": 1, "iso_params": [10.0],
        "aniso_type": 0, "n_fiber_fam": 0, "aniso_params": [],
        "network_type": 0, "network_params": [],
        "damage_type": 1, "damage_params": [5.0, 50.0],
        "n_visco": 0, "visco_params": [],
        "test": {"stretch_max": 1.5, "gamma_max": 0.6, "nsteps": 400, "dtime": 0.01},
    },
    "neo_hooke_visco": {
        "name": "neo_hooke_visco",
        "kbulk": 1000.0,
        "iso_type": 1, "iso_params": [10.0],
        "aniso_type": 0, "n_fiber_fam": 0, "aniso_params": [],
        "network_type": 0, "network_params": [],
        "damage_type": 0, "damage_params": [],
        "n_visco": 1, "visco_params": [0.5, 0.25],
        "test": {"stretch_max": 1.5, "gamma_max": 0.6, "nsteps": 400, "dtime": 0.01},
    },
    "affine_network": {
        "name": "affine_network",
        "kbulk": 500.0,
        "iso_type": 0, "iso_params": [],
        "aniso_type": 0, "n_fiber_fam": 0, "aniso_params": [],
        "network_type": 5, "network_params": [
            0.5, 1.0e6, 2.0, 1.0, 6, 1.0, 0.0, 0.0,
            1.0, 0.1, 0.01, 2.0, 0.1, 1.0, 0.0, 0.0,
        ],
        "damage_type": 0, "damage_params": [],
        "n_visco": 0, "visco_params": [],
        "test": {"stretch_max": 1.3, "gamma_max": 0.3, "nsteps": 200, "dtime": 0.01},
    },
    "humphrey_fiber": {
        "name": "humphrey_fiber",
        "kbulk": 500.0,
        "iso_type": 4, "iso_params": [2.0, 1.5],
        "aniso_type": 2, "n_fiber_fam": 1,
        "aniso_params": [100.0, 10.0, 1.0, 0.0, 0.0],
        "network_type": 0, "network_params": [],
        "damage_type": 0, "damage_params": [],
        "n_visco": 0, "visco_params": [],
        "test": {"stretch_max": 1.3, "gamma_max": 0.4, "nsteps": 400, "dtime": 0.01},
    },
    "humphrey_hgo_damage": {
        "name": "humphrey_hgo_damage",
        "kbulk": 500.0,
        "iso_type": 4, "iso_params": [2.0, 1.5],
        "aniso_type": 1, "n_fiber_fam": 1,
        "aniso_params": [100.0, 10.0, 0.226, 1.0, 0.0, 0.0],
        "network_type": 0, "network_params": [],
        "damage_type": 1, "damage_params": [5.0, 50.0],
        "n_visco": 0, "visco_params": [],
        "test": {"stretch_max": 1.3, "gamma_max": 0.4, "nsteps": 400, "dtime": 0.01},
    },
    "humphrey_fiber_damage": {
        "name": "humphrey_fiber_damage",
        "kbulk": 500.0,
        "iso_type": 4, "iso_params": [2.0, 1.5],
        "aniso_type": 2, "n_fiber_fam": 1,
        "aniso_params": [100.0, 10.0, 1.0, 0.0, 0.0],
        "network_type": 0, "network_params": [],
        "damage_type": 1, "damage_params": [5.0, 50.0],
        "n_visco": 0, "visco_params": [],
        "test": {"stretch_max": 1.3, "gamma_max": 0.4, "nsteps": 400, "dtime": 0.01},
    },
    "mooney_rivlin_visco": {
        "name": "mooney_rivlin_visco",
        "kbulk": 1000.0,
        "iso_type": 2, "iso_params": [6.3, 0.012],
        "aniso_type": 0, "n_fiber_fam": 0, "aniso_params": [],
        "network_type": 0, "network_params": [],
        "damage_type": 0, "damage_params": [],
        "n_visco": 1, "visco_params": [0.5, 0.25],
        "test": {"stretch_max": 1.5, "gamma_max": 0.6, "nsteps": 400, "dtime": 0.01},
    },
    "ogden_visco": {
        "name": "ogden_visco",
        "kbulk": 1000.0,
        "iso_type": 3,
        "iso_params": [3, 1.3, 5.0, 0.5, -2.0, 0.012, 2.0],
        "aniso_type": 0, "n_fiber_fam": 0, "aniso_params": [],
        "network_type": 0, "network_params": [],
        "damage_type": 0, "damage_params": [],
        "n_visco": 1, "visco_params": [0.5, 0.25],
        "test": {"stretch_max": 1.5, "gamma_max": 0.6, "nsteps": 400, "dtime": 0.01},
    },
    "humphrey_hgo_visco": {
        "name": "humphrey_hgo_visco",
        "kbulk": 500.0,
        "iso_type": 4, "iso_params": [2.0, 1.5],
        "aniso_type": 1, "n_fiber_fam": 1,
        "aniso_params": [100.0, 10.0, 0.226, 1.0, 0.0, 0.0],
        "network_type": 0, "network_params": [],
        "damage_type": 0, "damage_params": [],
        "n_visco": 1, "visco_params": [0.5, 0.25],
        "test": {"stretch_max": 1.3, "gamma_max": 0.4, "nsteps": 400, "dtime": 0.01},
    },
    "humphrey_fiber_visco": {
        "name": "humphrey_fiber_visco",
        "kbulk": 500.0,
        "iso_type": 4, "iso_params": [2.0, 1.5],
        "aniso_type": 2, "n_fiber_fam": 1,
        "aniso_params": [100.0, 10.0, 1.0, 0.0, 0.0],
        "network_type": 0, "network_params": [],
        "damage_type": 0, "damage_params": [],
        "n_visco": 1, "visco_params": [0.5, 0.25],
        "test": {"stretch_max": 1.3, "gamma_max": 0.4, "nsteps": 400, "dtime": 0.01},
    },
    "contractile_network": {
        "name": "contractile_network",
        "kbulk": 500.0,
        "iso_type": 0, "iso_params": [],
        "aniso_type": 0, "n_fiber_fam": 0, "aniso_params": [],
        "network_type": 4, "network_params": [
            0.5,                          # PHI
            0.2, 2.0, 1.0,               # N, B_orient, EFI
            11.0, 11.0,                   # FRIC, FFMAX
            6, 1.0, 0.0, 0.0,            # factor, prefdir
            0.988, 0.804, 38600.0, 0.438, # L, R0F, mu0, beta
            0.065, 1.007,                 # B0, lambda0
            0.014, 0.667,                 # R0C, ETAC
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,  # KCH(7) rate constants
        ],
        "damage_type": 0, "damage_params": [],
        "n_visco": 0, "visco_params": [],
        "test": {"stretch_max": 1.2, "gamma_max": 0.2, "nsteps": 100, "dtime": 0.01},
    },
    "mixed_network": {
        "name": "mixed_network",
        "kbulk": 500.0,
        "iso_type": 0, "iso_params": [],
        "aniso_type": 0, "n_fiber_fam": 0, "aniso_params": [],
        "network_type": 3, "network_params": [
            0.5,              # PHI
            5.0e5, 2.0,       # N_naff, PP
            5.0e5, 2.0, 1.0,  # N_aff, B_orient, EFI
            6, 1.0, 0.0, 0.0, # factor, prefdir
            1.0, 0.1, 0.01, 2.0, 0.1, 1.0, 0.0, 0.0,  # filprops(8)
        ],
        "damage_type": 0, "damage_params": [],
        "n_visco": 0, "visco_params": [],
        "test": {"stretch_max": 1.3, "gamma_max": 0.3, "nsteps": 200, "dtime": 0.01},
    },
    "humphrey_hgo_ai": {
        "name": "humphrey_hgo_ai",
        "kbulk": 500.0,
        "iso_type": 4, "iso_params": [2.0, 1.5],
        "aniso_type": 3, "n_fiber_fam": 1,
        "aniso_params": [100.0, 10.0, 5.0, 6, 1.0, 0.0, 0.0],
        "network_type": 0, "network_params": [],
        "damage_type": 0, "damage_params": [],
        "n_visco": 0, "visco_params": [],
        "test": {"stretch_max": 1.3, "gamma_max": 0.4, "nsteps": 200, "dtime": 0.01},
    },
    "affine_network_linkers": {
        "name": "affine_network_linkers",
        "kbulk": 500.0,
        "iso_type": 0, "iso_params": [],
        "aniso_type": 0, "n_fiber_fam": 0, "aniso_params": [],
        "network_type": 5, "network_params": [
            0.5, 1.0e6, 2.0, 1.0, 6, 1.0, 0.0, 0.0,
            1.0, 0.1, 0.01, 2.0, 0.1, 1.0, 0.014, 0.6667,
        ],
        "damage_type": 0, "damage_params": [],
        "n_visco": 0, "visco_params": [],
        "test": {"stretch_max": 1.3, "gamma_max": 0.3, "nsteps": 200, "dtime": 0.01},
    },
    "humphrey_muscle": {
        "name": "humphrey_muscle",
        "kbulk": 500.0,
        "iso_type": 4, "iso_params": [2.0, 1.5],
        "aniso_type": 5, "n_fiber_fam": 1,
        "aniso_params": [100.0, 10.0, 1.0, 0.0, 0.0, 50.0],
        "network_type": 0, "network_params": [],
        "damage_type": 0, "damage_params": [],
        "n_visco": 0, "visco_params": [],
        "test": {"stretch_max": 1.3, "gamma_max": 0.4, "nsteps": 400, "dtime": 0.01},
    },
    "nonaffine_network": {
        "name": "nonaffine_network",
        "kbulk": 500.0,
        "iso_type": 0, "iso_params": [],
        "aniso_type": 0, "n_fiber_fam": 0, "aniso_params": [],
        "network_type": 6, "network_params": [
            0.5, 1.0e6, 2.0, 6,
            1.0, 0.1, 0.01, 2.0, 0.1, 1.0, 0.0, 0.0,
        ],
        "damage_type": 0, "damage_params": [],
        "n_visco": 0, "visco_params": [],
        "test": {"stretch_max": 1.3, "gamma_max": 0.3, "nsteps": 200, "dtime": 0.01},
    },
}

# --- Helpers ------------------------------------------------------------------

def compute_nstatev(cfg):
    n = 1
    if cfg["damage_type"] > 0:
        n += 2
    n += 9 * cfg["n_visco"]
    if cfg["network_type"] == 4:  # contractile
        # Extract factor from network_params (position 6, 0-indexed)
        nparams = cfg["network_params"]
        factor = int(nparams[6])
        # Icosahedron: 20 faces, each subdivided into factor^2 subtriangles
        nwp = 20 * factor * factor
        n += 4 + nwp  # FRAC(4) + RU0(nwp)
    return n


def build_props(cfg):
    """Return flat list of all PROPS values."""
    p = [
        cfg["kbulk"],
        float(cfg["iso_type"]),
        float(cfg["aniso_type"]),
        float(cfg["n_fiber_fam"]),
        float(cfg["network_type"]),
        float(cfg["damage_type"]),
        float(cfg["n_visco"]),
    ]
    p.extend(cfg["iso_params"])
    for _ in range(cfg["n_fiber_fam"]):
        p.extend(cfg["aniso_params"])
    p.extend(cfg["network_params"])
    p.extend(cfg["damage_params"])
    p.extend(cfg["visco_params"])
    return p


def fmt_props_fortran(props):
    """Generate Fortran PROPS assignment lines."""
    lines = []
    for i, v in enumerate(props, 1):
        if v == int(v) and abs(v) < 1e10:
            lines.append(f"  props({i}) = {v:.1f}d0")
        else:
            lines.append(f"  props({i}) = {v:g}d0")
    return "\n".join(lines)


def fmt_props_abaqus(props):
    """Format PROPS values for ABAQUS *User Material card (8 per line)."""
    lines = []
    for i in range(0, len(props), 8):
        chunk = props[i : i + 8]
        line = ", ".join(f"{v:g}" for v in chunk)
        if i + 8 < len(props):
            line += ","
        lines.append(line)
    return "\n".join(lines)


# --- File generators ----------------------------------------------------------

def generate_umat_f90(outdir):
    """Concatenate all source modules into a single umat.f90."""
    out = outdir / "umat.f90"
    with open(out, "w") as f:
        f.write("! Auto-generated UMAT — do not edit manually.\n")
        f.write("! Regenerate with: python generate.py <config.json>\n\n")
        for src in SOURCE_FILES:
            path = SRC_DIR / src
            f.write(f"! {'='*72}\n")
            f.write(f"! SOURCE: {src}\n")
            f.write(f"! {'='*72}\n")
            f.write(path.read_text())
            f.write("\n\n")


def generate_aba_param(outdir):
    src = SRC_DIR / "aba_param.inc"
    (outdir / "aba_param.inc").write_text(src.read_text())


# --- UEL emission ---------------------------------------------------------

def element_cfg(cfg):
    """Merged element block (or None if UEL emission is not requested)."""
    if "element" not in cfg:
        return None
    elem = dict(DEFAULT_ELEMENT)
    elem.update(cfg["element"] or {})
    return elem


def subst_uel_config(text, elem):
    """Rewrite the parameter values in mod_uel_config.f90 from the element block."""
    import re
    text = re.sub(r"(numElem\s*=\s*)\d+", rf"\g<1>{elem['num_elem']}", text)
    text = re.sub(r"(ElemOffset\s*=\s*)\d+", rf"\g<1>{elem['elem_offset']}", text)
    fbar = ".true." if elem["fbar"] else ".false."
    text = re.sub(r"(use_fbar\s*=\s*)\.\w+\.", rf"\g<1>{fbar}", text)
    text = re.sub(r"(nIntPt\s*=\s*)\d+", rf"\g<1>{elem['nint']}", text)
    return text


def uel_source(elem):
    """Concatenated uel.f90 text: material modules + element layer."""
    parts = []
    for src in SOURCE_FILES + UEL_SOURCE_FILES:
        text = (SRC_DIR / src).read_text()
        if src.endswith("mod_uel_config.f90"):
            text = subst_uel_config(text, elem)
        parts.append(f"! {'='*72}\n! SOURCE: {src}\n! {'='*72}\n" + text)
    return ("! Auto-generated UEL (user element + material library) — do not edit.\n"
            "! Regenerate with: python generate.py <config.json>\n\n"
            + "\n\n".join(parts) + "\n")


def generate_uel_f90(outdir, cfg, elem):
    (outdir / "uel.f90").write_text(uel_source(elem))


# Unit-cube nodes in the element's local ordering (sketch in mod_uel_element):
# bottom face z=0: 1(0,0,0) 2(1,0,0) 3(1,1,0) 4(0,1,0); top face z=1: 5..8
UEL_CUBE_COORDS = [
    (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
    (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, 1.0), (0.0, 1.0, 1.0),
]


def generate_uel_test_driver(outdir, cfg, elem):
    """Standalone single-element driver: ramps the unit cube through affine
    deformation histories (same load cases as test_umat) and records the
    internal force on the x+ face (nodes 2,3,6,7)."""
    props = build_props(cfg)
    nprops = len(props)
    nstatev = compute_nstatev(cfg)
    nint = elem["nint"]
    test = cfg.get("test", {})
    stretch_max = test.get("stretch_max", 1.5)
    gamma_max = test.get("gamma_max", 0.6)
    nsteps = test.get("nsteps", 400)
    dtime = test.get("dtime", 0.01)

    coords_lines = "\n".join(
        f"  coords(1,{i+1}) = {x:.1f}d0; coords(2,{i+1}) = {y:.1f}d0; coords(3,{i+1}) = {z:.1f}d0"
        for i, (x, y, z) in enumerate(UEL_CUBE_COORDS))

    code = f"""\
! Auto-generated UEL test driver for material: {cfg["name"]}
! Drives one U3D8 element through affine deformation ramps (no ABAQUS).
! Output columns: step time, control parameter, x+ face force (fx, fy, fz).
program test_uel
  implicit none

  integer, parameter :: nnode = 8, ndofel = 24, mlvarx = 24, nrhs = 1
  integer, parameter :: nprops = {nprops}, nsdv = {nstatev}, nintp = {nint}
  integer, parameter :: nsvars = nintp*nsdv, njprop = 2
  integer, parameter :: mcrd = 3, jtype = 3, jelem = 1
  integer, parameter :: ndload = 0, mdload = 1, npredf = 1

  double precision :: rhs(mlvarx,1), amatrx(ndofel,ndofel), svars(nsvars)
  double precision :: energy(8), props(nprops), coords(mcrd,nnode)
  double precision :: u(ndofel), du(mlvarx,1), v(ndofel), a(ndofel)
  double precision :: time(2), dtime, params(1)
  double precision :: adlmag(mdload,1), ddlmag(mdload,1)
  double precision :: predef(2,npredf,nnode), pnewdt, period
  integer :: jdltyp(mdload,1), lflags(4), jprops(njprop)

  double precision :: Ftar(3,3)
  double precision :: stretch_max, gamma_max
  integer :: nsteps
  double precision, parameter :: zero = 0.0d0, one = 1.0d0

  ! --- Element/material setup ---
{fmt_props_fortran(props)}
  jprops(1) = nsdv   ! local SDVs per integration point
  jprops(2) = nsdv   ! global SDVs per integration point (UVARM)
  lflags = 0
  lflags(1) = 1      ! static general step
  lflags(2) = 1      ! nlgeom=yes
  dtime = {dtime}d0
  nsteps = {nsteps}
  stretch_max = {stretch_max}d0
  gamma_max = {gamma_max}d0
  v = zero; a = zero; params = zero; energy = zero
  adlmag = zero; ddlmag = zero; predef = zero; jdltyp = 0
  period = zero

{coords_lines}

  time = zero
  call uexternaldb(0, 0, time, zero, 0, 0)
  call execute_command_line('mkdir -p results')

  ! Uniaxial: F = diag(s, 1/sqrt(s), 1/sqrt(s))
  Ftar = zero
  Ftar(1,1) = stretch_max
  Ftar(2,2) = one/sqrt(stretch_max); Ftar(3,3) = one/sqrt(stretch_max)
  call run_case(Ftar, 'results/uel_uniaxial.dat', 'Uniaxial ')

  ! Biaxial: F = diag(s, s, 1/s^2)
  Ftar = zero
  Ftar(1,1) = stretch_max; Ftar(2,2) = stretch_max
  Ftar(3,3) = one/(stretch_max*stretch_max)
  call run_case(Ftar, 'results/uel_biaxial.dat', 'Biaxial  ')

  ! Pure shear: F12 = F21 = gamma
  Ftar = zero; Ftar(1,1) = one; Ftar(2,2) = one; Ftar(3,3) = one
  Ftar(1,2) = gamma_max; Ftar(2,1) = gamma_max
  call run_case(Ftar, 'results/uel_shear.dat', 'Shear    ')

  ! Simple shear: F12 = gamma
  Ftar = zero; Ftar(1,1) = one; Ftar(2,2) = one; Ftar(3,3) = one
  Ftar(1,2) = gamma_max
  call run_case(Ftar, 'results/uel_simple_shear.dat', 'Simple sh')

contains

  !> Ramp the element from I to Ftar in nsteps affine increments,
  !> carrying svars (state) across increments.
  subroutine run_case(Ft, fname, label)
    double precision, intent(in) :: Ft(3,3)
    character(*), intent(in) :: fname, label
    double precision :: F(3,3), uold(ndofel), fface(3), s
    integer :: i, n, k, kk, face_nodes(4)

    face_nodes = (/2, 3, 6, 7/)   ! x+ face (X=1)
    svars = zero; u = zero; uold = zero; time = zero
    pnewdt = one

    open(unit=21, file=fname, status='replace')
    do i = 1, nsteps
      s = dble(i)/dble(nsteps)
      F = identity3() + s*(Ft - identity3())
      ! Affine nodal displacements u_a = (F - I) X_a
      do n = 1, nnode
        do k = 1, 3
          u(3*(n-1)+k) = sum((F(k,:) - identity_row(k))*coords(:,n))
        end do
      end do
      du(:,1) = u - uold
      rhs = zero; amatrx = zero
      call uel(rhs, amatrx, svars, energy, ndofel, nrhs, nsvars, &
               props, nprops, coords, mcrd, nnode, u, du, v, a, jtype, &
               time, dtime, 1, i, jelem, params, ndload, jdltyp, adlmag, &
               predef, npredf, lflags, mlvarx, ddlmag, mdload, pnewdt, &
               jprops, njprop, period)
      ! Internal force on the x+ face: f = -sum(RHS) over face nodes
      fface = zero
      do kk = 1, 4
        n = face_nodes(kk)
        do k = 1, 3
          fface(k) = fface(k) - rhs(3*(n-1)+k, 1)
        end do
      end do
      write(21, '(5ES20.10)') time(1), s, fface(1), fface(2), fface(3)
      time(1) = time(1) + dtime
      uold = u
    end do
    close(21)
    write(*,'(A,A,A)') label, ' -> ', fname
  end subroutine run_case

  pure function identity3() result(iden)
    double precision :: iden(3,3)
    integer :: ii
    iden = 0.0d0
    do ii = 1, 3
      iden(ii,ii) = 1.0d0
    end do
  end function identity3

  pure function identity_row(k) result(row)
    integer, intent(in) :: k
    double precision :: row(3)
    row = 0.0d0
    row(k) = 1.0d0
  end function identity_row

end program test_uel

! Stub for ABAQUS-provided routine (standalone builds only)
subroutine getoutdir(outdir, lenoutdir)
  implicit none
  character(len=256), intent(out) :: outdir
  integer, intent(out) :: lenoutdir
  outdir = '.'
  lenoutdir = 1
end subroutine getoutdir
"""
    (outdir / "test_uel.f90").write_text(code)


def generate_uel_abaqus(outdir, cfg, elem):
    """ABAQUS single-element UEL deck: U3 real mesh + dummy mesh for UVARM
    visualization. Reuses the bcs_*.inp files written by generate_abaqus_dir
    (same node numbering and node sets)."""
    abq = outdir / "abaqus"
    abq.mkdir(exist_ok=True)

    props = build_props(cfg)
    nprops = len(props)
    nstatev = compute_nstatev(cfg)
    nvars = elem["nint"] * nstatev
    offset = elem["elem_offset"]
    dummy_type = "C3D8" if elem["nint"] == 8 else "C3D8R"

    # *UEL PROPERTY data: reals first, the two integer properties last
    uel_props = fmt_props_abaqus(list(props) + [nstatev, nstatev])

    deck = f"""\
*Heading
UEL single-element test — {cfg["name"]}
** Real mesh: user element U3 (8-node brick, F-bar={'on' if elem['fbar'] else 'off'}, {elem['nint']}-pt)
** Dummy mesh: {dummy_type} at element offset {offset}, carries UVARM output
*Node, nset=all_nodes
      1,           1.,           1.,           1.
      2,           1.,           0.,           1.
      3,           1.,           1.,           0.
      4,           1.,           0.,           0.
      5,           0.,           1.,           1.
      6,           0.,           0.,           1.
      7,           0.,           1.,           0.
      8,           0.,           0.,           0.
*User Element, type=U3, nodes=8, coordinates=3, properties={nprops}, iproperties=2, variables={nvars}, unsymm
1,2,3
*Element, type=U3, elset=main_element
1, 5, 6, 8, 7, 1, 2, 4, 3
*Element, type={dummy_type}, elset=dummy_mesh
{1 + offset}, 5, 6, 8, 7, 1, 2, 4, 3
*Nset, nset=Set-1, generate
 2,  8,  2
*Nset, nset=Set-2, generate
 1,  7,  2
*Nset, nset=Set-3
 1, 2, 5, 6
*Nset, nset=Set-4, generate
 5,  8,  1
*Nset, nset=Set-5
 2, 4, 6, 8
*Nset, nset=Set-6
 3, 4, 7, 8
*Nset, nset=Set-7, generate
 1, 4, 1
*Uel Property, elset=main_element
{uel_props}
*Solid Section, elset=dummy_mesh, material=dummy_material
*Material, name=dummy_material
*User output variables
{nstatev},
*Elastic
1.e-20
*Step, name=static, nlgeom=YES, unsymm=YES, inc=200
*Static
0.01, 1., 1e-05, 0.1
*INCLUDE, file=bcs_uni.inp
*OUTPUT,FIELD,VARIABLE=PRESELECT,FREQ=1
*ELEMENT OUTPUT, elset=dummy_mesh
UVARM
*OUTPUT,HISTORY,VARIABLE=PRESELECT,FREQ=1
*End Step
"""
    (abq / "uel_cube.inp").write_text(deck)

    run_sh = """\
#!/bin/bash
# Run ABAQUS single-element UEL test
# Usage: ./run_uel.sh [bcs_file]
BCS=${1:-bcs_uni.inp}
sed -i "s/INCLUDE, file=bcs_.*/INCLUDE, file=${BCS}/" uel_cube.inp
abaqus job=uel_cube user=../uel.f90 interactive
"""
    run_path = abq / "run_uel.sh"
    run_path.write_text(run_sh)
    run_path.chmod(run_path.stat().st_mode | stat.S_IEXEC)


def generate_test_driver(outdir, cfg):
    props = build_props(cfg)
    nprops = len(props)
    nstatev = compute_nstatev(cfg)
    test = cfg.get("test", {})
    stretch_max = test.get("stretch_max", 1.5)
    gamma_max = test.get("gamma_max", 0.6)
    nsteps = test.get("nsteps", 400)
    dtime = test.get("dtime", 0.01)

    code = f"""\
! Auto-generated test driver for material: {cfg["name"]}
! Runs uniaxial, biaxial, pure shear, and simple shear tests.
program test_umat
  implicit none

  integer, parameter :: ntens = 6, ndi = 3, nshr = 3
  integer, parameter :: nprops = {nprops}, nstatev = {nstatev}
  integer, parameter :: noel = 1, npt = 1

  double precision :: stress(ntens), statev(nstatev), ddsdde(ntens, ntens)
  double precision :: ddsddt(ntens), drplde(ntens)
  double precision :: stran(ntens), dstran(ntens)
  double precision :: time(2), predef(1), dpred(1)
  double precision :: props(nprops), coords(3), drot(3,3)
  double precision :: dfgrd0(3,3), dfgrd1(3,3)
  double precision :: sse, spd, scd, rpl, drpldt, dtime
  double precision :: temp, dtemp, pnewdt, celent
  character(len=8) :: cmname
  integer :: layer, kspt, kstep, kinc

  integer :: nsteps, i
  double precision :: stretch, dstretch, gamma_val, dgamma
  double precision, parameter :: zero = 0.0d0, one = 1.0d0

  ! --- Initialize all arrays ---
  stress = zero; statev = zero; ddsdde = zero
  stran = zero; dstran = zero; time = zero
  drot = zero; dfgrd0 = zero; dfgrd1 = zero
  coords = zero; predef = zero; dpred = zero
  temp = zero; dtemp = zero; pnewdt = one; celent = one
  rpl = zero; drpldt = zero; ddsddt = zero; drplde = zero
  layer = 1; kspt = 1; kinc = 1; cmname = 'UMAT'
  dfgrd0(1,1) = one; dfgrd0(2,2) = one; dfgrd0(3,3) = one

  ! --- Initialize external database (for RW network models) ---
  call uexternaldb(0, 0, time, zero, 0, 0)

  ! --- Material properties ---
{fmt_props_fortran(props)}

  ! --- Test parameters ---
  nsteps = {nsteps}
  dtime = {dtime}d0
  dstretch = ({stretch_max}d0 - one) / nsteps
  dgamma = (2.0d0 * {gamma_max}d0) / nsteps

  call execute_command_line('mkdir -p results')

  ! ========================== UNIAXIAL ==========================
  stress = zero; statev = zero; time = zero
  dfgrd1 = zero; dfgrd1(1,1) = one; dfgrd1(2,2) = one; dfgrd1(3,3) = one
  stretch = one
  open(unit=10, file='results/uniaxial.dat', status='replace')
  do i = 1, nsteps
    dfgrd1(1,1) = stretch
    dfgrd1(2,2) = one / sqrt(stretch)
    dfgrd1(3,3) = one / sqrt(stretch)
    call umat(stress, statev, ddsdde, sse, spd, scd, rpl, ddsddt, &
              drplde, drpldt, stran, dstran, time, dtime, temp, dtemp, &
              predef, dpred, cmname, ndi, nshr, ntens, nstatev, props, &
              nprops, coords, drot, pnewdt, celent, dfgrd0, dfgrd1, &
              noel, npt, layer, kspt, kstep, kinc)
    write(10, '(3ES20.10)') time(1), stretch, stress(1)
    time(1) = time(1) + dtime
    stretch = stretch + dstretch
  end do
  close(10)
  write(*, '(A)') 'Uniaxial  -> results/uniaxial.dat'

  ! ========================== BIAXIAL ==========================
  stress = zero; statev = zero; time = zero
  dfgrd1 = zero; dfgrd1(1,1) = one; dfgrd1(2,2) = one; dfgrd1(3,3) = one
  stretch = one
  open(unit=11, file='results/biaxial.dat', status='replace')
  do i = 1, nsteps
    dfgrd1(1,1) = stretch
    dfgrd1(2,2) = stretch
    dfgrd1(3,3) = one / (stretch * stretch)
    call umat(stress, statev, ddsdde, sse, spd, scd, rpl, ddsddt, &
              drplde, drpldt, stran, dstran, time, dtime, temp, dtemp, &
              predef, dpred, cmname, ndi, nshr, ntens, nstatev, props, &
              nprops, coords, drot, pnewdt, celent, dfgrd0, dfgrd1, &
              noel, npt, layer, kspt, kstep, kinc)
    write(11, '(3ES20.10)') time(1), stretch, stress(1)
    time(1) = time(1) + dtime
    stretch = stretch + dstretch
  end do
  close(11)
  write(*, '(A)') 'Biaxial   -> results/biaxial.dat'

  ! ========================== PURE SHEAR ==========================
  stress = zero; statev = zero; time = zero
  dfgrd1 = zero; dfgrd1(1,1) = one; dfgrd1(2,2) = one; dfgrd1(3,3) = one
  gamma_val = -{gamma_max}d0
  open(unit=12, file='results/shear.dat', status='replace')
  do i = 1, nsteps
    dfgrd1(1,2) = gamma_val
    dfgrd1(2,1) = gamma_val
    call umat(stress, statev, ddsdde, sse, spd, scd, rpl, ddsddt, &
              drplde, drpldt, stran, dstran, time, dtime, temp, dtemp, &
              predef, dpred, cmname, ndi, nshr, ntens, nstatev, props, &
              nprops, coords, drot, pnewdt, celent, dfgrd0, dfgrd1, &
              noel, npt, layer, kspt, kstep, kinc)
    write(12, '(3ES20.10)') time(1), gamma_val, stress(4)
    time(1) = time(1) + dtime
    gamma_val = gamma_val + dgamma
  end do
  close(12)
  write(*, '(A)') 'Shear     -> results/shear.dat'

  ! ======================== SIMPLE SHEAR ========================
  stress = zero; statev = zero; time = zero
  dfgrd1 = zero; dfgrd1(1,1) = one; dfgrd1(2,2) = one; dfgrd1(3,3) = one
  gamma_val = -{gamma_max}d0
  open(unit=13, file='results/simple_shear.dat', status='replace')
  do i = 1, nsteps
    dfgrd1(1,2) = gamma_val
    call umat(stress, statev, ddsdde, sse, spd, scd, rpl, ddsddt, &
              drplde, drpldt, stran, dstran, time, dtime, temp, dtemp, &
              predef, dpred, cmname, ndi, nshr, ntens, nstatev, props, &
              nprops, coords, drot, pnewdt, celent, dfgrd0, dfgrd1, &
              noel, npt, layer, kspt, kstep, kinc)
    write(13, '(3ES20.10)') time(1), gamma_val, stress(4)
    time(1) = time(1) + dtime
    gamma_val = gamma_val + dgamma
  end do
  close(13)
  write(*, '(A)') 'Simple sh -> results/simple_shear.dat'

end program test_umat

! Stub for ABAQUS-provided routine (standalone builds only)
subroutine getoutdir(outdir, lenoutdir)
  implicit none
  character(len=256), intent(out) :: outdir
  integer, intent(out) :: lenoutdir
  outdir = '.'
  lenoutdir = 1
end subroutine getoutdir
"""
    (outdir / "test_umat.f90").write_text(code)


def generate_makefile(outdir, elem=None):
    uel_all = " test_uel" if elem else ""
    mk = f"""\
FC      = gfortran
FFLAGS  = -O2 -ffree-form

all: test_umat{uel_all}

test_umat: umat.f90 test_umat.f90
\t$(FC) $(FFLAGS) -o $@ $^

run: test_umat
\t./test_umat
"""
    if elem:
        mk += """
test_uel: uel.f90 test_uel.f90
\t$(FC) $(FFLAGS) -o $@ $^

run_uel: test_uel
\t./test_uel
"""
    mk += """
clean:
\trm -f test_umat test_uel *.mod
\trm -rf results
"""
    (outdir / "Makefile").write_text(mk)


def generate_abaqus_dir(outdir, cfg):
    abq = outdir / "abaqus"
    abq.mkdir(exist_ok=True)

    props = build_props(cfg)
    nprops = len(props)
    nstatev = compute_nstatev(cfg)

    # --- cube.inp ---
    cube = """\
** Unit cube, single C3D8 element
*Node, nset=all_nodes
      1,           1.,           1.,           1.
      2,           1.,           0.,           1.
      3,           1.,           1.,           0.
      4,           1.,           0.,           0.
      5,           0.,           1.,           1.
      6,           0.,           0.,           1.
      7,           0.,           1.,           0.
      8,           0.,           0.,           0.
*Element, type=C3D8, elset=main_element
1, 5, 6, 8, 7, 1, 2, 4, 3
*Nset, nset=Set-1, generate
 2,  8,  2
*Nset, nset=Set-2, generate
 1,  7,  2
*Nset, nset=Set-3
 1, 2, 5, 6
*Nset, nset=Set-4, generate
 5,  8,  1
*Nset, nset=Set-5
 2, 4, 6, 8
*Nset, nset=Set-6
 3, 4, 7, 8
*Nset, nset=Set-7, generate
 1, 4, 1
*Nset, nset=Set-8
 1,3
*Nset, nset=Set-9
 6,8
*Elset, elset=Surf
 1,
*Surface, type=ELEMENT, name=Surf-1
Surf, S1
*Surface, type=ELEMENT, name=Surf-2
Surf, S2
*Surface, type=ELEMENT, name=Surf-3
Surf, S3
*Surface, type=ELEMENT, name=Surf-4
Surf, S4
*Surface, type=ELEMENT, name=Surf-5
Surf, S5
*Surface, type=ELEMENT, name=Surf-6
Surf, S6
*INCLUDE, file=sec.inp
*Step, name=static, nlgeom=Yes, inc=200
*Static
0.01, 1., 1e-05, 0.1
*INCLUDE, file=bcs_uni.inp
*OUTPUT,FIELD,VARIABLE=PRESELECT,FREQ=1
*ELEMENT OUTPUT, elset=main_element
SDV
*OUTPUT,HISTORY,VARIABLE=PRESELECT,FREQ=1
*End Step
"""
    (abq / "cube.inp").write_text(cube)

    # --- sec.inp ---
    sdv_lines = f"{nstatev},\n1, DET, \"DET\""
    if cfg["damage_type"] > 0:
        sdv_lines += "\n2, DMG, \"damage\"\n3, MAXSEF, \"max SEF\""
    sec = f"""\
*Solid Section, elset=main_element, material=UD
*Material, name=UD
*User Material, constants={nprops}
{fmt_props_abaqus(props)}
*DEPVAR
{sdv_lines}
"""
    (abq / "sec.inp").write_text(sec)

    # --- bcs_uni.inp ---
    bcs_uni = """\
*Boundary
Set-3, ZSYMM
*Boundary
Set-4, XSYMM
*Boundary
Set-1, YSYMM
*Boundary, type=displacement
Set-2, 2,2, 0.6
"""
    (abq / "bcs_uni.inp").write_text(bcs_uni)

    # --- bcs_bi.inp ---
    bcs_bi = """\
*Boundary
Set-1, YSYMM
*Boundary
Set-3, ZSYMM
*Boundary
Set-4, XSYMM
*Boundary
Set-2, 2,2, 0.6
*Boundary
Set-7, 1,1, 0.6
"""
    (abq / "bcs_bi.inp").write_text(bcs_bi)

    # --- bcs_sh.inp ---
    bcs_sh = """\
*Boundary
Set-3, ZSYMM
*Boundary
Set-1, 1,2, 0.0
*Boundary
Set-2, 1,1, 0.6
Set-2, 2,2, 0.0
"""
    (abq / "bcs_sh.inp").write_text(bcs_sh)

    # --- run.sh ---
    run_sh = """\
#!/bin/bash
# Run ABAQUS single-element test
# Usage: ./run.sh [bcs_file]
#   ./run.sh                  # uniaxial (default)
#   ./run.sh bcs_bi.inp       # biaxial
#   ./run.sh bcs_sh.inp       # shear
BCS=${1:-bcs_uni.inp}
# Update the included BCS file
sed -i "s/INCLUDE, file=bcs_.*/INCLUDE, file=${BCS}/" cube.inp
abaqus job=cube user=../umat.f90 interactive
"""
    run_path = abq / "run.sh"
    run_path.write_text(run_sh)
    run_path.chmod(run_path.stat().st_mode | stat.S_IEXEC)


# --- Main ---------------------------------------------------------------------

def list_models():
    print("""
Available model types:

  ISO_TYPE (isotropic):
    0  None
    1  Neo-Hookean          params: C10
    2  Mooney-Rivlin        params: C10, C01
    3  Ogden (N-term)       params: N, mu1, alpha1, ..., muN, alphaN
    4  Humphrey exponential params: C10, C01

  ANISO_TYPE (anisotropic, per fiber family):
    0  None
    1  HGO (dispersed)      params: K1, K2, kappa, fiber_x, fiber_y, fiber_z
    2  Humphrey fiber       params: K1, K2, fiber_x, fiber_y, fiber_z
    3  HGO (AI discrete)    params: K1, K2, bdisp, factor, fiber_x, fiber_y, fiber_z
    4  Humphrey (AI discrete) params: K1, K2, bdisp, factor, fiber_x, fiber_y, fiber_z
    5  Humphrey + activation params: K1, K2, fiber_x, fiber_y, fiber_z, T0M

  NETWORK_TYPE:
    0  None
    1  Affine (RW)          params: PHI, N, B_orient, EFI, pdir(3),
                                    L, R0F, mu0, beta, B0, lambda0, R0C, ETAC
    2  Non-affine (RW)      params: PHI, N, B_orient, EFI, PP,
                                    L, R0F, mu0, beta, B0, lambda0, R0C, ETAC
    3  Mixed (AI)           params: PHI, N_naff, PP, N_aff, B_orient, EFI, factor, pdir(3),
                                    L, R0F, mu0, beta, B0, lambda0, R0C, ETAC
    4  Contractile (AI)     params: PHI, N, B_orient, EFI, FRIC, FFMAX, factor, pdir(3),
                                    L, R0F, mu0, beta, B0, lambda0, R0C, ETAC, KCH(7)
    5  Affine (AI)          params: PHI, N, B_orient, EFI, factor, pdir(3),
                                    L, R0F, mu0, beta, B0, lambda0, R0C, ETAC
    6  Non-affine (AI)      params: PHI, N, PP, factor,
                                    L, R0F, mu0, beta, B0, lambda0, R0C, ETAC

  DAMAGE_TYPE:
    0  None
    1  Sigmoid              params: beta_d, psi_half

  VISCO (per Maxwell branch, max 3):
    params: tau, theta  (per branch)

Example configs: """ + ", ".join(EXAMPLES.keys()))


def sphere_quadrature(n=60):
    """Quasi-uniform unit-sphere quadrature in the format UEXTERNALDB reads
    ('x y z weight' per line, weights = 1/n). Fibonacci spiral; a generic
    orientation quadrature for the RW network types (1, 2). Replace with a
    higher-order spherical design for production accuracy if needed."""
    import math
    ga = math.pi * (3.0 - math.sqrt(5.0))
    w = 1.0 / n
    lines = []
    for i in range(n):
        z = 1.0 - 2.0 * (i + 0.5) / n
        r = math.sqrt(max(0.0, 1.0 - z * z))
        th = ga * i
        lines.append(f"{r*math.cos(th):.13f} {r*math.sin(th):.13f} {z:.13f} {w:.13f}")
    return "\n".join(lines) + "\n"


def generate_quadrature(outdir, cfg):
    """Ship a sphere quadrature for the RW network types (1, 2), which read it
    via UEXTERNALDB. AI types (3-6) integrate over an icosahedron and need none."""
    if cfg["network_type"] in (1, 2):
        content = sphere_quadrature(60)
        (outdir / "sphere_int60c.inp").write_text(content)
        (outdir / "abaqus" / "sphere_int60c.inp").write_text(content)


def generate(cfg):
    name = cfg["name"]
    outdir = SCRIPT_DIR / name
    outdir.mkdir(exist_ok=True)

    elem = element_cfg(cfg)

    generate_umat_f90(outdir)
    generate_aba_param(outdir)
    generate_test_driver(outdir, cfg)
    generate_makefile(outdir, elem)
    generate_abaqus_dir(outdir, cfg)
    generate_quadrature(outdir, cfg)
    if elem:
        generate_uel_f90(outdir, cfg, elem)
        generate_uel_test_driver(outdir, cfg, elem)
        generate_uel_abaqus(outdir, cfg, elem)

    # Save the config for reproducibility
    (outdir / "config.json").write_text(json.dumps(cfg, indent=2) + "\n")

    nprops = len(build_props(cfg))
    nstatev = compute_nstatev(cfg)

    print(f"\nGenerated material: {name}/")
    print(f"  umat.f90          Concatenated UMAT source ({nprops} PROPS, {nstatev} STATEV)")
    print(f"  test_umat.f90     Standalone test driver")
    print(f"  Makefile          Build & run: make run")
    print(f"  aba_param.inc     ABAQUS include")
    print(f"  config.json       Configuration (for regeneration)")
    print(f"  abaqus/           ABAQUS single-element test")
    print(f"    cube.inp        C3D8 mesh + step definition")
    print(f"    sec.inp         *User Material card")
    print(f"    bcs_uni.inp     Uniaxial boundary conditions")
    print(f"    bcs_bi.inp      Biaxial boundary conditions")
    print(f"    bcs_sh.inp      Shear boundary conditions")
    print(f"    run.sh          ABAQUS submission script")
    if elem:
        nvars = elem["nint"] * nstatev
        print(f"  uel.f90           User element + material library "
              f"(U3, {elem['nint']}-pt, F-bar={'on' if elem['fbar'] else 'off'}, Variables={nvars})")
        print(f"  test_uel.f90      Standalone single-element driver")
        print(f"    abaqus/uel_cube.inp + run_uel.sh   ABAQUS UEL test")

    print(f"\nStandalone test:  cd {name} && make run")
    if elem:
        print(f"UEL test:         cd {name} && make run_uel")
    print(f"ABAQUS test:      cd {name}/abaqus && ./run.sh")


def main():
    parser = argparse.ArgumentParser(
        description="Generate a self-contained UMAT material law directory."
    )
    parser.add_argument("config", nargs="?", help="JSON configuration file")
    parser.add_argument("--example", metavar="NAME",
                        help="Generate from built-in example: " + ", ".join(EXAMPLES.keys()))
    parser.add_argument("--list", action="store_true", help="List available model types")
    parser.add_argument("--uel", action="store_true",
                        help="Also emit a user element (uel.f90 + test_uel.f90 + ABAQUS "
                             "UEL deck) driven by this material; configs may instead "
                             "carry an explicit \"element\" block")
    args = parser.parse_args()

    if args.list:
        list_models()
        return

    if args.example:
        if args.example not in EXAMPLES:
            print(f"Unknown example: {args.example}")
            print(f"Available: {', '.join(EXAMPLES.keys())}")
            sys.exit(1)
        cfg = dict(EXAMPLES[args.example])
        if args.uel:
            cfg.setdefault("element", {})
        generate(cfg)
        return

    if args.config:
        with open(args.config) as f:
            cfg = json.load(f)
        if args.uel:
            cfg.setdefault("element", {})
        generate(cfg)
        return

    parser.print_help()


if __name__ == "__main__":
    main()

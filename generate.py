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


def generate_makefile(outdir):
    mk = """\
FC      = gfortran
FFLAGS  = -O2 -ffree-form

all: test_umat

test_umat: umat.f90 test_umat.f90
\t$(FC) $(FFLAGS) -o $@ $^

run: test_umat
\t./test_umat

clean:
\trm -f test_umat *.mod
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


def generate(cfg):
    name = cfg["name"]
    outdir = SCRIPT_DIR / name
    outdir.mkdir(exist_ok=True)

    generate_umat_f90(outdir)
    generate_aba_param(outdir)
    generate_test_driver(outdir, cfg)
    generate_makefile(outdir)
    generate_abaqus_dir(outdir, cfg)

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
    print(f"\nStandalone test:  cd {name} && make run")
    print(f"ABAQUS test:      cd {name}/abaqus && ./run.sh")


def main():
    parser = argparse.ArgumentParser(
        description="Generate a self-contained UMAT material law directory."
    )
    parser.add_argument("config", nargs="?", help="JSON configuration file")
    parser.add_argument("--example", metavar="NAME",
                        help="Generate from built-in example: " + ", ".join(EXAMPLES.keys()))
    parser.add_argument("--list", action="store_true", help="List available model types")
    args = parser.parse_args()

    if args.list:
        list_models()
        return

    if args.example:
        if args.example not in EXAMPLES:
            print(f"Unknown example: {args.example}")
            print(f"Available: {', '.join(EXAMPLES.keys())}")
            sys.exit(1)
        cfg = EXAMPLES[args.example]
        generate(cfg)
        return

    if args.config:
        with open(args.config) as f:
            cfg = json.load(f)
        generate(cfg)
        return

    parser.print_help()


if __name__ == "__main__":
    main()

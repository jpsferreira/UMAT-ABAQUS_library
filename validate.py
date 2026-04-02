#!/usr/bin/env python3
"""Validate new modular UMAT against legacy monolithic UMATs.

For each model pair, generates Fortran test drivers, compiles with the
respective UMAT source, runs both, and compares stress output.

Usage:
    uv run python validate.py
"""

import os
import shutil
import subprocess
import sys
import tempfile
import textwrap
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent

# ---------------------------------------------------------------------------
# Model pair definitions
# ---------------------------------------------------------------------------

PAIRS = {
    "neo_hooke": {
        "legacy_umat": "soft_tissues/1_neo_hooke/umat_nh.for",
        "legacy_inc": "soft_tissues/1_neo_hooke",
        "legacy_nprops": 9,
        "legacy_nstatev": 1,
        "legacy_props": [1000.0, 10.0, 0, 2.0, 0.835, 1.2, 7.0, 12.0, 2.0],
        "new_example": "neo_hooke",
        "new_nprops": 8,
        "new_nstatev": 1,
        "new_props": [1000.0, 1, 0, 0, 0, 0, 0, 10.0],
    },
    "mooney_rivlin": {
        "legacy_umat": "soft_tissues/2_mooney_rivlin/umat_mr.for",
        "legacy_inc": "soft_tissues/2_mooney_rivlin",
        "legacy_nprops": 10,
        "legacy_nstatev": 1,
        "legacy_props": [1000.0, 6.3, 0.012, 0, 2.0, 0.835, 1.2, 7.0, 12.0, 2.0],
        "new_example": "mooney_rivlin",
        "new_nprops": 9,
        "new_nstatev": 1,
        "new_props": [1000.0, 2, 0, 0, 0, 0, 0, 6.3, 0.012],
    },
    "neo_hooke_damage": {
        "legacy_umat": "soft_tissues/with_damage/1_neo_hooke_d/umat_nh_d.for",
        "legacy_inc": "soft_tissues/with_damage/1_neo_hooke_d",
        "legacy_nprops": 4,
        "legacy_nstatev": 3,
        "legacy_props": [1000.0, 10.0, 5.0, 50.0],
        "new_example": "neo_hooke_damage",
        "new_nprops": 10,
        "new_nstatev": 3,
        "new_props": [1000.0, 1, 0, 0, 0, 1, 0, 10.0, 5.0, 50.0],
    },
    "neo_hooke_visco": {
        "legacy_umat": "soft_tissues/with_visco/1_neo_hooke_v/umat_nh_v.for",
        "legacy_inc": "soft_tissues/with_visco/1_neo_hooke_v",
        "legacy_nprops": 9,
        "legacy_nstatev": 28,
        "legacy_props": [1000.0, 10.0, 1, 0.5, 0.25, 1.2, 7.0, 12.0, 2.0],
        "new_example": "neo_hooke_visco",
        "new_nprops": 10,
        "new_nstatev": 10,
        "new_props": [1000.0, 1, 0, 0, 0, 0, 1, 10.0, 0.5, 0.25],
    },
    # --- Ogden (legacy has KBULK as LAST prop, Ogden params first) ---
    "ogden": {
        "legacy_umat": "soft_tissues/3_ogden/umat_og.for",
        "legacy_inc": "soft_tissues/3_ogden",
        "legacy_nprops": 7,
        "legacy_nstatev": 1,
        "legacy_props": [1.3, 5.0, 0.5, -2.0, 0.012, 2.0, 1000.0],
        "new_example": "ogden_3term",
        "new_nprops": 15,
        "new_nstatev": 1,
        "new_props": [1000.0, 3, 0, 0, 0, 0, 0, 3, 1.3, 5.0, 0.5, -2.0, 0.012, 2.0],
    },
    # --- Mooney-Rivlin visco ---
    "mooney_rivlin_visco": {
        "legacy_umat": "soft_tissues/with_visco/2_mooney_rivlin_v/umat_mr_v.for",
        "legacy_inc": "soft_tissues/with_visco/2_mooney_rivlin_v",
        "legacy_nprops": 10,
        "legacy_nstatev": 28,
        "legacy_props": [1000.0, 6.3, 0.012, 1, 0.5, 0.25, 1.2, 7.0, 12.0, 2.0],
        "new_example": "mooney_rivlin_visco",
        "new_nprops": 11,
        "new_nstatev": 10,
        "new_props": [1000.0, 2, 0, 0, 0, 0, 1, 6.3, 0.012, 0.5, 0.25],
    },
    # --- Ogden visco (legacy: [mu1,a1,mu2,a2,mu3,a3,KBULK,V,tau1,theta1,...]) ---
    "ogden_visco": {
        "legacy_umat": "soft_tissues/with_visco/3_ogden_v/umat_og_v.for",
        "legacy_inc": "soft_tissues/with_visco/3_ogden_v",
        "legacy_nprops": 14,
        "legacy_nstatev": 28,
        "legacy_props": [1.3, 5.0, 0.5, -2.0, 0.012, 2.0, 1000.0, 1, 0.5, 0.25, 1.2, 7.0, 12.0, 2.0],
        "new_example": "ogden_visco",
        "new_nprops": 17,
        "new_nstatev": 10,
        "new_props": [1000.0, 3, 0, 0, 0, 0, 1, 3, 1.3, 5.0, 0.5, -2.0, 0.012, 2.0, 0.5, 0.25],
    },
}

# Test sweep parameters
NSTEPS = 400
STRETCH_MAX = 1.5
GAMMA_MAX = 0.6
DTIME = 0.01


def fmt_props(props):
    """Format a props list as Fortran assignments."""
    lines = []
    for i, v in enumerate(props, 1):
        if isinstance(v, float):
            lines.append(f"  props({i}) = {v}d0")
        else:
            lines.append(f"  props({i}) = {float(v)}d0")
    return "\n".join(lines)


def make_driver(nprops, nstatev, props, results_dir, call_uexternaldb=False):
    """Generate a Fortran test driver source string."""
    props_code = fmt_props(props)
    uextdb_call = "  call uexternaldb(0, 0, time, zero, 0, 0)" if call_uexternaldb else ""
    getoutdir_stub = textwrap.dedent("""\
        subroutine getoutdir(outdir, lenoutdir)
          implicit none
          character(len=256), intent(out) :: outdir
          integer, intent(out) :: lenoutdir
          outdir = '.'
          lenoutdir = 1
        end subroutine getoutdir
    """) if call_uexternaldb else ""

    return textwrap.dedent(f"""\
        program validate_driver
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

          stress = zero; statev = zero; ddsdde = zero
          stran = zero; dstran = zero; time = zero
          drot = zero; dfgrd0 = zero; dfgrd1 = zero
          coords = zero; predef = zero; dpred = zero
          temp = zero; dtemp = zero; pnewdt = one; celent = one
          rpl = zero; drpldt = zero; ddsddt = zero; drplde = zero
          layer = 1; kspt = 1; kinc = 1; cmname = 'UMAT'
          dfgrd0(1,1) = one; dfgrd0(2,2) = one; dfgrd0(3,3) = one

        {uextdb_call}

        {props_code}

          nsteps = {NSTEPS}
          dtime = {DTIME}d0
          dstretch = ({STRETCH_MAX}d0 - one) / nsteps
          dgamma = (2.0d0 * {GAMMA_MAX}d0) / nsteps

          call execute_command_line('mkdir -p {results_dir}')

          ! ========================== UNIAXIAL ==========================
          stress = zero; statev = zero; time = zero
          dfgrd1 = zero; dfgrd1(1,1) = one; dfgrd1(2,2) = one; dfgrd1(3,3) = one
          stretch = one
          open(unit=10, file='{results_dir}/uniaxial.dat', status='replace')
          do i = 1, nsteps
            dfgrd1(1,1) = stretch
            dfgrd1(2,2) = one / sqrt(stretch)
            dfgrd1(3,3) = one / sqrt(stretch)
            call umat(stress, statev, ddsdde, sse, spd, scd, rpl, ddsddt, &
                      drplde, drpldt, stran, dstran, time, dtime, temp, dtemp, &
                      predef, dpred, cmname, ndi, nshr, ntens, nstatev, props, &
                      nprops, coords, drot, pnewdt, celent, dfgrd0, dfgrd1, &
                      noel, npt, layer, kspt, kstep, kinc)
            write(10, '(8ES20.10)') time(1), stretch, stress(1), stress(2), &
                                     stress(3), stress(4), stress(5), stress(6)
            time(1) = time(1) + dtime
            stretch = stretch + dstretch
          end do
          close(10)

          ! ========================== BIAXIAL ==========================
          stress = zero; statev = zero; time = zero
          dfgrd1 = zero; dfgrd1(1,1) = one; dfgrd1(2,2) = one; dfgrd1(3,3) = one
          stretch = one
          open(unit=11, file='{results_dir}/biaxial.dat', status='replace')
          do i = 1, nsteps
            dfgrd1(1,1) = stretch
            dfgrd1(2,2) = stretch
            dfgrd1(3,3) = one / (stretch * stretch)
            call umat(stress, statev, ddsdde, sse, spd, scd, rpl, ddsddt, &
                      drplde, drpldt, stran, dstran, time, dtime, temp, dtemp, &
                      predef, dpred, cmname, ndi, nshr, ntens, nstatev, props, &
                      nprops, coords, drot, pnewdt, celent, dfgrd0, dfgrd1, &
                      noel, npt, layer, kspt, kstep, kinc)
            write(11, '(8ES20.10)') time(1), stretch, stress(1), stress(2), &
                                     stress(3), stress(4), stress(5), stress(6)
            time(1) = time(1) + dtime
            stretch = stretch + dstretch
          end do
          close(11)

          ! ========================== PURE SHEAR ==========================
          stress = zero; statev = zero; time = zero
          dfgrd1 = zero; dfgrd1(1,1) = one; dfgrd1(2,2) = one; dfgrd1(3,3) = one
          gamma_val = -{GAMMA_MAX}d0
          open(unit=12, file='{results_dir}/shear.dat', status='replace')
          do i = 1, nsteps
            dfgrd1(1,2) = gamma_val
            dfgrd1(2,1) = gamma_val
            call umat(stress, statev, ddsdde, sse, spd, scd, rpl, ddsddt, &
                      drplde, drpldt, stran, dstran, time, dtime, temp, dtemp, &
                      predef, dpred, cmname, ndi, nshr, ntens, nstatev, props, &
                      nprops, coords, drot, pnewdt, celent, dfgrd0, dfgrd1, &
                      noel, npt, layer, kspt, kstep, kinc)
            write(12, '(8ES20.10)') time(1), gamma_val, stress(1), stress(2), &
                                     stress(3), stress(4), stress(5), stress(6)
            time(1) = time(1) + dtime
            gamma_val = gamma_val + dgamma
          end do
          close(12)

          ! ======================== SIMPLE SHEAR ========================
          stress = zero; statev = zero; time = zero
          dfgrd1 = zero; dfgrd1(1,1) = one; dfgrd1(2,2) = one; dfgrd1(3,3) = one
          gamma_val = -{GAMMA_MAX}d0
          open(unit=13, file='{results_dir}/simple_shear.dat', status='replace')
          do i = 1, nsteps
            dfgrd1(1,2) = gamma_val
            call umat(stress, statev, ddsdde, sse, spd, scd, rpl, ddsddt, &
                      drplde, drpldt, stran, dstran, time, dtime, temp, dtemp, &
                      predef, dpred, cmname, ndi, nshr, ntens, nstatev, props, &
                      nprops, coords, drot, pnewdt, celent, dfgrd0, dfgrd1, &
                      noel, npt, layer, kspt, kstep, kinc)
            write(13, '(8ES20.10)') time(1), gamma_val, stress(1), stress(2), &
                                     stress(3), stress(4), stress(5), stress(6)
            time(1) = time(1) + dtime
            gamma_val = gamma_val + dgamma
          end do
          close(13)

        end program validate_driver
        {getoutdir_stub}
    """)


def compile_and_run(name, side, umat_path, driver_src, inc_dir, workdir, fflags):
    """Write driver, compile with UMAT, and run."""
    driver_file = workdir / f"driver_{side}.f90"
    driver_file.write_text(driver_src)

    exe = workdir / f"run_{side}"
    cmd = ["gfortran"] + fflags + ["-o", str(exe)]
    if inc_dir:
        cmd += ["-I", str(inc_dir)]
    cmd += [str(umat_path), str(driver_file)]

    result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(workdir))
    if result.returncode != 0:
        print(f"  COMPILE ERROR ({side}):")
        print(result.stderr)
        return False

    result = subprocess.run([str(exe)], capture_output=True, text=True, cwd=str(workdir))
    if result.returncode != 0:
        print(f"  RUNTIME ERROR ({side}):")
        print(result.stderr)
        return False

    return True


def compare_dat(legacy_path, new_path, label, atol=1e-8, rtol=1e-6):
    """Compare two .dat files. Returns True if they match within tolerance."""
    if not legacy_path.exists():
        print(f"    {label}: SKIP (legacy file missing)")
        return None
    if not new_path.exists():
        print(f"    {label}: SKIP (new file missing)")
        return None

    leg = np.loadtxt(str(legacy_path))
    new = np.loadtxt(str(new_path))

    if leg.shape != new.shape:
        print(f"    {label}: FAIL (shape mismatch: {leg.shape} vs {new.shape})")
        return False

    # Compare stress columns (columns 2-7, 0-indexed)
    stress_leg = leg[:, 2:]
    stress_new = new[:, 2:]

    if np.allclose(stress_leg, stress_new, atol=atol, rtol=rtol):
        max_err = np.max(np.abs(stress_leg - stress_new))
        print(f"    {label}: PASS (max err = {max_err:.2e})")
        return True
    else:
        diff = np.abs(stress_leg - stress_new)
        max_err = diff.max()
        # Find which component has max error
        idx = np.unravel_index(diff.argmax(), diff.shape)
        comp_names = ["S11", "S22", "S33", "S12", "S13", "S23"]
        print(f"    {label}: FAIL (max err = {max_err:.2e} at step {idx[0]}, "
              f"{comp_names[idx[1]]}, "
              f"legacy={stress_leg[idx]:.10e}, new={stress_new[idx]:.10e})")
        return False


def validate_pair(name, cfg):
    """Validate one legacy-vs-new pair."""
    print(f"\n{'='*60}")
    print(f"  {name}")
    print(f"{'='*60}")

    legacy_umat = ROOT / cfg["legacy_umat"]
    legacy_inc = ROOT / cfg["legacy_inc"]
    new_example_dir = ROOT / cfg["new_example"]
    new_umat = new_example_dir / "umat.f90"

    # Check legacy exists
    if not legacy_umat.exists():
        print(f"  SKIP: legacy UMAT not found at {legacy_umat}")
        return None

    # Generate new example if needed
    if not new_umat.exists():
        print(f"  Generating new example: {cfg['new_example']}")
        subprocess.run(
            [sys.executable, str(ROOT / "generate.py"), "--example", cfg["new_example"]],
            cwd=str(ROOT), check=True, capture_output=True,
        )

    if not new_umat.exists():
        print(f"  SKIP: could not generate new UMAT")
        return None

    # Work in a temp directory
    workdir = ROOT / "validation_tmp" / name
    workdir.mkdir(parents=True, exist_ok=True)

    # Copy aba_param.inc for legacy compilation
    aba_param = legacy_inc / "aba_param.inc"
    if aba_param.exists():
        shutil.copy2(aba_param, workdir / "aba_param.inc")
    param_umat = legacy_inc / "param_umat.inc"
    if param_umat.exists():
        shutil.copy2(param_umat, workdir / "param_umat.inc")

    # Generate and run legacy driver
    legacy_driver = make_driver(
        cfg["legacy_nprops"], cfg["legacy_nstatev"], cfg["legacy_props"],
        results_dir="results_legacy", call_uexternaldb=False,
    )
    ok_leg = compile_and_run(
        name, "legacy", legacy_umat, legacy_driver, workdir, workdir,
        fflags=["-O2"],
    )

    # Generate and run new driver
    new_driver = make_driver(
        cfg["new_nprops"], cfg["new_nstatev"], cfg["new_props"],
        results_dir="results_new", call_uexternaldb=True,
    )
    ok_new = compile_and_run(
        name, "new", new_umat, new_driver, None, workdir,
        fflags=["-O2", "-ffree-form"],
    )

    if not ok_leg or not ok_new:
        return False

    # Compare all load cases
    load_cases = ["uniaxial", "biaxial", "shear", "simple_shear"]
    results = []
    for lc in load_cases:
        leg_dat = workdir / "results_legacy" / f"{lc}.dat"
        new_dat = workdir / "results_new" / f"{lc}.dat"
        r = compare_dat(leg_dat, new_dat, lc)
        results.append(r)

    return all(r for r in results if r is not None)


def main():
    print("UMAT Validation: Legacy vs New Modular Code")
    print(f"Sweep: {NSTEPS} steps, stretch to {STRETCH_MAX}, gamma to {GAMMA_MAX}")

    all_pass = True
    summary = {}

    for name, cfg in PAIRS.items():
        result = validate_pair(name, cfg)
        summary[name] = result
        if result is False:
            all_pass = False

    # Summary
    print(f"\n{'='*60}")
    print("  SUMMARY")
    print(f"{'='*60}")
    for name, result in summary.items():
        status = "PASS" if result else ("SKIP" if result is None else "FAIL")
        print(f"  {name:25s} {status}")

    # Cleanup on full pass
    if all_pass:
        print("\nAll tests passed! Cleaning up validation_tmp/")
        shutil.rmtree(ROOT / "validation_tmp", ignore_errors=True)

    return 0 if all_pass else 1


if __name__ == "__main__":
    sys.exit(main())

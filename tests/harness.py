"""Shared test harness for the modular UMAT (standalone gfortran, no ABAQUS).

Reuses generate.py as the single source of truth for the source-file list and
PROPS/NSTATEV packing, so tests never re-encode the layout.

Provides:
  - build_umat_source(tmpdir)      -> writes umat.f90 (+ aba_param.inc)
  - compile_driver(tmpdir, driver) -> compiles driver + umat.f90 -> executable
  - run_exe(exe)                   -> (stdout, ok)
  - props_for(example_name)        -> (props, nstatev) from generate.EXAMPLES
  - F_* deformation-gradient helpers and fortran formatting helpers
"""
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import generate  # noqa: E402  (single source of truth for layout)

FFLAGS = ["-O2", "-ffree-form", "-ffpe-summary=none"]


# --- source assembly -------------------------------------------------------

def build_umat_source(tmpdir: Path) -> Path:
    """Concatenate the src/ modules into tmpdir/umat.f90 (same order as generate.py)."""
    out = tmpdir / "umat.f90"
    parts = []
    for src in generate.SOURCE_FILES:
        parts.append((generate.SRC_DIR / src).read_text())
    out.write_text("\n\n".join(parts) + "\n")
    (tmpdir / "aba_param.inc").write_text((generate.SRC_DIR / "aba_param.inc").read_text())
    return out


def compile_driver(tmpdir: Path, driver_src: str, name: str = "drv") -> Path:
    """Write a Fortran driver, compile it with umat.f90, return the executable path."""
    umat = build_umat_source(tmpdir)
    drv = tmpdir / f"{name}.f90"
    drv.write_text(driver_src)
    exe = tmpdir / name
    cmd = ["gfortran", *FFLAGS, "-o", str(exe), str(umat), str(drv)]
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=str(tmpdir))
    if r.returncode != 0:
        raise RuntimeError(f"compile failed:\n{r.stderr}")
    return exe


def run_exe(exe: Path) -> str:
    r = subprocess.run([str(exe)], capture_output=True, text=True, cwd=str(exe.parent))
    if r.returncode != 0:
        raise RuntimeError(f"run failed (rc={r.returncode}):\n{r.stdout}\n{r.stderr}")
    return r.stdout


# --- layout (delegated to generate.py) -------------------------------------

def props_for(example_name: str):
    cfg = generate.EXAMPLES[example_name]
    return generate.build_props(cfg), generate.compute_nstatev(cfg)


def props_from_cfg(cfg: dict):
    return generate.build_props(cfg), generate.compute_nstatev(cfg)


# --- single Cauchy-stress evaluation (shared by T2/T3) ----------------------

_STRESS_DRIVER = r"""
program stress_driver
  implicit none
  integer, parameter :: ntens=6, ndi=3, nshr=3, nprops={NP}, nstatev={NS}
  integer, parameter :: noel=1, npt=1
  double precision :: stress(ntens), statev(nstatev), ddsdde(ntens,ntens)
  double precision :: ddsddt(ntens), drplde(ntens), stran(ntens), dstran(ntens)
  double precision :: time(2), predef(1), dpred(1), props(nprops), coords(3), drot(3,3)
  double precision :: dfgrd0(3,3), dfgrd1(3,3)
  double precision :: sse,spd,scd,rpl,drpldt,dtime,temp,dtemp,pnewdt,celent
  character(len=8) :: cmname
  integer :: layer,kspt,kstep,kinc, ii
  double precision, parameter :: zero=0.d0, one=1.d0
  stress=zero; statev=zero; ddsdde=zero; stran=zero; dstran=zero; time=zero
  drot=zero; dfgrd0=zero; predef=zero; dpred=zero; coords=zero
  temp=zero; dtemp=zero; pnewdt=one; celent=one
  rpl=zero; drpldt=zero; ddsddt=zero; drplde=zero
  layer=1; kspt=1; kinc=1; kstep=1; cmname='UMAT'
  dtime={DTIME}d0
  dfgrd0(1,1)=one; dfgrd0(2,2)=one; dfgrd0(3,3)=one
{PROPS}
{FBASE}
  call uexternaldb(0,0,time,zero,0,0)
  call umat(stress, statev, ddsdde, sse, spd, scd, rpl, ddsddt, &
            drplde, drpldt, stran, dstran, time, dtime, temp, dtemp, &
            predef, dpred, cmname, ndi, nshr, ntens, nstatev, props, &
            nprops, coords, drot, pnewdt, celent, dfgrd0, dfgrd1, &
            noel, npt, layer, kspt, kstep, kinc)
  write(*,'(A)') 'STRESS'
  write(*,'(6ES24.14)') (stress(ii), ii=1,6)
  write(*,'(A)') 'DET'
  write(*,'(ES24.14)') statev(1)
  write(*,'(A)') 'SSE'
  write(*,'(ES24.14)') sse
contains
end program stress_driver

subroutine getoutdir(outdir, lenoutdir)
  implicit none
  character(len=256), intent(out) :: outdir
  integer, intent(out) :: lenoutdir
  outdir = '.'; lenoutdir = 1
end subroutine getoutdir
"""

import tempfile  # noqa: E402


def stress_at(props, nstatev, F, dtime=1e-2, setup_files=None):
    """Return (sigma[6], det, sse) for a single UMAT call at deformation gradient F.

    setup_files: optional {name: content} written into the run dir before exec
    (e.g. a sphere_int<N>c.inp quadrature file for RW network types).
    """
    src = (_STRESS_DRIVER
           .replace("{NP}", str(len(props)))
           .replace("{NS}", str(nstatev))
           .replace("{DTIME}d0", fdbl(dtime))
           .replace("{PROPS}", fmt_props(props))
           .replace("{FBASE}", fmt_F(F)))
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        for name, content in (setup_files or {}).items():
            (tmp / name).write_text(content)
        exe = compile_driver(tmp, src, "stress")
        out = run_exe(exe)
    lines = [ln for ln in out.splitlines() if ln.strip()]
    si = lines.index("STRESS")
    di = lines.index("DET")
    ei = lines.index("SSE")
    sigma = [float(x) for x in lines[si + 1].split()]
    det = float(lines[di + 1])
    sse = float(lines[ei + 1])
    return sigma, det, sse


# --- fortran formatting -----------------------------------------------------

def fdbl(v) -> str:
    """Format a Python float as a valid Fortran double-precision literal."""
    s = repr(float(v))
    if "e" in s or "E" in s:
        return s.replace("e", "d").replace("E", "d")
    if "inf" in s or "nan" in s:
        raise ValueError(f"non-finite literal: {v}")
    return s + "d0"


def fmt_props(props) -> str:
    return "\n".join(f"  props({i}) = {fdbl(v)}" for i, v in enumerate(props, 1))


def fmt_F(F, var="dfgrd1") -> str:
    lines = []
    for i in range(3):
        for j in range(3):
            lines.append(f"  {var}({i+1},{j+1}) = {fdbl(F[i][j])}")
    return "\n".join(lines)


# --- deformation gradients (row-major lists) --------------------------------

def F_identity():
    return [[1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0]]


def F_uniaxial(lam):
    s = lam ** -0.5
    return [[lam, 0, 0], [0, s, 0], [0, 0, s]]


def F_biaxial(lam):
    return [[lam, 0, 0], [0, lam, 0], [0, 0, 1.0 / (lam * lam)]]


def F_simple_shear(g):
    return [[1.0, g, 0], [0, 1.0, 0], [0, 0, 1.0]]


def F_general(lam=1.15, g=0.1):
    # off-axis: stretch + shear + lateral, isochoric-ish but generic
    return [[lam, g, 0.05],
            [0.0, 1.05, g],
            [0.02, 0.0, 1.0 / (lam * 1.05)]]


# --- sphere quadrature for RW network types (UEXTERNALDB input file) ---------

def fibonacci_sphere(n):
    """n quasi-uniform unit vectors on the sphere (Fibonacci spiral)."""
    import math
    ga = math.pi * (3.0 - math.sqrt(5.0))
    pts = []
    for i in range(n):
        z = 1.0 - 2.0 * (i + 0.5) / n
        r = math.sqrt(max(0.0, 1.0 - z * z))
        th = ga * i
        pts.append((r * math.cos(th), r * math.sin(th), z))
    return pts


def quadrature_file(n=60):
    """Return (filename, content) for a sphere_int<N>c.inp UEXTERNALDB reads.

    Lines are 'x y z weight'; weights are 1/n (normalized orientation average).
    UEXTERNALDB tries 720,300,60,21 in order, so n must be one of those.
    """
    assert n in (720, 300, 60, 21)
    w = 1.0 / n
    lines = [f"{x:.16e} {y:.16e} {z:.16e} {w:.16e}" for (x, y, z) in fibonacci_sphere(n)]
    return f"sphere_int{n}c.inp", "\n".join(lines) + "\n"

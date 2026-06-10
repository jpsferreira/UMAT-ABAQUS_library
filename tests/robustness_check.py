"""T8 - Robustness guards (no-NaN, graceful cutback on element inversion).

Verifies the UMAT degrades gracefully instead of poisoning the Newton iteration:
  E3: an inverted/degenerate element (det F <= 0) must return pnewdt < 1 (request
      a smaller increment) with finite stress, NOT NaN/Inf.
  control: a normal deformation must leave pnewdt unchanged (== 1) and finite.

Run: .venv/bin/python tests/robustness_check.py
"""
import math
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import harness as H  # noqa: E402

DRIVER = r"""
program robust_driver
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
  temp=zero; dtemp=zero; celent=one; rpl=zero; drpldt=zero; ddsddt=zero; drplde=zero
  layer=1; kspt=1; kinc=1; kstep=1; cmname='UMAT'
  dtime=1.d-2
  pnewdt=one                      ! ABAQUS passes a large value; 1.0 = "no change asked"
  dfgrd0(1,1)=one; dfgrd0(2,2)=one; dfgrd0(3,3)=one
{PROPS}
{FBASE}
  call uexternaldb(0,0,time,zero,0,0)
  call umat(stress, statev, ddsdde, sse, spd, scd, rpl, ddsddt, &
            drplde, drpldt, stran, dstran, time, dtime, temp, dtemp, &
            predef, dpred, cmname, ndi, nshr, ntens, nstatev, props, &
            nprops, coords, drot, pnewdt, celent, dfgrd0, dfgrd1, &
            noel, npt, layer, kspt, kstep, kinc)
  write(*,'(A)') 'PNEWDT'
  write(*,'(ES24.14)') pnewdt
  write(*,'(A)') 'STRESS'
  write(*,'(6ES24.14)') (stress(ii), ii=1,6)
contains
end program robust_driver

subroutine getoutdir(outdir, lenoutdir)
  implicit none
  character(len=256), intent(out) :: outdir
  integer, intent(out) :: lenoutdir
  outdir = '.'; lenoutdir = 1
end subroutine getoutdir
"""


def run(props, ns, F):
    src = (DRIVER.replace("{NP}", str(len(props))).replace("{NS}", str(ns))
           .replace("{PROPS}", H.fmt_props(props)).replace("{FBASE}", H.fmt_F(F)))
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        exe = H.compile_driver(tmp, src, "robust")
        out = H.run_exe(exe)
    lines = [ln for ln in out.splitlines() if ln.strip()]
    pnewdt = float(lines[lines.index("PNEWDT") + 1])
    stress = [float(x) for x in lines[lines.index("STRESS") + 1].split()]
    return pnewdt, stress


def _finite(xs):
    return all(math.isfinite(x) for x in xs)


def main(*_):
    print("T8 robustness (no-NaN, graceful cutback)\n")
    props, ns = H.props_for("neo_hooke")
    rc = 0

    # Control: a normal deformation leaves pnewdt == 1 with finite stress.
    pn, s = run(props, ns, H.F_uniaxial(1.2))
    if abs(pn - 1.0) < 1e-12 and _finite(s):
        print(f"[PASS] control: normal F -> pnewdt={pn:.3f}, finite stress")
    else:
        rc = 1
        print(f"[FAIL] control: pnewdt={pn}, finite={_finite(s)}")

    # E3: inverted element (det F < 0) -> cutback requested, finite (no NaN).
    F_inv = [[-0.5, 0, 0], [0, 1.0, 0], [0, 0, 1.0]]   # det = -0.5
    pn, s = run(props, ns, F_inv)
    if pn < 1.0 and _finite(s):
        print(f"[PASS] E3: inverted element (det<0) -> pnewdt={pn:.3f} (<1, cutback), finite stress")
    else:
        rc = 1
        print(f"[FAIL] E3: det<0 -> pnewdt={pn} (expected <1), finite={_finite(s)} (stress={s})")

    return rc


if __name__ == "__main__":
    sys.exit(main())

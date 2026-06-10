"""T4 - State-evolution tests (multi-step, carrying STATEV across increments).

The other suites do single-step-from-reference calls; these exercise the
time-integration / history that those cannot:
  - viscoelastic stress relaxation: hold a stretch; the Maxwell overstress must
    decay monotonically toward the elastic equilibrium.
  - damage ratchet: load / unload / reload; the damage variable must be
    non-decreasing and frozen during unload (irreversible damage).

Run: .venv/bin/python tests/evolution_check.py
"""
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import harness as H  # noqa: E402

DRIVER = r"""
program evolve
  implicit none
  integer, parameter :: ntens=6, ndi=3, nshr=3, nprops={NP}, nstatev={NS}
  integer, parameter :: noel=1, npt=1, nsteps={NSTEPS}
  double precision :: stress(ntens), statev(nstatev), ddsdde(ntens,ntens)
  double precision :: ddsddt(ntens), drplde(ntens), stran(ntens), dstran(ntens)
  double precision :: time(2), predef(1), dpred(1), props(nprops), coords(3), drot(3,3)
  double precision :: dfgrd0(3,3), dfgrd1(3,3)
  double precision :: sse,spd,scd,rpl,drpldt,dtime,temp,dtemp,pnewdt,celent
  character(len=8) :: cmname
  integer :: layer,kspt,kstep,kinc, ii, istep
  double precision :: lams(nsteps), lam
  double precision, parameter :: zero=0.d0, one=1.d0
  lams = (/ {LAMS} /)
  stress=zero; statev=zero; ddsdde=zero; stran=zero; dstran=zero; time=zero
  drot=zero; dfgrd0=zero; predef=zero; dpred=zero; coords=zero
  temp=zero; dtemp=zero; pnewdt=one; celent=one; rpl=zero; drpldt=zero
  ddsddt=zero; drplde=zero; layer=1; kspt=1; kinc=1; kstep=1; cmname='UMAT'
  dtime={DTIME}d0
  dfgrd0(1,1)=one; dfgrd0(2,2)=one; dfgrd0(3,3)=one
{PROPS}
  call uexternaldb(0,0,time,zero,0,0)
  do istep=1,nsteps
    lam = lams(istep)
    dfgrd1=zero
    dfgrd1(1,1)=lam; dfgrd1(2,2)=one/sqrt(lam); dfgrd1(3,3)=one/sqrt(lam)
    call umat(stress, statev, ddsdde, sse, spd, scd, rpl, ddsddt, &
              drplde, drpldt, stran, dstran, time, dtime, temp, dtemp, &
              predef, dpred, cmname, ndi, nshr, ntens, nstatev, props, &
              nprops, coords, drot, pnewdt, celent, dfgrd0, dfgrd1, &
              noel, npt, layer, kspt, kstep, kinc)
    write(*,'(I5,3ES22.12)') istep, lam, stress(1), statev(2)
    time(1)=time(1)+dtime
  end do
contains
end program evolve

subroutine getoutdir(outdir, lenoutdir)
  implicit none
  character(len=256), intent(out) :: outdir
  integer, intent(out) :: lenoutdir
  outdir = '.'; lenoutdir = 1
end subroutine getoutdir
"""


def run_evolution(props, ns, lams, dtime):
    src = (DRIVER.replace("{NP}", str(len(props))).replace("{NS}", str(ns))
           .replace("{NSTEPS}", str(len(lams)))
           .replace("{LAMS}", ", ".join(H.fdbl(x) for x in lams))
           .replace("{DTIME}d0", H.fdbl(dtime))
           .replace("{PROPS}", H.fmt_props(props)))
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        exe = H.compile_driver(tmp, src, "evolve")
        out = H.run_exe(exe)
    rows = []
    for ln in out.splitlines():
        p = ln.split()
        if len(p) == 4:
            rows.append((int(p[0]), float(p[1]), float(p[2]), float(p[3])))
    return rows  # (step, lambda, S11, damage)


def check_visco_relaxation():
    """Hold a stretch; viscoelastic overstress must relax toward the elastic value."""
    pv, nv = H.props_for("neo_hooke_visco")          # NH + 1 Maxwell branch
    pe, ne = H.props_for("neo_hooke")                # elastic equilibrium
    lam = 1.3
    rows = run_evolution(pv, nv, [lam] * 30, dtime=0.2)   # jump to 1.3, hold 30 steps
    s = [r[2] for r in rows]
    s_elastic = run_evolution(pe, ne, [lam], dtime=0.2)[0][2]
    fails = []
    if not (s[0] > s[-1]):
        fails.append(f"no relaxation: S11 first={s[0]:.4e} last={s[-1]:.4e}")
    # monotone non-increasing during hold
    if any(s[i + 1] > s[i] + 1e-6 * abs(s[0]) for i in range(len(s) - 1)):
        fails.append("S11 not monotonically decreasing during hold")
    # converges toward the elastic equilibrium
    if abs(s[-1] - s_elastic) > 0.05 * abs(s_elastic):
        fails.append(f"did not relax to elastic: last={s[-1]:.4e} elastic={s_elastic:.4e}")
    return fails, s[0], s[-1], s_elastic


def check_damage_ratchet():
    """Load / unload / reload; damage must be non-decreasing and frozen on unload."""
    import copy
    import generate
    cfg = copy.deepcopy(generate.EXAMPLES["neo_hooke_damage"])
    cfg["damage_params"] = [5.0, 3.0]   # psi_half=3 so damage activates by lam~1.4
    p, ns = H.props_from_cfg(cfg)
    # triangle: load 1.0->1.5, unload 1.5->1.0, reload 1.0->1.5
    up = [1.0 + 0.5 * i / 20 for i in range(21)]
    down = up[::-1][1:]
    rows = run_evolution(p, ns, up + down + up[1:], dtime=0.01)
    d = [r[3] for r in rows]
    fails = []
    # non-decreasing throughout
    if any(d[i + 1] < d[i] - 1e-12 for i in range(len(d) - 1)):
        fails.append("damage decreased somewhere (must be irreversible)")
    # damage must actually activate
    if d[-1] <= 1e-6:
        fails.append(f"damage never activated (max d={max(d):.2e}); tune psi_half")
    # frozen during unload: damage at end of unload == damage at peak load
    d_peak = d[20]
    d_unload_end = d[20 + len(down)]
    if abs(d_unload_end - d_peak) > 1e-9:
        fails.append(f"damage changed during unload: peak={d_peak:.4e} end={d_unload_end:.4e}")
    return fails, max(d)


def main(*_):
    print("T4 state evolution (visco relaxation, damage ratchet)\n")
    rc = 0

    fails, s0, sl, se = check_visco_relaxation()
    if fails:
        rc = 1
        print("[FAIL] viscoelastic relaxation")
        for f in fails:
            print(f"        {f}")
    else:
        print(f"[PASS] viscoelastic relaxation (hold: S11 {s0:.3e} -> {sl:.3e} -> elastic {se:.3e})")

    fails, dmax = check_damage_ratchet()
    if fails:
        rc = 1
        print("[FAIL] damage ratchet")
        for f in fails:
            print(f"        {f}")
    else:
        print(f"[PASS] damage ratchet (monotone, frozen on unload; max d={dmax:.3f})")

    return rc


if __name__ == "__main__":
    sys.exit(main())

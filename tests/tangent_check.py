"""T1 - Numerical consistency check of the UMAT material Jacobian (DDSDDE).

For a deformation gradient F, compares the analytic DDSDDE returned by the UMAT
against a central-difference tangent built with the standard Miehe symmetric
perturbation of F (the perturbation that reproduces ABAQUS's Jaumann-rate
material Jacobian):

    dF^(kl) = (eps/2) [ e_k (x) e_l . F + e_l (x) e_k . F ]
    C_num[:,kl] = ( sigma(F + dF) - sigma(F - dF) ) / (2 eps)   (Voigt columns)

Calibration: a verified-correct model (Neo-Hooke) MUST pass at tight tolerance.
Buggy configs (e.g. HGO with dispersion, active damage) are expected to FAIL
until the corresponding Theme-B/C fix lands.

Run:
    .venv/bin/python tests/tangent_check.py            # full battery
    .venv/bin/python tests/tangent_check.py neo_hooke  # one model, verbose
"""
import copy
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import harness as H  # noqa: E402
import generate  # noqa: E402

# Voigt (k,l) for columns 1..6 = [11,22,33,12,13,23]
VK = [1, 2, 3, 1, 1, 2]
VL = [1, 2, 3, 2, 3, 3]

DRIVER = r"""
program tangent_driver
  implicit none
  integer, parameter :: ntens=6, ndi=3, nshr=3, nprops={NP}, nstatev={NS}
  integer, parameter :: noel=1, npt=1

  double precision :: stress(ntens), statev(nstatev), ddsdde(ntens,ntens)
  double precision :: ddsddt(ntens), drplde(ntens), stran(ntens), dstran(ntens)
  double precision :: time(2), predef(1), dpred(1), props(nprops), coords(3), drot(3,3)
  double precision :: dfgrd0(3,3), dfgrd1(3,3)
  double precision :: sse,spd,scd,rpl,drpldt,dtime,temp,dtemp,pnewdt,celent
  character(len=8) :: cmname
  integer :: layer,kspt,kstep,kinc

  double precision :: F(3,3), Fp(3,3), s0(ntens), sp(ntens), sm(ntens)
  double precision :: cnum(ntens,ntens), dd0(ntens,ntens)
  double precision :: eps, j0, jp, jm
  integer :: vk(6), vl(6), kk, ll, jj, ii, mm
  double precision, parameter :: zero=0.d0, one=1.d0

  vk = (/1,2,3,1,1,2/)
  vl = (/1,2,3,2,3,3/)
  eps = {EPS}d0
  dtime = {DTIME}d0

{PROPS}

{FBASE}
  F = dfgrd1

  call uexternaldb(0,0,time,zero,0,0)

  ! --- base point: analytic tangent (J0 from statev(1)=det) ---
  call eval(F, s0, dd0, j0)

  ! --- central-difference tangent, column by column ---
  ! ABAQUS material Jacobian C_ijkl = (1/J) d(J*sigma_ij)/d(eps_kl), so we
  ! difference the KIRCHHOFF stress tau = J*sigma and divide by J at the base
  ! point. (Differencing Cauchy sigma alone omits the sigma_ij*delta_kl term.)
  do jj = 1, 6
    kk = vk(jj); ll = vl(jj)
    ! +dF
    Fp = F
    do ii = 1,3
      do mm = 1,3
        Fp(ii,mm) = F(ii,mm) + (eps/2.d0)*( delta(ii,kk)*F(ll,mm) + delta(ii,ll)*F(kk,mm) )
      end do
    end do
    call eval(Fp, sp, ddsdde, jp)
    ! -dF
    Fp = F
    do ii = 1,3
      do mm = 1,3
        Fp(ii,mm) = F(ii,mm) - (eps/2.d0)*( delta(ii,kk)*F(ll,mm) + delta(ii,ll)*F(kk,mm) )
      end do
    end do
    call eval(Fp, sm, ddsdde, jm)
    do ii = 1, 6
      cnum(ii,jj) = (jp*sp(ii) - jm*sm(ii)) / (2.d0*eps) / j0
    end do
  end do

  ! --- emit both matrices ---
  write(*,'(A)') 'ANALYTIC'
  do ii = 1, 6
    write(*,'(6ES24.14)') (dd0(ii,jj), jj=1,6)
  end do
  write(*,'(A)') 'NUMERICAL'
  do ii = 1, 6
    write(*,'(6ES24.14)') (cnum(ii,jj), jj=1,6)
  end do

contains

  pure double precision function delta(i,j)
    integer, intent(in) :: i,j
    delta = 0.d0
    if (i==j) delta = 1.d0
  end function delta

  subroutine eval(Fin, stress_out, dd_out, det_out)
    double precision, intent(in)  :: Fin(3,3)
    double precision, intent(out) :: stress_out(ntens), dd_out(ntens,ntens), det_out
    ! Reset to reference state so the tangent is the algorithmic d(stress)/d(strain)
    ! holding the incoming (reference) state fixed.
    stress = zero; statev = zero; ddsdde = zero
    stran = zero; dstran = zero; time = zero
    drot = zero; dfgrd0 = zero; predef = zero; dpred = zero; coords = zero
    temp = zero; dtemp = zero; pnewdt = one; celent = one
    rpl = zero; drpldt = zero; ddsddt = zero; drplde = zero
    layer = 1; kspt = 1; kinc = 1; kstep = 1; cmname = 'UMAT'
    dfgrd0(1,1)=one; dfgrd0(2,2)=one; dfgrd0(3,3)=one
    dfgrd1 = Fin
    call umat(stress, statev, ddsdde, sse, spd, scd, rpl, ddsddt, &
              drplde, drpldt, stran, dstran, time, dtime, temp, dtemp, &
              predef, dpred, cmname, ndi, nshr, ntens, nstatev, props, &
              nprops, coords, drot, pnewdt, celent, dfgrd0, dfgrd1, &
              noel, npt, layer, kspt, kstep, kinc)
    stress_out = stress
    dd_out = ddsdde
    det_out = statev(1)   ! umat writes statev(1) = det(F)
  end subroutine eval

end program tangent_driver

! Stub for ABAQUS-provided routine (standalone builds only)
subroutine getoutdir(outdir, lenoutdir)
  implicit none
  character(len=256), intent(out) :: outdir
  integer, intent(out) :: lenoutdir
  outdir = '.'
  lenoutdir = 1
end subroutine getoutdir
"""


def _parse_matrix(lines, start):
    m = []
    for r in range(start, start + 6):
        m.append([float(x) for x in lines[r].split()])
    return m


def tangent_error(props, nstatev, F, eps=1e-6, dtime=1e-2, setup_files=None):
    """Compile+run the tangent driver; return (analytic, numerical, max_abs_diff, scale)."""
    src = (DRIVER
           .replace("{NP}", str(len(props)))
           .replace("{NS}", str(nstatev))
           .replace("{EPS}d0", H.fdbl(eps))
           .replace("{DTIME}d0", H.fdbl(dtime))
           .replace("{PROPS}", H.fmt_props(props))
           .replace("{FBASE}", H.fmt_F(F)))
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        for name, content in (setup_files or {}).items():
            (tmp / name).write_text(content)
        exe = H.compile_driver(tmp, src, "tangent")
        out = H.run_exe(exe)
    lines = [ln for ln in out.splitlines() if ln.strip()]
    ai = lines.index("ANALYTIC")
    ni = lines.index("NUMERICAL")
    analytic = _parse_matrix(lines, ai + 1)
    numerical = _parse_matrix(lines, ni + 1)
    maxdiff = 0.0
    scale = 0.0
    for i in range(6):
        for j in range(6):
            maxdiff = max(maxdiff, abs(analytic[i][j] - numerical[i][j]))
            scale = max(scale, abs(analytic[i][j]))
    return analytic, numerical, maxdiff, scale


STATES = [
    ("uniaxial λ=1.2", H.F_uniaxial(1.2)),
    ("biaxial λ=1.2", H.F_biaxial(1.2)),
    ("simple shear γ=0.2", H.F_simple_shear(0.2)),
    ("general", H.F_general()),
]

# Fiber models gate stress/stiffness on the e1=0 (I4=1) tension threshold, which
# is non-smooth. A pure simple shear with the fiber along e1 lands EXACTLY on that
# kink (I4=C11=1), where any finite-difference tangent is ill-posed (see todo B4).
# Tangent verification must use states with the fiber clearly engaged (I4>1).
FIBER_STATES = [
    ("uniaxial λ=1.3 (fiber)", H.F_uniaxial(1.3)),
    ("biaxial λ=1.2", H.F_biaxial(1.2)),
    ("general (stretch+shear)", H.F_general()),
]

def _cfg(name, **over):
    c = copy.deepcopy(generate.EXAMPLES[name])
    c.update(over)
    return c


# HGO with zero dispersion: no I1-I4 coupling -> B1 cannot manifest (control).
_hgo_k0 = _cfg("humphrey_hgo", aniso_params=[100.0, 10.0, 0.0, 1.0, 0.0, 0.0])
# Damage tuned to be ON the steep part of the sigmoid at the test stretch
# (psi_half ~ W(lambda=1.2) ~ 1), so dmg_diff is O(1) and C1 shows up.
_dmg_active = _cfg("neo_hooke_damage", damage_params=[5.0, 1.0])

# Each entry: (label, config-dict, tag, states). tag 'good' gates the suite;
# 'suspect' is expected to FAIL and documents a known bug until its fix lands.
MODELS = [
    ("neo_hooke", generate.EXAMPLES["neo_hooke"], "good", STATES),
    ("mooney_rivlin", generate.EXAMPLES["mooney_rivlin"], "good", STATES),
    ("ogden_3term", generate.EXAMPLES["ogden_3term"], "good", STATES),
    ("neo_hooke_visco", generate.EXAMPLES["neo_hooke_visco"], "good", STATES),
    ("hgo_kappa0 (control)", _hgo_k0, "good", FIBER_STATES),
    ("humphrey_hgo κ=0.226", generate.EXAMPLES["humphrey_hgo"], "good", FIBER_STATES),  # B1 fixed
    ("neo_hooke_damage (inactive)", generate.EXAMPLES["neo_hooke_damage"], "good", STATES),
    ("damage_active", _dmg_active, "good", STATES),  # C1 fixed (tangent only; stress unchanged)
]

REL_TOL = 1e-4


def _fmt(m):
    return "\n".join("  " + "  ".join(f"{v: .6e}" for v in row) for row in m)


def main(argv):
    if len(argv) > 1:
        name = argv[1]
        props, ns = H.props_for(name)
        print(f"== {name}  (nprops={len(props)}, nstatev={ns}) ==")
        for label, F in STATES:
            a, n, md, sc = tangent_error(props, ns, F)
            rel = md / sc if sc > 0 else md
            verdict = "PASS" if rel < REL_TOL else "FAIL"
            print(f"  {label:22s} maxdiff={md:.3e} scale={sc:.3e} rel={rel:.3e}  {verdict}")
            if argv[-1] == "-v":
                print(" ANALYTIC\n" + _fmt(a) + "\n NUMERICAL\n" + _fmt(n))
        return 0

    # battery
    print(f"T1 numerical-tangent battery (rel tol {REL_TOL:g})\n")
    rc = 0
    for name, cfg, tag, states in MODELS:
        props, ns = H.props_from_cfg(cfg)
        worst = 0.0
        details = []
        for label, F in states:
            try:
                _, _, md, sc = tangent_error(props, ns, F)
                rel = md / sc if sc > 0 else md
            except Exception as e:  # compile/run failure is itself a finding
                details.append(f"{label}: ERROR {e}")
                rel = float("inf")
            worst = max(worst, rel)
            details.append(f"{label}: rel={rel:.2e}")
        ok = worst < REL_TOL
        flag = "PASS" if ok else "FAIL"
        note = "" if (ok == (tag == "good")) else "  <-- unexpected!"
        if tag == "suspect" and not ok:
            note = "  (expected: documents a known bug)"
        print(f"[{flag}] {name:20s} worst_rel={worst:.2e}  ({tag}){note}")
        for d in details:
            print(f"        {d}")
        # only 'good' models gate the suite
        if tag == "good" and not ok:
            rc = 1
    return rc


if __name__ == "__main__":
    sys.exit(main(sys.argv))

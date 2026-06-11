"""T9 - Numerical consistency check of the UEL assembly tangent (material_uel).

The displacement-based UEL assembles Kuu = int G^T A G dv, which requires the
push-forward of the FIRST elasticity tensor (no minor symmetry):

    a_ijkl = (1/J) F_jJ F_lL dP_iJ/dF_kL = c_ijkl + delta_ik sigma_jl

This is NOT the Jaumann-corrected modulus the UMAT path returns; getting this
term wrong leaves the residual intact but degrades Newton convergence — exactly
the latent bug the material_core/material_uel refactor fixed.

Check: central finite differences of the first Piola-Kirchhoff stress
P = J sigma F^{-T} under direct (unsymmetric) perturbation of single F entries,
pushed forward at the base point, compared against material_uel's analytic
tangent over all 81 components.

Calibration (lessons.md): a verified-good model (Neo-Hooke) must pass at tight
tolerance before any failure elsewhere is trusted. Fiber models are evaluated
with the fiber engaged (I4>1) to stay off the I4=1 activation kink.

Run:
    .venv/bin/python tests/uel_tangent_check.py            # full battery
    .venv/bin/python tests/uel_tangent_check.py neo_hooke  # one model, verbose
"""
import copy
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import harness as H  # noqa: E402
import generate  # noqa: E402

DRIVER = r"""
program uel_tangent_driver
  use mod_material, only: material_uel
  use mod_tensor,   only: matinv3d
  implicit none
  integer, parameter :: nprops={NP}, nstatev={NS}
  integer, parameter :: noel=1, npt=1

  double precision :: statev(nstatev), props(nprops)
  double precision :: time(2), predef(1), dpred(1)
  double precision :: dtime, pnewdt, sse
  double precision :: dfgrd1(3,3)
  double precision :: F(3,3), Fp(3,3)
  double precision :: sig0(3,3), a0(3,3,3,3)
  double precision :: sigp(3,3), sigm(3,3), adum(3,3,3,3)
  double precision :: Pp(3,3), Pm(3,3)
  double precision :: Anum(3,3,3,3), a_num(3,3,3,3)
  double precision :: j0, jp, jm, eps
  integer :: i, jj, kk, ll, mJ, mL
  double precision, parameter :: zero=0.d0, one=1.d0

  eps = {EPS}d0
  dtime = {DTIME}d0

{PROPS}

{FBASE}
  F = dfgrd1

  time = zero
  call uexternaldb(0,0,time,zero,0,0)

  ! --- base point: analytic UEL tangent ---
  call eval(F, sig0, a0, j0)

  ! --- central-difference dP/dF, one F entry at a time (unsymmetric) ---
  Anum = zero
  do kk = 1, 3
    do ll = 1, 3
      Fp = F
      Fp(kk,ll) = F(kk,ll) + eps
      call eval(Fp, sigp, adum, jp)
      call piola(Fp, sigp, jp, Pp)
      Fp = F
      Fp(kk,ll) = F(kk,ll) - eps
      call eval(Fp, sigm, adum, jm)
      call piola(Fp, sigm, jm, Pm)
      do i = 1, 3
        do mJ = 1, 3
          Anum(i,mJ,kk,ll) = (Pp(i,mJ) - Pm(i,mJ)) / (2.d0*eps)
        end do
      end do
    end do
  end do

  ! --- push forward: a_ijkl = (1/J) F_jJ F_lL A_iJkL ---
  a_num = zero
  do i = 1, 3
    do jj = 1, 3
      do kk = 1, 3
        do ll = 1, 3
          do mJ = 1, 3
            do mL = 1, 3
              a_num(i,jj,kk,ll) = a_num(i,jj,kk,ll) &
                + F(jj,mJ) * F(ll,mL) * Anum(i,mJ,kk,mL) / j0
            end do
          end do
        end do
      end do
    end do
  end do

  ! --- emit both tensors (81 components each, fixed order) ---
  write(*,'(A)') 'ANALYTIC'
  write(*,'(3ES24.14)') ((((a0(i,jj,kk,ll), ll=1,3), kk=1,3), jj=1,3), i=1,3)
  write(*,'(A)') 'NUMERICAL'
  write(*,'(3ES24.14)') ((((a_num(i,jj,kk,ll), ll=1,3), kk=1,3), jj=1,3), i=1,3)

contains

  !> P = J sigma F^{-T}  (first Piola-Kirchhoff from Cauchy)
  subroutine piola(Fin, sig, jdet, P)
    double precision, intent(in)  :: Fin(3,3), sig(3,3), jdet
    double precision, intent(out) :: P(3,3)
    double precision :: Finv(3,3)
    integer :: ii, jl, mm
    call matinv3d(Fin, Finv)
    P = zero
    do ii = 1, 3
      do jl = 1, 3
        do mm = 1, 3
          P(ii,jl) = P(ii,jl) + jdet * sig(ii,mm) * Finv(jl,mm)
        end do
      end do
    end do
  end subroutine piola

  !> One material_uel call from a fresh reference state (algorithmic tangent).
  subroutine eval(Fin, sig_out, a_out, det_out)
    double precision, intent(in)  :: Fin(3,3)
    double precision, intent(out) :: sig_out(3,3), a_out(3,3,3,3), det_out
    double precision :: Floc(3,3)
    statev = zero; time = zero; predef = zero; dpred = zero
    pnewdt = one
    Floc = Fin
    call material_uel(sig_out, a_out, sse, statev, nstatev, props, nprops, &
                      Floc, time, dtime, predef, dpred, pnewdt, &
                      noel, npt, 1, 1)
    det_out = statev(1)   ! material_core writes statev(1) = det(F)
  end subroutine eval

end program uel_tangent_driver

! Stub for ABAQUS-provided routine (standalone builds only)
subroutine getoutdir(outdir, lenoutdir)
  implicit none
  character(len=256), intent(out) :: outdir
  integer, intent(out) :: lenoutdir
  outdir = '.'
  lenoutdir = 1
end subroutine getoutdir
"""


def _parse_tensor(lines, start):
    """81 values, 3 per line, 27 lines."""
    vals = []
    for r in range(start, start + 27):
        vals.extend(float(x) for x in lines[r].split())
    return vals


def uel_tangent_error(props, nstatev, F, eps=1e-6, dtime=1e-2, setup_files=None):
    """Compile+run the driver; return (analytic, numerical, max_abs_diff, scale)."""
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
        exe = H.compile_driver(tmp, src, "uel_tangent")
        out = H.run_exe(exe)
    lines = [ln for ln in out.splitlines() if ln.strip()]
    ai = lines.index("ANALYTIC")
    ni = lines.index("NUMERICAL")
    analytic = _parse_tensor(lines, ai + 1)
    numerical = _parse_tensor(lines, ni + 1)
    maxdiff = max(abs(a - n) for a, n in zip(analytic, numerical))
    scale = max(abs(a) for a in analytic)
    return analytic, numerical, maxdiff, scale


STATES = [
    ("uniaxial λ=1.2", H.F_uniaxial(1.2)),
    ("biaxial λ=1.2", H.F_biaxial(1.2)),
    ("simple shear γ=0.2", H.F_simple_shear(0.2)),
    ("general", H.F_general()),
]

# Fiber models gate on the I4=1 tension kink — evaluate with the fiber engaged
# (see tangent_check.py and lessons.md).
FIBER_STATES = [
    ("uniaxial λ=1.3 (fiber)", H.F_uniaxial(1.3)),
    ("biaxial λ=1.2", H.F_biaxial(1.2)),
    ("general (stretch+shear)", H.F_general()),
]


def _cfg(name, **over):
    c = copy.deepcopy(generate.EXAMPLES[name])
    c.update(over)
    return c


_hgo_k0 = _cfg("humphrey_hgo", aniso_params=[100.0, 10.0, 0.0, 1.0, 0.0, 0.0])
_dmg_active = _cfg("neo_hooke_damage", damage_params=[5.0, 1.0])

MODELS = [
    ("neo_hooke", generate.EXAMPLES["neo_hooke"], "good", STATES),
    ("mooney_rivlin", generate.EXAMPLES["mooney_rivlin"], "good", STATES),
    ("ogden_3term", generate.EXAMPLES["ogden_3term"], "good", STATES),
    ("neo_hooke_visco", generate.EXAMPLES["neo_hooke_visco"], "good", STATES),
    ("hgo_kappa0 (control)", _hgo_k0, "good", FIBER_STATES),
    ("humphrey_hgo κ=0.226", generate.EXAMPLES["humphrey_hgo"], "good", FIBER_STATES),
    ("neo_hooke_damage (inactive)", generate.EXAMPLES["neo_hooke_damage"], "good", STATES),
    ("damage_active", _dmg_active, "good", STATES),
]

REL_TOL = 1e-4


def main(argv=None):
    argv = argv or [sys.argv[0]]
    if len(argv) > 1:
        name = argv[1]
        props, ns = H.props_for(name)
        print(f"== {name}  (nprops={len(props)}, nstatev={ns}) ==")
        for label, F in STATES:
            _, _, md, sc = uel_tangent_error(props, ns, F)
            rel = md / sc if sc > 0 else md
            verdict = "PASS" if rel < REL_TOL else "FAIL"
            print(f"  {label:22s} maxdiff={md:.3e} scale={sc:.3e} rel={rel:.3e}  {verdict}")
        return 0

    print(f"T9 UEL-tangent battery: a_ijkl = c + δ⊗σ vs FD of dP/dF (rel tol {REL_TOL:g})\n")
    rc = 0
    for name, cfg, tag, states in MODELS:
        props, ns = H.props_from_cfg(cfg)
        worst = 0.0
        details = []
        for label, F in states:
            try:
                _, _, md, sc = uel_tangent_error(props, ns, F)
                rel = md / sc if sc > 0 else md
            except Exception as e:
                details.append(f"{label}: ERROR {e}")
                rel = float("inf")
            worst = max(worst, rel)
            details.append(f"{label}: rel={rel:.2e}")
        ok = worst < REL_TOL
        flag = "PASS" if ok else "FAIL"
        note = "" if (ok == (tag == "good")) else "  <-- unexpected!"
        print(f"[{flag}] {name:20s} worst_rel={worst:.2e}  ({tag}){note}")
        for d in details:
            print(f"        {d}")
        if tag == "good" and not ok:
            rc = 1
    return rc


if __name__ == "__main__":
    sys.exit(main(sys.argv))

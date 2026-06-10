"""T10 - End-to-end checks of the generated UEL (element + material library).

Two checks on a single U3D8 element (unit cube, affine nodal displacements
u_a = (F - I) X_a, one increment from the virgin state):

(a) Traction oracle — under affine loading the element's deformation gradient
    is spatially constant and equals the prescribed F (F-bar inert: Fc = F),
    so the internal force on the x+ face must equal the UMAT answer exactly:
        f_pred = sigma_umat . (J F^{-T} e1)        (Nanson, reference area 1)
    This ties the WHOLE element chain (kinematics, material call, B-matrix,
    Gauss loop) to the independently verified UMAT stress.

(b) Stiffness consistency — AMATRX vs central finite differences of RHS over
    all 24 DOFs (state reset per call):  K_fd = -dRHS/du.
    With fbar=off this must match tightly (consistent tangent through the
    assembly — validates the geometric term a = c + δ⊗σ end to end).
    With fbar=on it validates the F-bar tangent (Qmat) at an affine base
    state, where de Souza Neto's q-matrix linearization is exact.

Run:
    .venv/bin/python tests/uel_check.py
"""
import sys
import tempfile
import subprocess
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import harness as H  # noqa: E402
import generate  # noqa: E402

DRIVER = r"""
program uel_fd_driver
  implicit none
  integer, parameter :: nnode = 8, ndofel = 24, mlvarx = 24, nrhs = 1
  integer, parameter :: nprops = {NP}, nsdv = {NS}, nintp = {NINT}
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

  double precision :: Ftar(3,3), u0(ndofel), up(ndofel)
  double precision :: rhs0(ndofel), rhsp(ndofel), rhsm(ndofel)
  double precision :: k0(ndofel,ndofel), kfd(ndofel,ndofel), fface(3)
  double precision :: h
  integer :: i, j, n, k, face_nodes(4)
  double precision, parameter :: zero = 0.0d0, one = 1.0d0

{PROPS}
  jprops(1) = nsdv
  jprops(2) = nsdv
  lflags = 0; lflags(1) = 1; lflags(2) = 1
  dtime = {DTIME}d0
  h = {H}d0
  v = zero; a = zero; params = zero
  adlmag = zero; ddlmag = zero; predef = zero; jdltyp = 0
  period = zero
  face_nodes = (/2, 3, 6, 7/)

  ! Unit cube, element-local node ordering (z=0 face: 1-4, z=1 face: 5-8)
  coords(:,1) = (/zero, zero, zero/)
  coords(:,2) = (/one,  zero, zero/)
  coords(:,3) = (/one,  one,  zero/)
  coords(:,4) = (/zero, one,  zero/)
  coords(:,5) = (/zero, zero, one /)
  coords(:,6) = (/one,  zero, one /)
  coords(:,7) = (/one,  one,  one /)
  coords(:,8) = (/zero, one,  one /)

{FBASE}

  time = zero
  call uexternaldb(0, 0, time, zero, 0, 0)

  ! Affine nodal displacements for the base state
  do n = 1, nnode
    do k = 1, 3
      u0(3*(n-1)+k) = sum(Ftar(k,:)*coords(:,n)) - coords(k,n)
    end do
  end do

  ! --- base call: RHS, AMATRX, face force ---
  call eval(u0, rhs0, k0)
  fface = zero
  do i = 1, 4
    n = face_nodes(i)
    do k = 1, 3
      fface(k) = fface(k) - rhs0(3*(n-1)+k)
    end do
  end do

  ! --- central FD of RHS over all DOFs (fresh state each call) ---
  do j = 1, ndofel
    up = u0; up(j) = u0(j) + h
    call eval(up, rhsp, amatrx)
    up = u0; up(j) = u0(j) - h
    call eval(up, rhsm, amatrx)
    kfd(:,j) = -(rhsp - rhsm) / (2.0d0*h)
  end do

  write(*,'(A)') 'FFACE'
  write(*,'(3ES24.14)') fface
  write(*,'(A)') 'RHS'
  write(*,'(3ES24.14)') (rhs0(i), i=1,ndofel)
  write(*,'(A)') 'AMATRX'
  do i = 1, ndofel
    write(*,'(4ES24.14)') (k0(i,j), j=1,ndofel)
  end do
  write(*,'(A)') 'KFD'
  do i = 1, ndofel
    write(*,'(4ES24.14)') (kfd(i,j), j=1,ndofel)
  end do

contains

  !> One UEL call from the virgin state (single increment, kinc=1).
  subroutine eval(uin, rhs_out, k_out)
    double precision, intent(in)  :: uin(ndofel)
    double precision, intent(out) :: rhs_out(ndofel), k_out(ndofel,ndofel)
    svars = zero; time = zero; pnewdt = one
    rhs = zero; amatrx = zero; energy = zero
    u = uin
    du(:,1) = uin   ! single increment from u=0
    call uel(rhs, amatrx, svars, energy, ndofel, nrhs, nsvars, &
             props, nprops, coords, mcrd, nnode, u, du, v, a, jtype, &
             time, dtime, 1, 1, jelem, params, ndload, jdltyp, adlmag, &
             predef, npredf, lflags, mlvarx, ddlmag, mdload, pnewdt, &
             jprops, njprop, period)
    rhs_out = rhs(:,1)
    k_out = amatrx
  end subroutine eval

end program uel_fd_driver

! Stub for ABAQUS-provided routine (standalone builds only)
subroutine getoutdir(outdir, lenoutdir)
  implicit none
  character(len=256), intent(out) :: outdir
  integer, intent(out) :: lenoutdir
  outdir = '.'
  lenoutdir = 1
end subroutine getoutdir
"""

FFLAGS = ["-O2", "-ffree-form", "-ffpe-summary=none"]


def _compile_and_run(tmp, driver_src, elem):
    (tmp / "uel.f90").write_text(generate.uel_source(elem))
    (tmp / "aba_param.inc").write_text((generate.SRC_DIR / "aba_param.inc").read_text())
    drv = tmp / "uel_drv.f90"
    drv.write_text(driver_src)
    exe = tmp / "uel_drv"
    r = subprocess.run(["gfortran", *FFLAGS, "-o", str(exe), str(tmp / "uel.f90"), str(drv)],
                       capture_output=True, text=True, cwd=str(tmp))
    if r.returncode != 0:
        raise RuntimeError(f"compile failed:\n{r.stderr}")
    return H.run_exe(exe)


def _parse_vec(lines, start, n, per_line):
    vals = []
    rows = (n + per_line - 1) // per_line
    for r in range(start, start + rows):
        vals.extend(float(x) for x in lines[r].split())
    return vals[:n]


def _parse_mat(lines, start, n, per_line):
    rows_per = (n + per_line - 1) // per_line
    m = []
    for i in range(n):
        m.append(_parse_vec(lines, start + i * rows_per, n, per_line))
    return m


def uel_run(props, nstatev, F, elem, h=1e-6, dtime=1e-2, setup_files=None):
    """Return (fface[3], rhs[24], K[24][24], Kfd[24][24]) for one increment to F."""
    src = (DRIVER
           .replace("{NP}", str(len(props)))
           .replace("{NS}", str(nstatev))
           .replace("{NINT}", str(elem["nint"]))
           .replace("{DTIME}d0", H.fdbl(dtime))
           .replace("{H}d0", H.fdbl(h))
           .replace("{PROPS}", H.fmt_props(props))
           .replace("{FBASE}", H.fmt_F(F, var="Ftar")))
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        for name, content in (setup_files or {}).items():
            (tmp / name).write_text(content)
        out = _compile_and_run(tmp, src, elem)
    lines = [ln for ln in out.splitlines() if ln.strip()]
    fface = _parse_vec(lines, lines.index("FFACE") + 1, 3, 3)
    rhs = _parse_vec(lines, lines.index("RHS") + 1, 24, 3)
    K = _parse_mat(lines, lines.index("AMATRX") + 1, 24, 4)
    Kfd = _parse_mat(lines, lines.index("KFD") + 1, 24, 4)
    return fface, rhs, K, Kfd


# --- small 3x3 helpers (no numpy in the test env) ---------------------------

def det3(F):
    return (F[0][0]*(F[1][1]*F[2][2] - F[1][2]*F[2][1])
            - F[0][1]*(F[1][0]*F[2][2] - F[1][2]*F[2][0])
            + F[0][2]*(F[1][0]*F[2][1] - F[1][1]*F[2][0]))


def inv3(F):
    d = det3(F)
    c = [[(F[(i+1) % 3][(j+1) % 3]*F[(i+2) % 3][(j+2) % 3]
           - F[(i+1) % 3][(j+2) % 3]*F[(i+2) % 3][(j+1) % 3]) / d
          for j in range(3)] for i in range(3)]
    # c[i][j] = cofactor -> inverse = transpose of cofactor matrix / det
    return [[c[j][i] for j in range(3)] for i in range(3)]


def predicted_face_force(sigma_voigt, F):
    """sigma . (J F^{-T} e1): Nanson area vector of the X=1 face (ref area 1)."""
    s = sigma_voigt
    sig = [[s[0], s[3], s[4]],
           [s[3], s[1], s[5]],
           [s[4], s[5], s[2]]]
    J = det3(F)
    Finv = inv3(F)
    avec = [J * Finv[0][i] for i in range(3)]   # (F^{-T})_{i1} = Finv[0][i]
    return [sum(sig[i][j] * avec[j] for j in range(3)) for i in range(3)]


# --- checks ------------------------------------------------------------------

ORACLE_TOL = 1e-8
FD_TOL = 1e-6

CASES = [
    # (label, example name, fbar, states)
    ("neo_hooke fbar=on", "neo_hooke", True,
     [("general", H.F_general()), ("uniaxial λ=1.3", H.F_uniaxial(1.3))]),
    ("neo_hooke fbar=off", "neo_hooke", False,
     [("general", H.F_general())]),
    ("humphrey_hgo fbar=on", "humphrey_hgo", True,
     [("uniaxial λ=1.3 (fiber)", H.F_uniaxial(1.3))]),
    ("neo_hooke_visco fbar=on", "neo_hooke_visco", True,
     [("general", H.F_general())]),
]


def main(argv=None):
    print(f"T10 UEL end-to-end (oracle tol {ORACLE_TOL:g}, FD-stiffness tol {FD_TOL:g})\n")
    rc = 0
    for label, name, fbar, states in CASES:
        cfg = generate.EXAMPLES[name]
        props, ns = H.props_from_cfg(cfg)
        elem = dict(generate.DEFAULT_ELEMENT)
        elem["fbar"] = fbar
        for slabel, F in states:
            try:
                fface, rhs, K, Kfd = uel_run(props, ns, F, elem)
                sigv, _, _ = H.stress_at(props, ns, F)
                fpred = predicted_face_force(sigv, F)

                fscale = max(abs(v) for v in fpred)
                ferr = max(abs(a - b) for a, b in zip(fface, fpred)) / fscale

                kscale = max(abs(v) for row in K for v in row)
                kerr = max(abs(a - b) for ra, rb in zip(K, Kfd)
                           for a, b in zip(ra, rb)) / kscale

                ok = ferr < ORACLE_TOL and kerr < FD_TOL
                flag = "PASS" if ok else "FAIL"
                print(f"[{flag}] {label:24s} {slabel:22s} "
                      f"oracle_rel={ferr:.2e}  K_fd_rel={kerr:.2e}")
                if not ok:
                    rc = 1
            except Exception as e:
                print(f"[FAIL] {label:24s} {slabel:22s} ERROR: {e}")
                rc = 1
    return rc


if __name__ == "__main__":
    sys.exit(main(sys.argv))

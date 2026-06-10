"""T3 - Closed-form analytic oracles for the UMAT Cauchy stress.

These pin absolute correctness (not just self-consistency) against textbook
results (Holzapfel), so a wrong-but-self-consistent implementation can't hide.

  1. Volumetric: vol-only material at pure dilatation F = a*I gives
     sigma = K/2 (J - 1/J) * I, deviator = 0.
  2. Neo-Hooke isochoric uniaxial (F = diag(λ, λ^-1/2, λ^-1/2), J=1):
     sigma11 = (4/3) C10 (λ^2 - 1/λ),  sigma22 = sigma33 = -(2/3) C10 (λ^2 - 1/λ).
  3. One-term Ogden with α=2 is identically Neo-Hooke with C10 = μ/2.

Run: .venv/bin/python tests/oracle_check.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import harness as H  # noqa: E402

TOL = 1e-6  # relative


def _rel(a, b, scale):
    return abs(a - b) / scale if scale else abs(a - b)


def check_volumetric():
    K = 1000.0
    cfg = {"name": "vol", "kbulk": K, "iso_type": 0, "iso_params": [],
           "aniso_type": 0, "n_fiber_fam": 0, "aniso_params": [],
           "network_type": 0, "network_params": [], "damage_type": 0,
           "damage_params": [], "n_visco": 0, "visco_params": []}
    props, ns = H.props_from_cfg(cfg)
    fails = []
    for a in (1.05, 1.1, 0.95):
        F = [[a, 0, 0], [0, a, 0], [0, 0, a]]
        sigma, det, _ = H.stress_at(props, ns, F)
        J = a ** 3
        p = 0.5 * K * (J - 1.0 / J)
        scale = max(abs(p), 1.0)
        for i, comp in enumerate(("11", "22", "33")):
            if _rel(sigma[i], p, scale) > TOL:
                fails.append(f"a={a} σ{comp}={sigma[i]:.6e} expected {p:.6e}")
        for i, comp in zip((3, 4, 5), ("12", "13", "23")):
            if abs(sigma[i]) > TOL * scale:
                fails.append(f"a={a} σ{comp}={sigma[i]:.3e} expected 0")
    return fails


def check_neo_hooke_uniaxial():
    K, C10 = 1.0e6, 10.0  # high bulk -> J stays ~1; F is isochoric so pv≈0 anyway
    cfg = {"name": "nh", "kbulk": K, "iso_type": 1, "iso_params": [C10],
           "aniso_type": 0, "n_fiber_fam": 0, "aniso_params": [],
           "network_type": 0, "network_params": [], "damage_type": 0,
           "damage_params": [], "n_visco": 0, "visco_params": []}
    props, ns = H.props_from_cfg(cfg)
    fails = []
    for lam in (1.1, 1.3, 0.9):
        F = H.F_uniaxial(lam)
        sigma, det, _ = H.stress_at(props, ns, F)
        d = lam * lam - 1.0 / lam
        s11 = (4.0 / 3.0) * C10 * d
        s22 = -(2.0 / 3.0) * C10 * d
        scale = max(abs(s11), 1.0)
        if _rel(sigma[0], s11, scale) > TOL:
            fails.append(f"λ={lam} σ11={sigma[0]:.6e} expected {s11:.6e}")
        if _rel(sigma[1], s22, scale) > TOL or _rel(sigma[2], s22, scale) > TOL:
            fails.append(f"λ={lam} σ22={sigma[1]:.6e}/σ33={sigma[2]:.6e} expected {s22:.6e}")
    return fails


def check_ogden_equals_neohooke():
    K, mu = 1.0e6, 20.0  # Ogden(N=1, μ, α=2) == NH(C10=μ/2)
    nh = {"name": "nh", "kbulk": K, "iso_type": 1, "iso_params": [mu / 2.0],
          "aniso_type": 0, "n_fiber_fam": 0, "aniso_params": [],
          "network_type": 0, "network_params": [], "damage_type": 0,
          "damage_params": [], "n_visco": 0, "visco_params": []}
    og = dict(nh, iso_type=3, iso_params=[1, mu, 2.0])
    pnh, ns1 = H.props_from_cfg(nh)
    pog, ns2 = H.props_from_cfg(og)
    fails = []
    for label, F in (("uniaxial", H.F_uniaxial(1.3)),
                     ("biaxial", H.F_biaxial(1.2)),
                     ("shear", H.F_simple_shear(0.2))):
        snh, _, _ = H.stress_at(pnh, ns1, F)
        sog, _, _ = H.stress_at(pog, ns2, F)
        scale = max(max(abs(x) for x in snh), 1.0)
        for i in range(6):
            if _rel(snh[i], sog[i], scale) > TOL:
                fails.append(f"{label} comp{i}: NH={snh[i]:.6e} Ogden={sog[i]:.6e}")
    return fails


def check_ogden_isochoric_j_independence():
    """Deviatoric KIRCHHOFF stress (tau = J*sigma) depends only on Cbar, not J.

    The isochoric Cauchy stress scales as 1/J (correct); the J-invariant measure
    is the Kirchhoff stress. Near-isotropic Cbar triggers the Ogden all-equal
    (neig==3) branch, where the C-vs-Cbar invariant bug (B2) scales the response
    with J. Same Cbar (=> same Fbar, same isochoric Kirchhoff stress) at J=1 and
    J!=1 must give the same deviatoric Kirchhoff stress.
    """
    props, ns = H.props_for("ogden_3term")
    s = 1.0 + 1e-6  # all three Cbar eigenvalues within neig tol (1e-5) -> neig==3
    base = [[s, 0, 0], [0, 1.0, 0], [0, 0, 1.0 / s]]

    def devK(sig, J):
        t = [J * x for x in sig]
        tr = (t[0] + t[1] + t[2]) / 3.0
        return [t[0] - tr, t[1] - tr, t[2] - tr, t[3], t[4], t[5]]

    sa, Ja, _ = H.stress_at(props, ns, base)
    da = devK(sa, Ja)
    scale = max(max(abs(x) for x in da), 1e-30)
    fails = []
    for c in (1.3, 0.8):
        Fb = [[c * v for v in row] for row in base]
        sb, Jb, _ = H.stress_at(props, ns, Fb)
        db = devK(sb, Jb)
        rel = max(abs(da[i] - db[i]) for i in range(6)) / scale
        if rel > 1e-6:
            fails.append(f"J={c**3:.2f}: dev Kirchhoff stress drifts rel={rel:.2e} from J=1 (must be ~0)")
    return fails


def _fiber_only(aniso_type, aniso_params):
    return {"name": "fib", "kbulk": 500.0, "iso_type": 0, "iso_params": [],
            "aniso_type": aniso_type, "n_fiber_fam": 1, "aniso_params": aniso_params,
            "network_type": 0, "network_params": [], "damage_type": 0,
            "damage_params": [], "n_visco": 0, "visco_params": []}


def check_fiber_compression_energy():
    """Fiber strain energy must vanish when the fiber is in compression.

    The models already return zero fiber stress/stiffness for a slack fiber
    (e1 <= 0); the strain energy must switch off consistently, otherwise it
    reports energy with no conjugate stress and pollutes the damage criterion
    (sef = sseiso + sseaniso).
    """
    # Fiber along e1, compressed: F = diag(0.8, 1/sqrt(0.8), 1/sqrt(0.8)), J=1
    # so I4bar = 0.64 < 1 (slack fiber). With J=1 the volumetric energy is 0.
    a = 0.8
    F = [[a, 0, 0], [0, 1.0 / a**0.5, 0], [0, 0, 1.0 / a**0.5]]
    vol = _fiber_only(0, [])
    _, _, sse_vol = H.stress_at(*H.props_from_cfg(vol), F)
    fails = []
    cases = [("Humphrey-fiber", _fiber_only(2, [100.0, 10.0, 1.0, 0.0, 0.0])),
             ("HGO κ=0", _fiber_only(1, [100.0, 10.0, 0.0, 1.0, 0.0, 0.0]))]
    for name, cfg in cases:
        _, _, sse = H.stress_at(*H.props_from_cfg(cfg), F)
        if abs(sse - sse_vol) > 1e-8:
            fails.append(f"{name}: spurious fiber energy {sse - sse_vol:.3e} in compression (must be 0)")
    return fails


def main():
    print(f"T3 closed-form oracle checks (rel tol {TOL:g})\n")
    rc = 0
    for name, fn in (("volumetric pressure p(J)", check_volumetric),
                     ("Neo-Hooke uniaxial closed form", check_neo_hooke_uniaxial),
                     ("Ogden(α=2) == Neo-Hooke(C10=μ/2)", check_ogden_equals_neohooke),
                     ("Ogden isochoric stress is J-independent (neig==3)", check_ogden_isochoric_j_independence),
                     ("Fiber strain energy vanishes in compression", check_fiber_compression_energy)):
        fails = fn()
        if fails:
            rc = 1
            print(f"[FAIL] {name}")
            for f in fails[:6]:
                print(f"        {f}")
        else:
            print(f"[PASS] {name}")
    return rc


if __name__ == "__main__":
    sys.exit(main())

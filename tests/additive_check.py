"""T2 - Additive-contribution decomposition: volumetric + isotropic + anisotropic.

The framework assembles sigma = svol + PE:(sfic_iso + sfic_aniso) (+ network).
Running each block in isolation (others disabled via the PROPS header) and
sharing the same KBULK, the total must equal the sum of contributions:

    sigma(vol+iso+aniso) == sigma(iso) + sigma(aniso) - sigma(vol)

because the volumetric part is counted once in each of the iso-only and
aniso-only runs. This holds by construction UNLESS a block leaks into another,
so it is a regression guard against future cross-contamination between blocks.

Note: this invariant holds for ANY kappa (verified below for kappa=0 and
kappa=0.226). The B1 dispersion-coupling defect is a *dead write* (diso is
mutated but never read once the iso fictitious tensors are built), so it drops
the coupling consistently from both the aniso-only and combined runs and does
NOT break additivity. B1 is a missing-physics bug, caught by the T1 tangent
check (and would be caught by an absolute HGO-dispersion stress oracle), not here.

Also checks contribution isolation:
  - volumetric stress is always spherical (deviator = 0),
  - isotropic and anisotropic contributions vanish at F = I.

Run: .venv/bin/python tests/additive_check.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import harness as H  # noqa: E402

TOL = 1e-7   # relative; additivity holds to machine precision when clean
KBULK = 500.0
C10 = 10.0
FIBER = [1.0, 0.0, 0.0]          # along e1
HGO = [100.0, 10.0]             # k1, k2


def _base():
    return {"name": "x", "kbulk": KBULK, "iso_type": 0, "iso_params": [],
            "aniso_type": 0, "n_fiber_fam": 0, "aniso_params": [],
            "network_type": 0, "network_params": [], "damage_type": 0,
            "damage_params": [], "n_visco": 0, "visco_params": []}


def _cfgs(kappa):
    vol = _base()
    iso = dict(_base(), iso_type=1, iso_params=[C10])
    ani = dict(_base(), aniso_type=1, n_fiber_fam=1,
               aniso_params=HGO + [kappa] + FIBER)
    allc = dict(iso, aniso_type=1, n_fiber_fam=1,
                aniso_params=HGO + [kappa] + FIBER)
    return vol, iso, ani, allc


def _sig(cfg, F):
    p, ns = H.props_from_cfg(cfg)
    s, _, _ = H.stress_at(p, ns, F)
    return s


def check_additivity(kappa):
    """Return worst relative additivity residual over a few engaged states."""
    vol, iso, ani, allc = _cfgs(kappa)
    worst = 0.0
    for label, F in (("uniaxial λ=1.3", H.F_uniaxial(1.3)),
                     ("general", H.F_general())):
        sv = _sig(vol, F)
        si = _sig(iso, F)
        sa = _sig(ani, F)
        st = _sig(allc, F)
        scale = max(max(abs(x) for x in st), 1.0)
        for k in range(6):
            resid = abs(st[k] - (si[k] + sa[k] - sv[k]))
            worst = max(worst, resid / scale)
    return worst


def check_isolation():
    fails = []
    vol, iso, ani, _ = _cfgs(0.226)
    # volumetric stress is spherical at a generic (non-dilatational) F
    F = H.F_general()
    sv = _sig(vol, F)
    if max(abs(sv[0] - sv[1]), abs(sv[0] - sv[2])) > 1e-6 * max(abs(sv[0]), 1.0):
        fails.append(f"vol not spherical: {sv[:3]}")
    if max(abs(sv[3]), abs(sv[4]), abs(sv[5])) > 1e-6 * max(abs(sv[0]), 1.0):
        fails.append(f"vol has shear: {sv[3:]}")
    # iso and aniso contributions vanish at F = I
    for name, cfg in (("iso", iso), ("aniso", ani)):
        s = _sig(cfg, H.F_identity())
        if max(abs(x) for x in s) > 1e-8:
            fails.append(f"{name} stress at F=I not zero: {s}")
    return fails


def main():
    print(f"T2 additive decomposition (vol + iso + aniso, rel tol {TOL:g})\n")
    rc = 0

    fails = check_isolation()
    if fails:
        rc = 1
        print("[FAIL] contribution isolation")
        for f in fails:
            print(f"        {f}")
    else:
        print("[PASS] contribution isolation (vol spherical; iso/aniso vanish at F=I)")

    # Additivity is a structural invariant that must hold for every kappa.
    for kappa in (0.0, 0.226):
        w = check_additivity(kappa)
        ok = w < TOL
        if not ok:
            rc = 1
        print(f"[{'PASS' if ok else 'FAIL'}] additivity (vol+iso+aniso), HGO κ={kappa:<5} worst_rel={w:.2e}")

    return rc


if __name__ == "__main__":
    sys.exit(main())

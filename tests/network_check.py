"""T6 - Filament-network battery (standalone, no ABAQUS).

Covers the network models that have ZERO reference coverage in validate.py:
  - smoke: every network type (1..6) compiles, runs, and returns finite,
    non-NaN Cauchy stress + symmetric DDSDDE;
  - RW quadrature loading: types 1/2 produce a non-zero network with a
    sphere_int<N>c.inp present, and SILENTLY ZERO without it (documents E6);
  - tangent consistency: analytic DDSDDE vs finite-difference per type (GATING;
    all network types are consistent to < 1e-4 after the C5 fix).

Run: .venv/bin/python tests/network_check.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import harness as H          # noqa: E402
import tangent_check as TC   # noqa: E402
import generate              # noqa: E402

QUAD = dict([H.quadrature_file(60)])  # {sphere_int60c.inp: content}


def _net_cfg(network_type, network_params):
    return {"name": "net", "kbulk": 500.0, "iso_type": 0, "iso_params": [],
            "aniso_type": 0, "n_fiber_fam": 0, "aniso_params": [],
            "network_type": network_type, "network_params": network_params,
            "damage_type": 0, "damage_params": [], "n_visco": 0, "visco_params": []}


# RW configs (not in generate.EXAMPLES). Counts pinned by T5/umat_builder.
FIL = [1.0, 0.1, 0.01, 2.0, 0.1, 1.0, 0.0, 0.0]   # L,R0F,mu0,beta,B0,lambda0,R0C,ETAC
RW_AFFINE = _net_cfg(1, [0.5, 1.0e6, 2.0, 1.0, 1.0, 0.0, 0.0] + FIL)       # +pdir
RW_NONAFF = _net_cfg(2, [0.5, 1.0e6, 2.0, 1.0, 2.0, 1.0, 0.0, 0.0] + FIL)  # +PP +pdir

# (label, cfg, needs_quadrature)
MODELS = [
    ("affine-RW (1)", RW_AFFINE, True),
    ("nonaffine-RW (2)", RW_NONAFF, True),
    ("mixed-AI (3)", generate.EXAMPLES["mixed_network"], False),
    ("contractile-AI (4)", generate.EXAMPLES["contractile_network"], False),
    ("affine-AI (5)", generate.EXAMPLES["affine_network"], False),
    ("nonaffine-AI (6)", generate.EXAMPLES["nonaffine_network"], False),
]


def _finite(xs):
    return all(x == x and abs(x) < 1e30 for x in xs)  # x==x rejects NaN


def smoke():
    """Every type: finite, non-NaN stress + symmetric DDSDDE at a moderate stretch."""
    F = H.F_uniaxial(1.2)
    fails = []
    for label, cfg, needs_q in MODELS:
        p, ns = H.props_from_cfg(cfg)
        sf = QUAD if needs_q else None
        try:
            s, det, _ = H.stress_at(p, ns, F, setup_files=sf)
        except Exception as e:
            fails.append(f"{label}: ERROR {str(e)[:70]}")
            continue
        if not _finite(s):
            fails.append(f"{label}: non-finite stress {s}")
        # DDSDDE symmetry from the analytic tangent matrix
        a, _, _, _ = TC.tangent_error(p, ns, F, setup_files=sf)
        asym = max(abs(a[i][j] - a[j][i]) for i in range(6) for j in range(6))
        sc = max(max(abs(a[i][j]) for j in range(6)) for i in range(6)) or 1.0
        if asym / sc > 1e-6:
            fails.append(f"{label}: DDSDDE not symmetric (rel {asym/sc:.2e})")
    return fails


def rw_quadrature_loading():
    """RW affine: non-zero with quadrature present; zero without (documents E6)."""
    p, ns = H.props_from_cfg(RW_AFFINE)
    F = H.F_uniaxial(1.3)
    s_with, _, _ = H.stress_at(p, ns, F, setup_files=QUAD)
    s_without, _, _ = H.stress_at(p, ns, F, setup_files=None)
    fails = []
    if max(abs(x) for x in s_with) < 1.0:
        fails.append(f"RW affine with quadrature gave ~zero stress {s_with}")
    if max(abs(x) for x in s_without) > 1e-30:
        fails.append(f"RW affine WITHOUT quadrature gave non-zero stress {s_without} "
                     f"(expected the current silent-zero)")
    return fails, s_with, s_without


def b3_prefdir_sensitivity():
    """B3: non-affine RW must respond to the preferred direction.

    With b_orient>0, prefdir=e1 vs e2 must give a different stress. The old code
    hardwired angle=0 and ignored prefdir entirely, so this was identically 0.
    """
    import copy

    def run(pdir):
        cfg = copy.deepcopy(RW_NONAFF)
        cfg["network_params"][2] = 5.0       # b_orient
        cfg["network_params"][5:8] = pdir    # preferred direction
        p, ns = H.props_from_cfg(cfg)
        s, _, _ = H.stress_at(p, ns, H.F_uniaxial(1.3), setup_files=QUAD)
        return s

    e1 = run([1.0, 0.0, 0.0])
    e2 = run([0.0, 1.0, 0.0])
    rel = max(abs(e1[i] - e2[i]) for i in range(6)) / max(max(abs(x) for x in e1), 1.0)
    return rel


def tangent_consistency():
    """Report analytic-vs-FD tangent error per type (documented; does not gate)."""
    F = H.F_uniaxial(1.2)
    rows = []
    for label, cfg, needs_q in MODELS:
        p, ns = H.props_from_cfg(cfg)
        sf = QUAD if needs_q else None
        try:
            _, _, md, scl = TC.tangent_error(p, ns, F, setup_files=sf)
            rows.append((label, md / scl if scl else md))
        except Exception as e:
            rows.append((label, f"ERR {str(e)[:40]}"))
    return rows


def main(*_):
    print("T6 network battery\n")
    rc = 0

    fails = smoke()
    if fails:
        rc = 1
        print("[FAIL] smoke (finite stress + symmetric DDSDDE)")
        for f in fails:
            print(f"        {f}")
    else:
        print("[PASS] smoke: all 6 network types finite, non-NaN, symmetric DDSDDE")

    qf, s_with, s_without = rw_quadrature_loading()
    if qf:
        rc = 1
        print("[FAIL] RW quadrature loading")
        for f in qf:
            print(f"        {f}")
    else:
        print(f"[PASS] RW quadrature loading (|σ| with file={max(abs(x) for x in s_with):.2e}, "
              f"without file={max(abs(x) for x in s_without):.2e} -> silent zero, see E6)")

    b3 = b3_prefdir_sensitivity()
    if b3 > 0.1:
        print(f"[PASS] B3 non-affine RW responds to preferred direction (rel e1-vs-e2={b3:.2f})")
    else:
        rc = 1
        print(f"[FAIL] B3 non-affine RW ignores prefdir (rel e1-vs-e2={b3:.2e})")

    tc = tangent_consistency()
    tc_fail = [lab for lab, rel in tc if not (isinstance(rel, float) and rel < 1e-4)]
    if tc_fail:
        rc = 1
        print("[FAIL] network tangent consistency (analytic DDSDDE vs FD):")
    else:
        print("[PASS] network tangent consistency (analytic DDSDDE vs FD, all types < 1e-4):")
    for label, rel in tc:
        val = f"{rel:.2e}" if isinstance(rel, float) else str(rel)
        mark = "ok" if (isinstance(rel, float) and rel < 1e-4) else "INCONSISTENT"
        print(f"        {label:20s} rel={val}  {mark}")

    return rc


if __name__ == "__main__":
    sys.exit(main())

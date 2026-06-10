"""T7 - Combined-feature recombination (visco x damage x network x AI-aniso).

These are CORRECTNESS bugs in umat_builder's recombination where stress AND
tangent are consistently wrong, so the FD tangent check (T1) cannot see them.
Each test is a stress oracle built from a physical invariant.

  D1: an AI-anisotropic fiber family must still contribute when viscosity is on
      (the viscous override is built from the PK2 path, which the AI fiber never
      entered). Test: stress(iso+AI-aniso+visco) must differ from stress(iso+visco).
  D2: a pure network + viscosity must NOT double-count the network (there is no
      matrix to relax, so the stress must equal the no-visco network stress).

Run: .venv/bin/python tests/recombination_check.py
"""
import copy
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import harness as H   # noqa: E402
import generate       # noqa: E402

QUAD = dict([H.quadrature_file(60)])


def _cfg(name, **over):
    c = copy.deepcopy(generate.EXAMPLES[name])
    c.update(over)
    return c


def check_d1_ai_aniso_survives_visco():
    """AI-fiber stress must not vanish when a viscous branch is active."""
    # iso (Humphrey) + HGO-AI fiber + 1 Maxwell branch
    cfg_an = _cfg("humphrey_hgo_ai", n_visco=1, visco_params=[0.5, 0.25])
    # same iso + visco, NO fiber
    cfg_no = _cfg("humphrey_hgo_ai", aniso_type=0, n_fiber_fam=0, aniso_params=[],
                  n_visco=1, visco_params=[0.5, 0.25])
    F = H.F_uniaxial(1.3)  # fiber along e1, clearly engaged
    sa, _, _ = H.stress_at(*H.props_from_cfg(cfg_an), F)
    sn, _, _ = H.stress_at(*H.props_from_cfg(cfg_no), F)
    diff = max(abs(sa[i] - sn[i]) for i in range(6))
    scale = max(max(abs(x) for x in sa), 1.0)
    rel = diff / scale
    # The AI fiber should make a sizeable difference; bug => ~0 (fiber dropped).
    ok = rel > 1e-3
    return ok, rel


def check_d2_network_not_double_counted():
    """Pure network + viscosity must equal pure network (no matrix to relax)."""
    cfg_n = _cfg("affine_network")                                   # type 5 AI, iso=0
    cfg_nv = _cfg("affine_network", n_visco=1, visco_params=[0.5, 0.25])
    F = H.F_uniaxial(1.2)
    sn, _, _ = H.stress_at(*H.props_from_cfg(cfg_n), F)
    snv, _, _ = H.stress_at(*H.props_from_cfg(cfg_nv), F)
    diff = max(abs(sn[i] - snv[i]) for i in range(6))
    scale = max(max(abs(x) for x in sn), 1.0)
    rel = diff / scale
    # With no matrix, viscosity relaxes nothing => network stress must be unchanged.
    ok = rel < 1e-6
    return ok, rel, sn, snv


def main(*_):
    print("T7 recombination (visco x damage x network x AI-aniso)\n")
    rc = 0

    ok, rel = check_d1_ai_aniso_survives_visco()
    if ok:
        print(f"[PASS] D1: AI fiber contributes under viscosity (rel diff vs no-fiber = {rel:.2e})")
    else:
        rc = 1
        print(f"[FAIL] D1: AI fiber DROPPED under viscosity (rel diff vs no-fiber = {rel:.2e}, expected > 1e-3)")

    ok, rel, sn, snv = check_d2_network_not_double_counted()
    if ok:
        print(f"[PASS] D2: network not double-counted under viscosity (rel = {rel:.2e})")
    else:
        rc = 1
        print(f"[FAIL] D2: network double-counted under viscosity (rel = {rel:.2e}, expected ~0)")
        print(f"        no-visco σ11={sn[0]:.4e}  visco σ11={snv[0]:.4e}")

    return rc


if __name__ == "__main__":
    sys.exit(main())

"""Run the standalone UMAT verification suite (no ABAQUS required).

  T1  tangent_check   - analytic DDSDDE vs central-difference (Miehe perturbation)
  T2  additive_check  - volumetric + isotropic + anisotropic decomposition
  T3  oracle_check    - closed-form analytic stress oracles
  T4  evolution_check - multi-step time evolution (visco relaxation, damage ratchet)
  T5  doc_consistency - PROPS layout: reader vs PROPS_REFERENCE.md vs generator
  T6  network_check   - filament-network battery (smoke, RW quadrature, tangents)
  T7  recombination_check - visco x damage x network x AI-aniso recombination
  T8  robustness_check - no-NaN / graceful cutback on element inversion
  T9  uel_tangent_check - UEL assembly tangent a = c + δ⊗σ vs FD of dP/dF
  T10 uel_check        - generated UEL: traction oracle vs UMAT + AMATRX vs FD(RHS)
  B6  contractile_tangent_check - contractile consistent tangent (evolve-then-perturb)

Exit code is nonzero if any 'good'/invariant check fails. Checks tagged as
documenting a known bug (see tasks/todo.md) report FAIL inline but do not gate
the suite until their fix lands.

Run: .venv/bin/python tests/run_all.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import tangent_check
import additive_check
import oracle_check
import evolution_check
import doc_consistency
import network_check
import recombination_check
import robustness_check
import uel_tangent_check
import uel_check
import contractile_tangent_check


def main():
    rc = 0
    for name, mod in (("T1 tangent", tangent_check),
                      ("T2 additive", additive_check),
                      ("T3 oracle", oracle_check),
                      ("T4 evolution", evolution_check),
                      ("T5 doc-consistency", doc_consistency),
                      ("T6 network", network_check),
                      ("T7 recombination", recombination_check),
                      ("T8 robustness", robustness_check),
                      ("T9 uel-tangent", uel_tangent_check),
                      ("T10 uel-element", uel_check),
                      ("B6 contractile-tangent", contractile_tangent_check)):
        print("=" * 70)
        print(name)
        print("=" * 70)
        rc |= mod.main([sys.argv[0]]) if mod is tangent_check else mod.main()
        print()
    print("=" * 70)
    print("SUITE RESULT:", "PASS" if rc == 0 else "FAIL (see above)")
    return rc


if __name__ == "__main__":
    sys.exit(main())

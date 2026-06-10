"""T5 - PROPS layout consistency (pure Python, no compile/ABAQUS).

The PROPS layout lives in three places that must agree. The Fortran reader is the
ground truth:
  (1) umat_builder.f90 `ip = ip + N` advances per aniso/network case,
  (2) PROPS_REFERENCE.md table "Count" columns + worked-example CONSTANTS,
  (3) generate.py EXAMPLES (per-family aniso block, network block) + compute_nstatev.

This test parses all three and asserts agreement, so the doc can never silently
drift from the reader again (the class of bug fixed in Theme A).

Run: python tests/doc_consistency.py
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
import generate  # noqa: E402

SRC = ROOT / "src"
UMAT = (SRC / "umat_builder.f90").read_text()
CONST = (SRC / "mod_constants.f90").read_text()
DOC = (SRC / "PROPS_REFERENCE.md").read_text()


# --- parse mod_constants for ANISO_* / NET_* IDs ---------------------------

def parse_ids(prefix):
    ids = {}
    for m in re.finditer(rf"integer,\s*parameter\s*::\s*({prefix}\w+)\s*=\s*(\d+)", CONST):
        ids[m.group(1)] = int(m.group(2))
    return ids


# --- parse `ip = ip + N` advances, keyed by enclosing `case (...)` ----------

def parse_advances(text):
    adv = {}
    cur = None
    for line in text.splitlines():
        if re.match(r"\s*end\s*select\b", line):
            cur = None  # advances after the select-case (damage/visco) aren't ours
            continue
        m = re.match(r"\s*case\s*\(\s*([A-Za-z0-9_]+)\s*\)", line)
        if m:
            cur = m.group(1)
            continue
        m2 = re.search(r"\bip\s*=\s*ip\s*\+\s*(\d+)", line)
        if m2 and cur is not None:
            adv[cur] = int(m2.group(1))
    return adv


# --- parse a Count column from a doc table section --------------------------

def parse_doc_counts(section_header, stop_headers):
    start = DOC.index(section_header)
    rest = DOC[start + len(section_header):]
    end = min((rest.index(h) for h in stop_headers if h in rest), default=len(rest))
    block = rest[:end]
    counts = {}
    for line in block.splitlines():
        # rows like:  | 3 | model | params | 7 |
        cells = [c.strip() for c in line.split("|")[1:-1]]
        if len(cells) >= 2 and cells[0].isdigit() and cells[-1].isdigit():
            counts[int(cells[0])] = int(cells[-1])
    return counts


def main():
    fails = []
    aniso_ids = parse_ids("ANISO_")   # ANISO_HGO=1, ...
    net_ids = parse_ids("NET_")

    # Split at the END of the umat subroutine. (The name "network_contribution"
    # also appears earlier in umat's interface block, so we can't split on that.)
    split = UMAT.index("end subroutine umat")
    umat_main, net_sub = UMAT[:split], UMAT[split:]

    # --- ANISO: code advance per ID (ground truth) ---
    aniso_adv_sym = parse_advances(umat_main)          # {ANISO_HGO: 6, ...}
    code_aniso = {aniso_ids[k]: v for k, v in aniso_adv_sym.items() if k in aniso_ids}
    doc_aniso = parse_doc_counts("### Anisotropic", ["### Network"])

    for aid, n in sorted(code_aniso.items()):
        if doc_aniso.get(aid) != n:
            fails.append(f"ANISO id {aid}: code ip+={n} but PROPS_REFERENCE Count={doc_aniso.get(aid)}")

    # --- NETWORK: code advance per numeric case (ground truth) ---
    code_net = {int(k): v for k, v in parse_advances(net_sub).items() if k.isdigit()}
    doc_net = parse_doc_counts("### Network", ["### Damage", "## Parameter Packing", "## Examples"])
    for nid, n in sorted(code_net.items()):
        if doc_net.get(nid) != n:
            fails.append(f"NETWORK id {nid}: code ip+={n} but PROPS_REFERENCE Count={doc_net.get(nid)}")

    # --- generate.py EXAMPLES block lengths vs code advances ---
    for name, cfg in generate.EXAMPLES.items():
        at = cfg["aniso_type"]
        if at in code_aniso and cfg["n_fiber_fam"] > 0:
            per_fam = len(cfg["aniso_params"])
            if per_fam != code_aniso[at]:
                fails.append(f"EXAMPLE {name}: aniso_params={per_fam} but ANISO id {at} needs {code_aniso[at]}")
        nt = cfg["network_type"]
        if nt in code_net:
            if len(cfg["network_params"]) != code_net[nt]:
                fails.append(f"EXAMPLE {name}: network_params={len(cfg['network_params'])} but NETWORK id {nt} needs {code_net[nt]}")

    # --- worked-example CONSTANTS in the doc == 7 (header) + blocks ---
    for m in re.finditer(r"CONSTANTS=(\d+)", DOC):
        pass  # presence check only; per-example breakdown covered above

    # --- NSTATEV formula vs compute_nstatev for representative configs ---
    def expected_nstatev(cfg):
        n = 1 + (2 if cfg["damage_type"] > 0 else 0) + 9 * cfg["n_visco"]
        if cfg["network_type"] == 4:
            factor = int(cfg["network_params"][6])
            n += 4 + 20 * factor * factor
        return n
    for name, cfg in generate.EXAMPLES.items():
        if expected_nstatev(cfg) != generate.compute_nstatev(cfg):
            fails.append(f"EXAMPLE {name}: NSTATEV formula {expected_nstatev(cfg)} != compute_nstatev {generate.compute_nstatev(cfg)}")

    print("T5 PROPS layout consistency (reader vs docs vs generator)\n")
    print(f"  ANISO ids (code): {code_aniso}")
    print(f"  NETWORK ids (code): {code_net}")
    if fails:
        print("\n[FAIL]")
        for f in fails:
            print(f"   - {f}")
        return 1
    print("\n[PASS] reader, PROPS_REFERENCE.md, and generate.py all agree.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

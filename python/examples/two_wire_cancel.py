"""Generate the two-wire cancellation scenario JSON.

Run this script from anywhere; it writes to ``inputs/two_wire_cancel.json``
relative to the repository root and prints the path for convenience.
"""

from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
PYTHON_DIR = REPO_ROOT / "python"
if str(PYTHON_DIR) not in sys.path:
    sys.path.insert(0, str(PYTHON_DIR))

from scenario_api import Domain, Material, Scenario, UniformRegion, Wire


def main() -> None:
    dom = Domain(Lx=0.20, Ly=0.20, nx=201, ny=201)
    air = Material(name="air", mu_r=1.0)
    regions = [UniformRegion(material="air")]
    wires = [
        Wire(x=-0.03, y=0.0, radius=0.003, I=10.0),
        Wire(x=0.03, y=0.0, radius=0.003, I=-10.0),
    ]

    scenario = Scenario(domain=dom, materials=[air], regions=regions, sources=wires)

    output_path = REPO_ROOT / "inputs" / "two_wire_cancel.json"
    scenario.save_json(output_path)
    print(f"Wrote scenario to {output_path}")


if __name__ == "__main__":
    main()

"""Generate the two-wire scenario JSON used by regression tests.

Run this script from anywhere. By default it overwrites
``inputs/two_wire_cancel.json`` relative to the repository root and prints the
resolved path. Pass ``--output`` to write to an alternate location, which is
handy for CI smoke tests that should not touch the checked-in fixture.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
PYTHON_DIR = REPO_ROOT / "python"
if str(PYTHON_DIR) not in sys.path:
    sys.path.insert(0, str(PYTHON_DIR))

from scenario_api import (
    Domain,
    FieldMapOutput,
    LineProbeOutput,
    Material,
    Scenario,
    UniformRegion,
    Wire,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=REPO_ROOT / "inputs" / "two_wire_cancel.json",
        help=(
            "Destination JSON file. Defaults to the repository fixture "
            "(inputs/two_wire_cancel.json)."
        ),
    )
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    dom = Domain(Lx=0.20, Ly=0.20, nx=201, ny=201)
    air = Material(name="air", mu_r=1.0)
    regions = [UniformRegion(material="air")]
    wires = [
        Wire(x=-0.03, y=0.0, radius=0.003, I=10.0),
        Wire(x=0.03, y=0.0, radius=0.003, I=-10.0),
    ]

    outputs = [
        FieldMapOutput(id="domain_field", path="outputs/two_wire_field_map.csv"),
        LineProbeOutput(id="midline", axis="x", value=0.0, path="outputs/two_wire_midline.csv"),
    ]

    scenario = Scenario(
        domain=dom,
        materials=[air],
        regions=regions,
        sources=wires,
        outputs=outputs,
    )

    output_path = scenario.save_json(args.output)
    print(f"Wrote scenario to {output_path}")


if __name__ == "__main__":
    main()

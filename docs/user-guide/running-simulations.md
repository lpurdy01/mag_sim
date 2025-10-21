# Running Simulations

Follow this guide to execute bundled scenarios, manage timelines, and export outputs.

## 1. Generate a scenario

Use the provided Python helpers to emit JSON inputs. For example, the three-phase stator generator produces both CI-sized and high-resolution variants:

```bash
python3 python/gen_three_phase_stator.py --profile ci --out inputs/three_phase_stator_ci.json
```

Switch to the high-resolution preset via `--profile hires`. Each profile controls grid resolution, frame count, and output cadence while keeping the downstream pipeline identical.

## 2. Solve from the command line

```bash
./build/motor_sim --scenario inputs/three_phase_stator_ci.json --solve --vtk-series outputs/three_phase_ci.pvd
```

Add `--parallel-frames` to solve timeline entries concurrently. Warm starts (`--warm-start`) reuse the previous frame’s solution, typically halving CG iterations after the first frame. Prolongation (`--use-prolongation` with `--coarse-nx/--coarse-ny`) seeds fine grids from a coarse solve to accelerate convergence.

Outputs inherit `_frame_###` suffixes when a `timeline` is present, avoiding collisions across frames.

## 3. Timeline authoring

Augment a scenario with a `timeline` array to prescribe per-frame changes:

```json
"timeline": [
  {"t": 0.0, "wire_currents": [10.0, -10.0]},
  {"t": 0.001, "wire_currents": [5.0, -5.0], "rotor_angles": {"rotor": 90.0}},
  {
    "t": 0.002,
    "wires": [{"index": 0, "current": 0.0}, {"index": 1, "scale": 0.25}],
    "magnets": [{"index": 0, "angle_deg": 45.0}],
    "rotor_angles": [{"name": "aux_rotor", "angle_deg": -15.0}]
  }
]
```

Each entry expands into an independent frame. Declare rotors under the top-level `rotors` key to group polygons, magnets, and wires that move together. Magnet overrides accept either rotation angles or explicit magnetisation vectors.

## 4. Inspecting results

Open the emitted `.pvd` series in ParaView to explore fields and outlines. The generator emits:

- `outputs/three_phase_frame_###.vti`: cell-centred B/H fields per frame.
- `outputs/three_phase_ci.pvd`: ParaView time-series index (from `--vtk-series`).
- `outputs/three_phase_outlines.vtp`: geometry overlays for slots, stator surfaces, and wires.
- `outputs/bore_angle.csv`: bore-average flux magnitude and angle traces.

Use `python/animate_three_phase.py` to create annotated animations, or the `python/visualize_scenario_field.py` helper for static renders with streamlines and overlays.

[Open in GUI](../developer-guide/dev-environment.md){ .md-button }

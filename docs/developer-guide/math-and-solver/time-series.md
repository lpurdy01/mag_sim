# Time-Series Frames and Parallel Solving

The solver can execute quasi-static studies by expanding a scenario's `timeline`
array into independent frames. Each frame is a snapshot of the problem with its
own wire currents and magnetisation overrides. The solver rasterises and solves
these frames sequentially by default, or concurrently when `--parallel-frames`
is supplied. For an end-to-end walkthrough see the [Running Simulations guide](../../user-guide/running-simulations.md);
this page focuses on schema details and developer diagnostics.

## Authoring a Timeline

Add a `timeline` array to a scenario JSON file. Each entry is an object that can
contain:

- `t` or `time`: Optional timestamp stored for reference and reporting.
- `wire_currents`: Array of absolute currents matching the order of the `sources`
  list. Alternatively, provide a `wires` array with objects containing
  `index` and either `current`, `I`, or `scale` (relative multiplier).
- `rotor_angle` / `rotor_angle_deg`: Legacy shorthand that rotates the first
  rotor when `rotors` are declared, or all magnet regions when no rotors exist.
- `rotor_angles`: Fine-grained rotor overrides. Accepts either an array of
  angles (aligned with the order of the `rotors` list), an array of objects with
  `index`/`name` and `angle`/`angle_deg`, or an object mapping rotor names to
  angles in degrees.
- `magnets`: Array of per-region overrides, each with an `index` and either an
  `angle_deg` rotation or a `magnetization` vector (`[Mx, My]` array or object
  with `Mx`/`My`).

### Rotor definitions

Declare rigid rotor groups in the root of the scenario JSON under `"rotors"`.
Each entry can specify:

- `name`: Optional identifier used by `rotor_angles` and metadata.
- `pivot`: Two-element array `[x, y]` describing the rotation centre (defaults
  to the domain origin if omitted).
- `polygon_indices`: Zero-based indices into the polygon region list that
  should rotate with the rotor.
- `magnet_indices`: Magnet region indices that move and spin with the rotor.
- `wire_indices`: Source indices (wires) that orbit the pivot.

Polygons and magnets carry their geometry through the rotation, ensuring both
field rasterisation and exported outlines match the transformed shape. Wires are
repositioned before the solve so current deposition tracks the new location.

Example:

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

Each timeline entry produces a `ScenarioFrame` with its own copy of the base
specification. The `expandScenarioTimeline` helper returns the expanded list of
frames for use in tests and tooling.

## CLI Usage

Run the solver across all frames with:

```bash
./build/motor_sim --scenario inputs/your_scenario.json --solve
```

When a timeline is present, outputs automatically gain a `_frame_XXX` suffix to
avoid overwriting results. Enable multi-threaded solving with:

```bash
./build/motor_sim --scenario inputs/your_scenario.json --solve --parallel-frames
```

The solver launches up to `hardware_concurrency() - 1` worker threads, capped by
the number of frames.

### Warm starts and prolongation

Timeline solves can reuse previous solutions to accelerate convergence. Passing `--warm-start` retains the final \(A_z\) field
from the preceding frame and feeds it to the next CG call as the initial guess. For sudden changes (e.g., when switching from a
coarse exploratory solve to a finer production grid), combine `--use-prolongation` with optional `--coarse-nx/--coarse-ny`
overrides to seed the fine grid from a cheap coarse solve. These options materially reduce CG iteration counts while matching
the SOR baseline to within the regression tolerances.

## Output Naming

For single-frame scenarios (no `timeline` key), output paths are respected as
written. When a timeline is present, each field map, line probe, and midline CSV
appends `_frame_###` (three digits by default) before the extension, e.g.
`outputs/flux_map_frame_002.csv`.

## Testing Support

`tests/timeline_frames_test.cpp` ensures timeline parsing, wire overrides,
rotor-driven geometry motion, and magnet rotations behave predictably. The
accompanying JSON fixture lives in
`inputs/tests/timeline_frames_test.json`.

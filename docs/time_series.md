# Time-Series Frames and Parallel Solving

The solver can execute quasi-static studies by expanding a scenario's `timeline`
array into independent frames. Each frame is a snapshot of the problem with its
own wire currents and magnetisation overrides. The solver rasterises and solves
these frames sequentially by default, or concurrently when `--parallel-frames`
is supplied.

## Authoring a Timeline

Add a `timeline` array to a scenario JSON file. Each entry is an object that can
contain:

- `t` or `time`: Optional timestamp stored for reference and reporting.
- `wire_currents`: Array of absolute currents matching the order of the `sources`
  list. Alternatively, provide a `wires` array with objects containing
  `index` and either `current`, `I`, or `scale` (relative multiplier).
- `rotor_angle` / `rotor_angle_deg`: Rotate all magnet regions by the given
  angle in degrees relative to their base magnetisation vector.
- `magnets`: Array of per-region overrides, each with an `index` and either an
  `angle_deg` rotation or a `magnetization` vector (`[Mx, My]` array or object
  with `Mx`/`My`).

Example:

```json
"timeline": [
  {"t": 0.0, "wire_currents": [10.0, -10.0]},
  {"t": 0.001, "wire_currents": [5.0, -5.0], "rotor_angle": 90.0},
  {
    "t": 0.002,
    "wires": [{"index": 0, "current": 0.0}, {"index": 1, "scale": 0.25}],
    "magnets": [{"index": 0, "angle_deg": 45.0}]
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

## Output Naming

For single-frame scenarios (no `timeline` key), output paths are respected as
written. When a timeline is present, each field map, line probe, and midline CSV
appends `_frame_###` (three digits by default) before the extension, e.g.
`outputs/flux_map_frame_002.csv`.

## Testing Support

`tests/timeline_frames_test.cpp` ensures timeline parsing, wire overrides, and
magnet rotations behave predictably. The accompanying JSON fixture lives in
`inputs/tests/timeline_frames_test.json`.

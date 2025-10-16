# Back-EMF Probes

Back-electromotive-force (back-EMF) measurements estimate the voltage induced in
coils by changing magnetic flux across timeline frames. A back-EMF probe samples
the field map inside a user-defined loop, records the flux for each solved frame,
and post-processes consecutive pairs to produce \(-\Delta \Phi / \Delta t\).

## JSON schema

Declare probes in the scenario `outputs` array:

```json
{
  "type": "back_emf_probe",
  "id": "stator_phase_a",
  "component": "Bx",
  "region": {
    "type": "polygon",
    "vertices": [
      [-0.02, -0.01],
      [0.02, -0.01],
      [0.02, 0.01],
      [-0.02, 0.01]
    ]
  },
  "frames": [0, 1, 2, 3]
}
```

* `component` chooses which magnetic flux-density component to integrate.
  Supported values are `"Bx"`, `"By"`, and `"Bmag"` (the default). Use `Bx`
  or `By` when the coil surface is predominantly pierced by one in-plane field
  component.
* `region` describes the integration area. Provide either a polygon (three or
  more vertices) or a rectangle via `{"type": "rect", "x_range": [...],
  "y_range": [...]}`.
* `frames` is optional. When omitted the probe spans every frame in the
  timeline. Provide at least two frame indices when you need a sparse subset.
* `path` is optional; by default results are written to
  `outputs/<id>_emf.csv`.

## Output format

Back-EMF probes write a single CSV containing one row per interval between
frames:

```text
frame_start,frame_end,t_start,t_end,delta_t,flux_start,flux_end,emf
0,1,0.000000000000e+00,1.000000000000e-03,1.000000000000e-03,0.000000000000e+00,2.500000000000e-03,-2.500000000000e+00
```

Columns describe the contributing frames, their timestamps, the integrated flux
(values are in webers per metre), and the induced voltage per unit length in
volts per metre. When timelines are absent or only one frame is available, the
solver reports an error because \(\Delta t\) would be undefined.

## Python API

The `BackEmfProbeOutput` helper mirrors the JSON schema:

```python
from scenario_api import BackEmfProbeOutput

output = BackEmfProbeOutput(
    id="phase_a",
    component="By",
    vertices=[(-0.01, -0.02), (0.01, -0.02), (0.01, 0.02), (-0.01, 0.02)],
    frames=[0, 2, 4],
)
```

Combine with `Scenario.outputs.append(output)` before calling `save_json`.

## Numerical notes

* Flux is approximated by averaging the requested component over cell centres
  inside the region and multiplying by the cell area (per-unit-length measure).
* When `component="Bmag"` the solver integrates the magnitude of **B**, which is
  helpful when the dominant component varies with rotor angle but introduces an
  approximation.
* Timeline timestamps should be monotonic and distinct; the solver raises an
  error if two consecutive frames share the same time value.
* Back-EMF probes currently operate on field-map data. They do not capture coil
  turn counts or end-winding effects — multiply the reported voltage per unit
  length by the active length and number of turns to obtain line voltages. The
  three-phase PM motor demo integrates the flux magnitude over each phase’s
  positive slot to produce a sinusoidal reference waveform that lines up with
  the rotating air-gap field.
* `tests/back_emf_probe_test.cpp` drives the integration helper with synthetic
  three-phase flux waveforms, verifying that the induced EMFs follow the
  expected sinusoidal shape, maintain 120° phase separation, and scale with both
  flux amplitude and electrical speed.

# Force and Torque Probes

Magnetic machines often require estimates of the mechanical force or torque that
the air-gap fields apply to a rotor. `mag_sim` now exposes a **stress tensor
probe** that integrates the Maxwell stress tensor along a closed contour, giving
per-unit-length force and torque values directly from the solved field map.

## Maxwell Stress Tensor refresher

For magnetostatics the Maxwell stress tensor in Cartesian coordinates is

\[
T_{ij} = \frac{1}{\mu_0} \left(B_i B_j - \tfrac{1}{2} \delta_{ij} B^2\right),
\]

with the traction on a surface of unit normal **n** given by

\[
\mathbf{f} = T \cdot \hat{\mathbf{n}}.
\]

Integrating the traction around a closed 2D loop yields the net force per unit
length, while the torque about the out-of-plane axis follows from

\[
\tau_z = \oint (x f_y - y f_x)\,\mathrm{d}s.
\]

The solver evaluates these expressions numerically using midpoint sampling on
polygon edges. Bilinear interpolation recovers **B** at arbitrary contour points
so probes are not restricted to grid-aligned loops.

## JSON schema

Declare probes inside the scenario `outputs` array:

```json
{
  "type": "probe",
  "id": "rotor_mst",
  "probe_type": "force_and_torque",
  "method": "stress_tensor",
  "loop": {
    "type": "polygon",
    "vertices": [
      [-0.03, -0.02],
      [0.03, -0.02],
      [0.03, 0.02],
      [-0.03, 0.02]
    ]
  },
  "path": "outputs/rotor_mst.csv"
}
```

* `probe_type` chooses which quantities are of interest: `"force"`,
  `"torque"`, or `"force_and_torque"` (all three values are written either way).
* `method` currently supports only `"stress_tensor"`.
* `loop` supplies at least three vertices describing the closed contour. Provide
  either an array of `[x, y]` pairs or an object with `type="polygon"` and a
  `vertices` array as shown above.

The solver writes a CSV with a single row:

```text
Fx,Fy,Tz
-1.234567890123e+02,5.678901234567e+01,-2.468013579240e-03
```

Values represent force (newtons per metre) and torque (newton-metres per metre)
about the global origin in SI units. When timelines are active the file name is
extended with the usual `_frame_###` suffix.

## Usage tips

* Keep probe contours well inside the simulation domain so bilinear
  interpolation never samples outside the grid.
* Dense grids yield smoother estimates. For coarse meshes consider slightly
  inflating the loop to avoid sampling immediately adjacent to discretisation
  artefacts.
* Combine with timeline frames to capture torque ripple across electrical angles
  or to cross-check against virtual-work calculations in future releases.

The `tests/probe_output_test.cpp` fixture exercises the ingestion and evaluation
path with a synthetic field that generates a known downward force, providing a
regression guard for the new feature.

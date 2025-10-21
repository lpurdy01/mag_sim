# VTK Output

`motor_sim` writes VTK ImageData (`.vti`) for field maps and optional PolyData (`.vtp`) outlines.

## Image data (`.vti`)

| Array | Location | Description |
| ----- | -------- | ----------- |
| `B` | Cell data | Magnetic flux density vector magnitude and components. |
| `H` | Cell data | Magnetic field intensity derived from `B` and material properties. |
| `J` | Cell data | Impressed current density (per frame). |
| `region_index` | Cell data | Material region identifier used for visual overlays. |

Open the sequence through the generated `.pvd` index when `--vtk-series` is passed on the CLI.

## Geometry outlines (`.vtp`)

PolyData outlines contain per-cell attributes for `kind`, `loop_index`, and optional rotor or stator labels. Join the accompanying `*_outlines_labels.csv` file in ParaView to expose human-readable names when colouring or filtering loops.

For visualization workflows see [Visualization Overview](../../user-guide/visualization/overview.md).

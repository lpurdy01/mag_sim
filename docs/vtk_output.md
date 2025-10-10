# VTK Field Map Exports

The solver can emit full-field snapshots in VTK XML ImageData (`.vti`) format.
VTK files provide:

- Compact binary storage suitable for large grids
- Direct visualisation in ParaView, VisIt, and similar tools
- Straightforward ingestion from Python via `vtk` and NumPy

## Requesting a VTK export

Add a field-map output with `"format": "vti"` to a scenario JSON (or use
`scenario_api.FieldMapOutput(format="vti")`). When no explicit `path` is
provided the ingestor writes to `outputs/<id>.vti`.

```json
{
  "type": "field_map",
  "id": "field_frame_000",
  "quantity": "B",
  "format": "vti"
}
```

Each VTK file stores cell-centred datasets named `Bx`, `By`, `|B|`, `Hx`, `Hy`,
`|H|`, and `energy_density`. Values are computed by averaging the four
surrounding nodes before evaluating magnitudes and the energy density
(\(\tfrac{1}{2}\,\mathbf{B}\cdot\mathbf{H}\)).

## Verifying the output

Use the helper script to confirm ParaView-readable structure and non-empty
arrays:

```bash
python python/verify_vtk.py outputs/field_frame_000.vti
```

On success the script prints `VTK verification: PASSED`. The implementation uses
raw-appended binary payloads with 64-bit length headers, so files remain small
for large grids while keeping the XML metadata legible.


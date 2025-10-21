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

Each VTK file stores cell-centred datasets named `B`, `Bx`, `By`, `|B|`, `H`,
`Hx`, `Hy`, `|H|`, and `energy_density`. The `B` and `H` arrays are true vector
fields (three components with a zero Z entry) so ParaView can render streamlines
or glyphs without manual component selection. Scalar components are still
available for quick inspection. Values are computed by averaging the four
surrounding nodes before evaluating magnitudes and the energy density
(\(\tfrac{1}{2}\,\mathbf{B}\cdot\mathbf{H}\)).

Alongside each `.vti` export the solver writes a
`*_outlines.vtp` PolyData companion that contains the domain rectangle, material
polygons, magnet outlines, and wire circles as closed polylines. ParaView can
load the outline file as a second source and overlay it on the field map using a
`Glyph`/`Tube` filter or simple line rendering. Cell data include a `kind`
integer (`0=domain`, `1=material`, `2=magnet`, `3=wire`) and a `loop_index`
identifier. A CSV with the same stem and `_labels.csv` suffix records the
mapping from `loop_index` to human-readable labels **and** a `group` column used
for rotor assemblies or other geometry groupings, making it easy to drive
ParaView filters or colour maps via simple joins in the Spreadsheet view or
external tooling.

## Verifying the output

Use the helper script to confirm ParaView-readable structure and non-empty
arrays (including the combined vector fields):

```bash
python python/verify_vtk.py outputs/field_frame_000.vti
```

On success the script prints `VTK verification: PASSED`. The implementation uses
raw-appended binary payloads with 64-bit length headers, so files remain small
for large grids while keeping the XML metadata legible.


[Open in GUI](../../developer-guide/dev-environment.md){ .md-button }

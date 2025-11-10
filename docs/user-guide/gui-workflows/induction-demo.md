# Induction motor GUI workflow

This walkthrough demonstrates how to assemble the bundled three-phase induction
motor demo inside the Flask GUI, run the transient solve, and review the
resulting field maps and animation. It assumes you have already built
`motor_sim` and started the Flask server as described in the
[GUI overview](../gui_flask.md).

## Assets

| Purpose | File | Notes |
| --- | --- | --- |
| Scenario JSON | [`inputs/induction_gui_demo.json`](../../inputs/induction_gui_demo.json) | Defines the transient spin-up with field-map, VTK series, and probe outputs. |
| Domain outline | [`docs/assets/dxf/induction_demo/domain.dxf`](../../assets/dxf/induction_demo/domain.dxf) | Axis-aligned boundary for the rectangular workspace. |
| Stator laminations | [`docs/assets/dxf/induction_demo/stator_only.dxf`](../../assets/dxf/induction_demo/stator_only.dxf) | Single layer containing the stator steel outline. |
| Rotor back-iron | [`docs/assets/dxf/induction_demo/rotor_core.dxf`](../../assets/dxf/induction_demo/rotor_core.dxf) | Rotor core polygon. |
| Rotor bars | [`docs/assets/dxf/induction_demo/rotor_bars.dxf`](../../assets/dxf/induction_demo/rotor_bars.dxf) | Conductor bars, one polyline per bar. |
| Air references | [`docs/assets/dxf/induction_demo/air_regions.dxf`](../../assets/dxf/induction_demo/air_regions.dxf) | Convenience overlay to confirm the air gap and slot cavities. |

> Tip: if you tweak the scenario and want to regenerate the DXF set, choose
> **File → Export DXF set…** after the run. The GUI will emit fresh DXFs for the
> loaded scenario using the same layering rules described above.

## Step 1 – Prime the workspace

1. Launch the server with `python -m python.gui.app_flask` and open the
   interface in your browser.
2. In the **File** menu choose **New project** if a previous session is still
   loaded.
3. Download the assets listed above (or keep the repository checked out so the
   relative links work directly).

## Step 2 – Import CAD layers

Use the **CAD workspace** panel to upload each DXF. The recommended order helps
keep the combined preview tidy, but the GUI will happily accept any sequence.

| Upload order | Layer selection | Suggested category |
| --- | --- | --- |
| `domain.dxf` | Single `domain` layer | `domain` |
| `stator_only.dxf` | Entire layer | `material` (stator steel) |
| `rotor_core.dxf` | Entire layer | `material` (rotor core) |
| `rotor_bars.dxf` | Entire layer | `wire` |
| `air_regions.dxf` | Optional visual guide | `structural` |

After each upload, expand the card, tick the `Include` checkboxes, and assign
categories from the drop-down. Click **Update mapping** to persist the
selection, then confirm the coloured overlay in the **Combined preview** on the
right. When everything is aligned the stator/rotor teeth should nest cleanly
within the domain boundary.

## Step 3 – Load the scenario

1. Scroll to **Simulation setup** and use **Browse…** to select
   `inputs/induction_gui_demo.json`.
2. Leave the solver defaults (CG, `1e-6` tolerance, `20000` max iterations)
   intact. The scenario already requests the `induction_motor_field` field-map
   along with the VTK series and probes, so the `--outputs` box can remain
   blank.
3. Click **Run simulation**. The progress card and live log will appear and the
   **Run** button will grey out while the solver works.

## Step 4 – Monitor progress and animation

- The live log shows each timeline frame as it converges. When the solver writes
  a field map (e.g. `Frame 0: wrote field_map 'induction_motor_field' …`) the
  GUI renders a preview immediately.
- The **Results & visuals** card keeps a gallery of these frames. Use the
  **Play** button to loop through the snapshots, **Prev/Next** to step manually,
  or drag the slider to inspect a specific frame. The caption updates with the
  frame index so you can correlate it with the log.
- You can re-render the final field map with different overlays at any time by
  adjusting the controls (vector scaling, log floor, streamlines) and clicking
  **Update visualisation**.

## Step 5 – Collect artefacts

When the run completes successfully the **Downloads** list exposes:

- The uploaded scenario JSON (with any edits from the DXF workspace preserved).
- The solver log.
- `outputs/induction_motor_field.csv` plus the other CSV/VTP artefacts declared
  in the scenario.
- A PNG of the default field-map render.

These files, along with the animation controls, are enough to document the
induction motor run or to feed the CSVs into the downstream plotting tools such
as `python/animate_three_phase.py` if you need a high-resolution movie.

## Next steps

- Save the current project through **File → Save scenario…** if you tweak the
  geometry or solver settings. The saved JSON can be re-uploaded later to keep
  iterating.
- Use **File → Export DXF set…** to share the geometry with teammates or to
  seed another CAD package.
- For heavier visualisations (rotor slot colouring, full GIF exports) hand the
  downloads to the existing scripts in `python/`—the GUI keeps the workspace
  clean while providing a quick feedback loop.

# Iron ring GUI workflow

This walkthrough shows how to recreate and run the bundled
[`iron_ring_demo.json`](https://raw.githubusercontent.com/mag-co/mag_sim/main/inputs/iron_ring_demo.json)
scenario directly from the Flask GUI. It focuses on three goals:

1. Importing DXF geometry (domain, iron body, windings) and keeping the layers
   organised even when your CAD exporter splits stator, rotor, and coil data
   into separate files.
2. Building the scenario JSON graphically via the **Materials & regions**,
   **Windings**, and **Timeline** stages so you can tweak permeability and
   stimulus settings without leaving the browser.
3. Launching the solve, reviewing the automatically generated field-map frames,
   and downloading both the JSON and DXF exports for future iterations.

> **Tip:** The interface lives at
> <https://lpurdy01.github.io/mag_sim/user-guide/gui_flask/>. The workflow
> pages (including this one) render on the same documentation site so you can
> open them next to the GUI while you work.

## Prerequisites

- A `motor_sim` build in `./build/motor_sim`.
- The GUI dependencies installed via `./scripts/setup_gui_env.sh`.
- DXF exports of your ring: either a single all-in-one file such as
  `custom_target_v0.dxf` or multiple files (domain, iron, windings). Upload in
  the order shown below so the combined preview layers stack predictably.

| Upload order | Example layers                       | Recommended category |
| ------------ | ------------------------------------ | -------------------- |
| 1            | `DOMAIN`, `REFERENCE`                | Domain boundary      |
| 2            | `IRON_RING`, `CORE`, `BACKIRON`      | Material region      |
| 3            | `AIR_GAP`, `SLOT_AIR` (if exported)  | Material region      |
| 4            | `COIL_A+`, `COIL_A-`, `COIL_B+`, …   | Conductor path       |

If your CAD software emits one file per subsystem (stator, rotor, windings),
import the stator/core first, then any air-gap cut-outs, and finally the wire
sets. Every upload can mix layers from different categories—you simply map each
layer to the desired role via the table in the CAD workspace.

## Step 1 – Import and map DXF layers

1. Open the GUI and use the **File → New project** menu entry to start fresh.
2. In the **CAD workspace**, upload your domain DXF. Keep the layer selected and
   set its category to **Domain boundary**. Click **Update mapping** to redraw
   the combined preview.
3. Upload the iron body DXF (or the same file again if everything lives in
   `custom_target_v0.dxf`). Mark the iron layer(s) as **Material region** and
   deselect any helper layers you do not want in the simulation.
4. Upload the winding DXF(s) last. Mark each conductor layer as
   **Conductor path** so the **Windings** stage can reference them later.
5. Refresh the combined preview after each change to confirm the stack order
   looks correct (domain outline, iron ring, then the copper overlays).

## Step 2 – Define the material palette

Switch to **Materials & regions**. The iron ring demo only needs two entries,
matching the JSON snippet below:

```json
"materials": [
  {"name": "air", "mu_r": 1.0},
  {"name": "iron", "mu_r": 1000.0}
]
```

The GUI ships with `air` and `steel` placeholders; rename the second row to
`iron` and adjust the permeability (or add extra rows for laminations,
aluminium, etc.). Click **Save materials** to persist the palette. The DXF layer
category pickers immediately use the updated names.

## Step 3 – Map conductor layers to windings

Open the **Windings & conductor mapping** stage. Each row becomes a
`current_region` in the scenario JSON.

1. Click **Add winding** until you have a row for each phase (e.g. Phase A, B,
   C). Name the row, set the phase ID, and leave turns/fill at their defaults or
   edit them to match your stack-up.
2. Use the wire-layer multi-select to assign the DXF layers you tagged as
   **Conductor path**. Hold Ctrl/Cmd to pick both the positive and negative
   sides of a phase if they were exported as separate layers.
3. Choose the orientation so clockwise windings are “into” the page and
   counter-clockwise windings are “out”. Orientation simply multiplies the
   generated `current_region` orientation value.
4. Click **Apply windings**. The GUI rasterises the DXF geometry, rewrites the
   `sources` array, and stores the winding configuration in the project. Use
   **File → Save scenario…** at any time to download the JSON with the generated
   current regions.

## Step 4 – Configure the drive timeline

The iron ring demo uses six DC wires, so you can either keep the default static
setup (leave the timeline empty) or experiment with a transient drive:

- **Balanced three-phase:** Set the peak current (e.g. 30 A), frequency (e.g.
  50 Hz), steps per cycle, and the phase sequence that matches your winding
  order. Click **Save timeline** and the GUI will generate a `timeline` array
  plus matching `transient.dt`/`n_steps` fields.
- **Manual JSON:** Paste the frames produced by a script or spreadsheet. Each
  entry must include a time stamp `t` and a `phase_currents` mapping.

The alert inside the card summarises both modes and links back to this page if
you need to cross-reference while editing.

## Step 5 – Run, animate, and download

1. In **Simulation setup**, either upload `iron_ring_demo.json` or click
   **Save scenario…** and re-upload the project’s JSON. Pick the solver (`cg` is
   typical for this case), tolerance, and iteration cap.
2. Click **Run simulation**. The log and progress panels stream updates; every
   `field_map` output reported by the solver automatically lands in the
   **Results visualisation** gallery.
3. Use the playback controls to loop the frames or drag the slider to inspect a
   single snapshot. The **Update visualisation** form regenerates the PNG with
   different arrow density, colour scales, or overlays without rerunning the
   solver.
4. Download artefacts from the **Results & downloads** panel. The solver log,
   field-map CSVs, DXF exports, and the project JSON are staged into
   `python/gui/results/` so the **Download** buttons work even when the original
   files lived elsewhere.

## Troubleshooting checklist

- If the windings form is empty, ensure at least one DXF layer is marked as
  **Conductor path** in the CAD workspace.
- If the timeline inputs are disabled, upload or preview a scenario first so
  the project has somewhere to store the generated frames.
- When a solver log reports `outputs/iron_ring_field.csv` but the gallery stays
  empty, confirm the CSV exists under the repository’s `outputs/` directory. The
  GUI resolves relative paths against both the scenario folder and the solver’s
  working directory, so a simple `ls outputs` is often enough to confirm.
- Use **File → New project** to reset the workspace if you want to experiment
  with a different DXF set without carrying over the previous windings or
  timeline.

# Web GUI (Flask)

The Flask-based GUI provides a lightweight browser front-end for running
`mag_sim` scenarios. It mirrors the functionality explored in the Streamlit
prototype while staying close to the underlying command-line workflow.

## Requirements

- A built `motor_sim` executable in `./build/motor_sim`.
- Flask, Matplotlib, NumPy, and ezdxf available in the active Python
  environment. Execute `./scripts/setup_gui_env.sh` to install the toolchain
  (and create the GUI runtime folders) when a pre-provisioned devcontainer is
  not available.

## Starting the server

From the repository root:

```bash
python -m python.gui.app_flask
```

The development server listens on `http://127.0.0.1:5000` by default. When
running inside GitHub Codespaces or a similar environment, forward port 5000 to
access the interface from your browser.

Alternatively, set `FLASK_APP=python.gui.app_flask` and run `flask run` if you
prefer Flask's CLI wrapper.

## Using the interface

The CAD-inspired layout splits the page into a top navigation bar, a workflow
sidebar, and content panels for each stage. The **File** menu exposes
project-wide actions (**New project** clears the workspace, **Save scenario…**
downloads the current JSON, **Export DXF set…** produces layered DXFs for the
current scenario), while the sidebar links jump to the key panels without
leaving the page.

Looking for a guided example? Follow the hosted walkthroughs:

- [Induction motor GUI workflow](https://lpurdy01.github.io/mag_sim/user-guide/gui-workflows/induction-demo/)
  – multi-file DXF import and transient animation.
- [Iron ring GUI workflow](https://lpurdy01.github.io/mag_sim/user-guide/gui-workflows/iron-ring-demo/)
  – mapping custom DXF layers to materials, windings, and timelines.

Work through the panels in order:

1. **CAD workspace:** Import DXF files, mark each layer with a category (domain,
   material, magnet, conductor), and update the combined preview. Upload domain
   references first, then materials/magnets, and finally wire layers. Each file
   can contain multiple categories; the form lets you toggle inclusion per
   layer.
2. **Materials & regions:** Define the permeability palette that the scenario’s
   `materials` block uses. The table feeds the DXF layer drop-downs so your CAD
   mapping stays in sync with the JSON.
3. **Windings & conductor mapping:** Create one row per phase/coil group,
   configure turns/fill/orientation, and assign the conductor layers that were
   tagged in the CAD workspace. Saving the form rasterises the DXF geometry into
   `current_region` entries and stores the winding metadata in the project.
4. **Timeline designer:** Generate balanced three-phase drive patterns or paste
   manual timeline JSON. The GUI updates both the `timeline` array and the
   `transient` parameters so subsequent CLI runs remain consistent.
5. **Geometry preview:** Render the scenario’s domain/magnets without solving to
   confirm the JSON matches the CAD stack-up.
6. **Simulation setup:** Upload or reuse the scenario JSON, tweak solver
   settings (`cg`/`sor`, tolerance, iteration cap, `--outputs`), and launch the
   solve. Use the alert inside the card as a reminder that the sections above
   control the graphical scenario composer.
7. **Progress, results, and downloads:** Watch the SSE-driven progress bar and
   log, review the live field-map gallery, and download the scenario/log/DXF
   exports/field-map CSVs that the run produced.

Only one solve is allowed at a time; the interface reports an error if you try
to launch a second run before the first finishes. Use **Stop** to request
termination—`terminate()` is sent to the subprocess and the log records the
request before the child exits.

### Environment upkeep

Generated artefacts accumulate under `python/gui/uploads` and
`python/gui/results`. Run `./scripts/maintain_gui_env.sh` periodically (set
`KEEP_DAYS` or pass a day count to keep more history) to prune stale files and
surface any outdated Python packages.

### Visualisation controls

The results card includes a lightweight form for regenerating the field-map
image without leaving the browser:

- **Vector scaling**: switch between linear arrows, logarithmic scaling, or hide
  vectors entirely.
- **Colour scale**: toggle between linear and logarithmic magnitude mapping.
- **Arrow skip**: control the quiver density (defaults to every fourth sample).
- **Log floor / Vector floor**: clamp low values when using logarithmic modes.
- **Region outlines / Streamlines**: enable or disable material/magnet borders
  and streamline overlays.

Choose the desired parameters and click **Update visualisation**. The server
re-renders the image with Matplotlib and returns a fresh PNG via the
`/visualization.png` endpoint, keeping the UI responsive with minimal
JavaScript.

During a simulation the GUI also builds a frame gallery from each reported
`field_map`. Use the playback controls above the image to loop the snapshots,
step frame-by-frame, or drag the slider to inspect a single moment. The caption
updates with the frame index, making it easy to cross-reference the solver log
without leaving the page.

### Materials, windings, and timeline composer

The three new stages share the same persistence model:

- Saving the **Materials** form rewrites the `materials` array in the current
  project JSON and immediately feeds those entries back into the DXF layer
  pickers.
- The **Windings** form stores a project-local winding list and generates
  `current_region` sources by rasterising the selected DXF layers. Every selected
  layer can be reused across multiple windings (the UI simply records the
  token). Orientation switches flip the sign of the generated `orientation`
  field.
- The **Timeline** designer stores both the generated frames and the inputs used
  to create them so you can revisit the balanced parameters or edit the raw JSON
  later.

Whenever you change any of these composer stages, use **File → Save scenario…**
to export the updated JSON, or simply press **Run simulation** to drive
`motor_sim` with the newly generated spec.

## Known limitations

- DXF ingestion currently supports LINE/POLYLINE/LWPOLYLINE and CIRCLE entities.
  Splines, text, and block references are ignored. Use your CAD tool to explode
  complex geometry before uploading.
- Preview and visualisation rendering rely on Matplotlib, NumPy, and ezdxf.
  These dependencies ship with the devcontainer, but a custom environment must
  provide the same stack.
- The progress parser uses a simple percentage heuristic. If the solver log does
  not emit percentage tokens the bar will stay at 0 % until the run completes.
- The current prototype keeps global state and targets single-user usage. A
  production deployment should move to per-session queues and authentication.

## Related work

A Streamlit-based GUI prototype lives on the `feat/gui-streamlit` branch. Both
paths explore similar workflows; use this Flask version when you need explicit
control over routing and templating, or switch branches to compare the
Streamlit-powered experience.

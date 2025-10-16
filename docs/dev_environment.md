# Developer Environment & VS Code Workflow

This guide captures the day-to-day tooling expected for contributors. It assumes
GitHub Codespaces, VS Code Remote (WSL/SSH), or a local Linux workstation with
CMake ≥ 3.16 and GCC ≥ 11.

## 1. Command-line quickstart

```bash
scripts/setup_env.sh           # optional helper for Codespaces/WSL
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build -j
ctest --test-dir build --output-on-failure
```

The build directory is never committed (`.gitignore` keeps it clean). Tests live
under `tests/` and are surfaced through CTest so they can be executed both from
the terminal and IDE integrations.

### 1.1 Iteration-friendly test cadence

For day-to-day work, favour the smallest slice of the suite that exercises your
changes before falling back to the full battery:

* Use `ctest --output-on-failure -R <regex>` (or invoke an individual binary
  such as `./build/torque_validation_test`) while iterating on a specific
  feature. This keeps turnaround low and mirrors the CI harness.
* When touching the solver or scenario ingestion, rerun the relevant timeline
  with `./build/motor_sim ...` to double-check CSV/VTK outputs without waiting
  for unrelated regressions.
* Before publishing a branch or opening a PR, always follow up with
  `ctest --output-on-failure --parallel $(nproc)` (or
  `scripts/run_ci_checks.sh`) so the code you push matches the CI matrix.

`ctest --parallel` honours the `CTEST_PARALLEL_LEVEL` environment variable; the
CI helper script defaults it to `$(nproc)` so local runs use all available
hardware. Set it explicitly (e.g. `CTEST_PARALLEL_LEVEL=4`) if you need to
reserve cores for other tasks.

Key runtime flags for `motor_sim`:

* `--solver {sor|cg}` toggles between the legacy Gauss–Seidel solver and the new preconditioned conjugate gradient (default is
  SOR for continuity).
* `--warm-start` reuses the previous frame's field as the initial guess when traversing a timeline, dramatically cutting CG
  iterations.
* `--use-prolongation` seeds the fine-grid solve from an automatically selected coarse solve; adjust with `--coarse-nx/--coarse-ny`.
* `--progress-every <seconds>` controls the live progress cadence (default 2 s). Set it to `0` to emit a sample every iteration
  when collecting detailed residual histories.
* `--snapshot-every <iters>` enables downsampled field dumps requested by progress sinks (handy for spotting spatial stagnation).
* `--progress-history <path>` writes the emitted residual samples to a CSV. Timeline runs append `_frame_###` to the basename;
  combine with `--progress-every 0` for a full per-iteration log.
* `--quiet` suppresses progress output when scripting multiple runs.

## 2. VS Code configuration

The repo ships with `.vscode/` settings tuned for the `CMake Tools` and `C/C++`
extensions.

### 2.1 Tasks

Use `Terminal → Run Task` (or the command palette) to access the predefined
tasks:

| Task label               | Purpose                                                                  |
| ------------------------ | ------------------------------------------------------------------------ |
| `cmake-configure`        | Configure the out-of-source build directory in Debug mode.               |
| `cmake-build`            | Build all targets; depends on configuration.                             |
| `ctest`                  | Run the full test suite with verbose failure output.                     |
| `run-analytic-wire-test` | Compile + launch the analytic validation binary in place.                |
| `run-solver-benchmark`   | Compile + execute the benchmark helper with a lean 129×129 grid preset. |
| `run-two-wire-scenario`  | Build + execute `motor_sim` against the bundled JSON scenario.           |

These tasks chain automatically (for example, running the benchmark will trigger
a build if required). Feel free to duplicate them locally for custom grids or
release builds.

### 2.2 Debugging

The `Run and Debug` panel exposes launch configurations for four executables:

1. **motor_sim** — loads JSON scenarios; defaults to printing usage when no args are supplied.
2. **Run two-wire scenario** — launches `motor_sim` with the bundled cancellation case and writes the midline CSV.
3. **analytic_wire_test** — deterministic regression for the single-wire case.
4. **solver_benchmark** — convenience harness for measuring throughput.

All launchers invoke the `cmake-build` task beforehand to guarantee the binary is
up to date. Adjust command-line arguments via the `args` array in
`.vscode/launch.json` if you need alternative solver tolerances or grids.

### 2.3 IntelliSense hints

* Compiler path defaults to `/usr/bin/g++` (set in `.vscode/settings.json`).
* `compile_commands.json` is generated automatically inside `build/` by CMake;
  the C/C++ extension will pick it up to power IntelliSense if you add the build
  directory to VS Code's workspace settings (the default setup already does this
  via the CMake Tools extension).

## 3. Codespaces resource guidance

Codespaces typically provide 2 vCPUs and 4 GB RAM on the basic tier. To keep
turnaround snappy:

* Prefer grid sizes ≤ 257² for exploratory runs. Larger domains converge but can
  exceed a minute when `tol = 1e-6`.
* Disable verbose logging (`SolveOptions::verbose = false`) unless debugging
  convergence, as streaming residuals slows execution considerably.
* Run `solver_benchmark` (see Section 4) before committing new solver tweaks to
  record performance numbers for regressions.

## 4. Solver benchmarking toolkit

`build/solver_benchmark` exercises the same setup as the analytic regression but
allows custom grid resolutions, tolerances, and repetition counts. Example:

```bash
./build/solver_benchmark --nx 129 --ny 129 --max-iters 5000 --tol 1e-6 --repeats 5
```

The tool reports:

* convergence outcome and residual,
* average/min/max wall-clock time,
* throughput expressed as millions of cell-iterations per second,
* milliseconds per million cell-iterations (handy for quick mental estimates).

Add `--csv outputs/solver_benchmarks.csv` to append results for later plotting or
diffing. The CSV header is written automatically the first time.

## 5. Updating documentation with new metrics

When solver changes land, re-run the benchmark on representative grids (e.g.,
129² and 257²) and update `docs/solver_performance.md` with the new numbers.
Include the machine profile (Codespaces, local workstation, etc.) so future
contributors can compare apples to apples.

## 6. Visualisation helpers

The regression test and scenario runner emit optional CSV artefacts under
`outputs/` when they finish. Three Python scripts consume them:

* `python/visualize_wire.py` compares the centreline magnitude against the
  analytic solution for a finite-radius wire. The dip to zero at the origin is a
  physical consequence of the uniform current distribution inside the conductor.
* `python/visualize_wire_field.py` renders a top-down map with quiver arrows so
  you can verify the azimuthal direction of **B** around the wire.
* `python/examples/two_wire_cancel.py` writes `inputs/two_wire_cancel.json` using
  the lightweight scenario DSL. Regenerate it if you tweak parameters during
  debugging.

## 7. Other IDEs

CLion, Qt Creator, or bare-terminal workflows operate via the same CMake entry
points. The only project-specific requirement is that headers under `include/`
remain on the compiler include path. The `scripts/setup_env.sh` helper installs
the minimal toolchain on Debian/Ubuntu-based systems if you need a quick start.


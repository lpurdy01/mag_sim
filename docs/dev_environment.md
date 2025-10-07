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

## 2. VS Code configuration

The repo ships with `.vscode/` settings tuned for the `CMake Tools` and `C/C++`
extensions.

### 2.1 Tasks

Use `Terminal → Run Task` (or the command palette) to access the predefined
tasks:

| Task label              | Purpose                                                                  |
| ----------------------- | ------------------------------------------------------------------------ |
| `cmake-configure`       | Configure the out-of-source build directory in Debug mode.               |
| `cmake-build`           | Build all targets; depends on configuration.                             |
| `ctest`                 | Run the full test suite with verbose failure output.                     |
| `run-analytic-wire-test`| Compile + launch the analytic validation binary in place.                |
| `run-solver-benchmark`  | Compile + execute the benchmark helper with a lean 129×129 grid preset. |

These tasks chain automatically (for example, running the benchmark will trigger
a build if required). Feel free to duplicate them locally for custom grids or
release builds.

### 2.2 Debugging

The `Run and Debug` panel exposes launch configurations for three executables:

1. **motor_sim** — the placeholder main application.
2. **analytic_wire_test** — deterministic regression for the single-wire case.
3. **solver_benchmark** — convenience harness for measuring throughput.

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

## 6. Other IDEs

CLion, Qt Creator, or bare-terminal workflows operate via the same CMake entry
points. The only project-specific requirement is that headers under `include/`
remain on the compiler include path. The `scripts/setup_env.sh` helper installs
the minimal toolchain on Debian/Ubuntu-based systems if you need a quick start.


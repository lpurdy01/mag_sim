# Command-Line Interface

The `motor_sim` executable accepts the following high-level flags:

| Flag | Description |
| ---- | ----------- |
| `--scenario PATH` | Load a scenario JSON file. Required for all solves. |
| `--solve` | Execute the magnetostatic solve after rasterising the scenario. |
| `--solver {sor,cg}` | Select the iterative solver (default: `cg`). |
| `--tol VALUE` | Relative residual tolerance (default: `1e-6`). |
| `--max-iters N` | Iteration limit for the selected solver. |
| `--warm-start` | Reuse the previous frame’s solution as the initial guess. |
| `--use-prolongation` | Seed a fine grid from a coarse solve (pair with `--coarse-nx/--coarse-ny`). |
| `--vtk-series PATH` | Emit a ParaView `.pvd` time-series descriptor alongside per-frame `.vti` outputs. |
| `--outputs LIST` | Restrict emitted outputs by ID (`--outputs none` disables emission). |
| `--list-outputs` | Print the outputs declared by the scenario without solving. |
| `--parallel-frames` | Solve timeline frames concurrently (up to hardware concurrency minus one). |
| `--progress-every SECONDS` | Adjust progress reporting cadence. |
| `--quiet` | Suppress progress output. |
| `--snapshot-every N` | Request downsampled diagnostics every `N` iterations. |

Run `./build/motor_sim --help` to inspect the complete flag list, including regression and development switches.

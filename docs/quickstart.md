# Quickstart

Follow these steps to build the solver, generate a bundled scenario, and inspect the resulting fields.

## 1. Build the project

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

## 2. Generate a scenario

The repository ships a generator for the three-phase stator demo:

```bash
python3 python/gen_three_phase_stator.py --profile ci --out inputs/three_phase_stator_ci.json
```

## 3. Solve the field

```bash
./build/motor_sim --scenario inputs/three_phase_stator_ci.json --solve --vtk-series outputs/three_phase_ci.pvd
```

Use `--parallel-frames` to solve timeline frames concurrently and `--warm-start` to reuse the previous solution when stepping through the series.

## 4. Visualize the results

Open the emitted `.pvd` time series in ParaView or your preferred VTK viewer. See [Visualization Overview](user-guide/visualization/overview.md) for guidance and the [VTK Output reference](reference/file-formats/vtk.md) for file structure details.

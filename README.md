# 2D Electromagnetic Motor Simulator (Bootstrap)

This is the initial scaffold: CMake + C++17 + VS Code + Devcontainer + CI.

## Build (WSL/Linux)

```bash
scripts/setup_env.sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build -j
./build/motor_sim
```

## Test

```bash
cd build && ctest --output-on-failure
```

## Next

* Add field data structures and a simple iterative solver.
* Add CSV/VTK writers and Python notebooks for plots.

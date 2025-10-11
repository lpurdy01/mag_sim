#!/usr/bin/env bash

# Run a close approximation of the GitHub Actions CI locally so changes to the
# workflow can be validated before pushing.

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${BUILD_DIR:-${ROOT_DIR}/build}"
CTEST_PARALLEL_LEVEL="${CTEST_PARALLEL_LEVEL:-}" 
if [[ -z "$CTEST_PARALLEL_LEVEL" ]]; then
  if command -v nproc >/dev/null 2>&1; then
    CTEST_PARALLEL_LEVEL="$(nproc)"
  else
    CTEST_PARALLEL_LEVEL="1"
  fi
fi

cd "$ROOT_DIR"

echo "[ci-check] Configuring project (build dir: $BUILD_DIR)"
cmake -S . -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE="${CMAKE_BUILD_TYPE:-Release}"

echo "[ci-check] Building"
cmake --build "$BUILD_DIR" -j

echo "[ci-check] Running full test suite"
(cd "$BUILD_DIR" && ctest --output-on-failure --parallel "$CTEST_PARALLEL_LEVEL")

echo "[ci-check] Running analytic wire regression"
(cd "$BUILD_DIR" && ctest --output-on-failure -R analytic_wire)

if [[ "${SKIP_PY_DEPS:-0}" != "1" ]]; then
  echo "[ci-check] Ensuring visualization Python dependencies are available"
  python3 -m pip install --user --quiet numpy matplotlib vtk
fi

echo "[ci-check] Preparing artifact directories"
rm -rf outputs ci_artifacts
mkdir -p outputs ci_artifacts

declare -a regression_bins=(
  analytic_wire_test
  two_wire_cancel_test
  line_current_interface_test
  magnet_strip_test
  probe_output_test
  torque_validation_test
  back_emf_probe_test
  rotor_ripple_test
  warm_start_prolongation_test
  progress_snapshot_test
)

for bin in "${regression_bins[@]}"; do
  echo "[ci-check] Running $bin"
  "$BUILD_DIR/$bin" 2>&1 | tee "ci_artifacts/${bin}.txt"
done

{
  echo "# Regression accuracy summary"
  echo
  for bin in "${regression_bins[@]}"; do
    report="ci_artifacts/${bin}.txt"
    name="${bin}"
    echo "## ${name}"
    cat "$report"
    echo
  done
} > ci_artifacts/test_accuracy_report.txt

echo "[ci-check] Solving scenarios"
"$BUILD_DIR/motor_sim" --scenario inputs/two_wire_cancel.json --solve --outputs all --solver sor
"$BUILD_DIR/motor_sim" --scenario inputs/two_wire_cancel.json --solve --outputs none --solver cg --quiet
"$BUILD_DIR/motor_sim" --scenario inputs/line_current_interface.json --solve --outputs all --solver sor
"$BUILD_DIR/motor_sim" --scenario inputs/line_current_interface.json --solve --outputs none --solver cg --quiet
"$BUILD_DIR/motor_sim" --scenario inputs/iron_ring_demo.json --solve --outputs all --tol 5e-6 --max-iters 40000 --solver sor
"$BUILD_DIR/motor_sim" --scenario inputs/tests/magnet_strip_test.json --solve --outputs all --solver sor
"$BUILD_DIR/motor_sim" --scenario inputs/tests/magnet_strip_test.json --solve --outputs none --solver cg --quiet
"$BUILD_DIR/motor_sim" --scenario inputs/tests/rotor_ripple_test.json --solve --outputs all --parallel-frames --max-iters 40000 --tol 2.5e-6 --solver sor
"$BUILD_DIR/motor_sim" --scenario inputs/tests/rotor_ripple_test.json --solve --outputs none --parallel-frames --max-iters 40000 --tol 2.5e-6 --solver cg --quiet --warm-start

echo "[ci-check] Verifying VTK export"
python3 python/verify_vtk.py outputs/rotor_ripple_field_frame_000.vti

echo "[ci-check] Rendering visualisations"
export MPLBACKEND=Agg
python3 python/visualize_scenario_field.py \
  --scenario inputs/two_wire_cancel.json \
  --field-map outputs/two_wire_field_map.csv \
  --save ci_artifacts/two_wire_field.png \
  --draw-boundaries \
  --streamlines
python3 python/visualize_scenario_field.py \
  --scenario inputs/line_current_interface.json \
  --field-map outputs/interface_field.csv \
  --save ci_artifacts/line_current_interface_field.png \
  --draw-boundaries \
  --streamlines \
  --overlay-analytic interface \
  --analytic-contours 10
python3 python/visualize_scenario_field.py \
  --scenario inputs/iron_ring_demo.json \
  --field-map outputs/iron_ring_field.csv \
  --save ci_artifacts/iron_ring_field.png \
  --draw-boundaries \
  --streamlines \
  --color-scale log \
  --vector-mode log
python3 python/visualize_scenario_field.py \
  --scenario inputs/tests/magnet_strip_test.json \
  --field-map outputs/magnet_strip_field.csv \
  --save ci_artifacts/magnet_strip_field.png \
  --draw-boundaries \
  --streamlines \
  --color-scale log \
  --vector-mode log
python3 python/visualize_scenario_field.py \
  --scenario inputs/tests/rotor_ripple_test.json \
  --field-map outputs/rotor_ripple_field_frame_000.csv \
  --save ci_artifacts/rotor_ripple_frame0.png \
  --outline-vtp outputs/rotor_ripple_field_frame_000_outlines.vtp \
  --draw-boundaries \
  --streamlines \
  --color-scale log \
  --vector-mode log
python3 python/visualize_scenario_field.py \
  --scenario inputs/tests/rotor_ripple_test.json \
  --field-map outputs/rotor_ripple_field_frame_002.csv \
  --save ci_artifacts/rotor_ripple_frame2.png \
  --outline-vtp outputs/rotor_ripple_field_frame_002_outlines.vtp \
  --draw-boundaries \
  --streamlines \
  --color-scale log \
  --vector-mode log
python3 python/visualize_wire.py \
  --csv outputs/validation_wire_line.csv \
  --save ci_artifacts/analytic_wire_line.png \
  --no-show

echo "[ci-check] Collecting CSV artifacts"
shopt -s nullglob
for file in outputs/*.csv outputs/*.vti outputs/*_outlines.vtp outputs/*_outlines_labels.csv; do
  cp "$file" ci_artifacts/
done

echo "[ci-check] Done. Artefacts are available in ci_artifacts/"

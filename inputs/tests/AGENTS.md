# Test scenarios

- Every regression test in `tests/` should load its scenario from this folder to
  exercise the JSON ingestion path. When you update a test expectation, update
  the corresponding JSON here in the same commit.
- Keep filenames aligned with the test binary (e.g.,
  `magnet_strip_test.cpp` ↔ `magnet_strip_test.json`) so CI scripts can discover
  them mechanically.
- Scenarios under `inputs/tests/` are referenced by the CI artifact workflow;
  avoid renaming IDs or outputs without updating `.github/workflows/ci.yml` and
  `scripts/run_ci_checks.sh` accordingly.

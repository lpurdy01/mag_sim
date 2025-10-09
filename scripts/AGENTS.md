# Scripts directory tips

- Keep automation helpers POSIX shell compatible unless there is a compelling
  reason to depend on bash-specific features. New scripts should start with
  `#!/usr/bin/env bash` and `set -euo pipefail` if they require bash.
- `run_ci_checks.sh` mirrors `.github/workflows/ci.yml`; update them together so
  contributors can catch regressions locally before pushing.
- Do not add Python virtual environment directories or generated artefacts
  (e.g., `ci_artifacts/`) under `scripts/`—they belong at the repository root
  and should remain gitignored.

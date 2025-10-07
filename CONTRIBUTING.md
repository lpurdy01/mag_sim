# Contributor Onboarding Quickstart

1. **Clone & enter** the repository: `git clone https://github.com/<your-org>/mag_sim && cd mag_sim`.
2. **Review the charter** below to understand goals, physics scope, and expectations.
3. **Set up the toolchain** using `scripts/setup_env.sh` (WSL/Codespaces) or install GCC, CMake ≥ 3.16, and clang-format manually.
4. **Configure the build**: `cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug`.
5. **Build the simulator**: `cmake --build build -j`.
6. **Run smoke tests**: `ctest --test-dir build` (or execute specific tests under `tests/`).
7. **Explore docs** in `docs/` for physics references and implementation notes.
8. **Create a branch** using the naming conventions in Section 8 (e.g., `docs/contributing`).
9. **Implement changes** adhering to style, documentation, and architectural rules.
10. **Verify the QA checklist** before opening your pull request.

---

# Contributing & Development Standards
*(2D Electromagnetic Motor Simulator — C/C++ Learning Project)*

---

## 1. Guiding Principles

This project exists primarily as a **learning and experimentation platform** for 2D magnetostatic simulation, not as production software.  
All decisions—technical or organizational—should prioritize:

- **Transparency over abstraction**: favor explicit, well-commented code over opaque libraries.  
- **Educational clarity**: every function should help reveal how a simulator works.  
- **Reproducibility**: anyone (human or AI agent) should be able to clone, build, and run within minutes.  
- **Incremental iteration**: commit working slices of functionality—don’t over-engineer ahead of learning needs.  

---

## 2. Repository Structure & Scope

The repository should maintain the following structure and purpose:

```

src/             → Core C++ sources (simulation, solver, IO)
include/         → Public headers for importable modules
tests/           → Minimal regression & smoke tests
python/          → Optional post-processing & visualization scripts
scripts/         → Setup & automation (e.g., setup_env.sh)
.vscode/         → IDE integration tasks, launch, settings
.devcontainer/   → Codespaces & WSL reproducible environment
.github/workflows/ → CI definitions (Ubuntu runner)
docs/            → Theory notes, equations, and developer documentation

```

- Keep **each directory independent and self-documenting** with README.md files as they mature.  
- Avoid monolithic files; instead, isolate concerns (e.g., geometry parsing, solver core, I/O writer).  

---

## 3. Language & Toolchain

- **Language standard:** C++17 (minimum).  
- **Compiler:** GCC / G++ on Ubuntu 22.04 (default via WSL or Codespaces).  
- **Build system:** CMake ≥ 3.16, out-of-source builds only (`build/`).  
- **IDE:** VS Code + Remote-WSL (C/C++ + CMake Tools extensions).  
- **Debugger:** GDB integrated via VS Code.  
- **Formatting:** `clang-format` (LLVM/Allman, 4-space indents).  
- **Static analysis:** enable `-Wall -Wextra -Wpedantic`; treat warnings as errors in CI builds.

---

## 4. Coding Conventions

### 4.1 File Organization
Each `.hpp` and `.cpp` file must start with:
```cpp
// filename: <name>.cpp
// part of 2D Electromagnetic Motor Simulator
```

### 4.2 Naming & Style

* **Classes**: `PascalCase` (e.g., `Grid2D`, `FieldSolver`)
* **Functions & variables**: `camelCase`
* **Constants/macros**: `ALL_CAPS`
* **Namespaces**: always under `motorsim::`

### 4.3 Documentation

Use **Sphinx-style docstrings** for all public functions and classes:

```cpp
/**
 * @brief Compute B-field from A-field.
 * @param Az Grid of vector potential values.
 * @return Vector field of Bx, By.
 */
```

Add **inline math or references** for physics formulas where relevant.

### 4.4 Includes

* Use **forward declarations** where possible.
* Include headers only where necessary.
* Avoid circular dependencies (see LabEx best practices).

### 4.5 Units & Constants

* Work exclusively in **SI units**: meters, amperes, tesla, henry.
* Keep `MU0` and other constants defined in `types.hpp`.

---

## 5. Simulation Architecture Rules

1. **Single Responsibility** — each component should do one thing:

   * `Grid2D`: domain discretization
   * `Solver`: field computation
   * `Writer`: output formats (CSV/VTK)
   * `Main`: orchestration

2. **Transparency over performance** — implement math explicitly before optimizing.

3. **Minimal dependencies** — standard library only, unless absolutely necessary.

   * Acceptable optional: header-only JSON or Eigen (with justification).
   * No GUI or visualization libraries inside the C++ core.

4. **Separation of Concerns** — physics logic (solver) must not depend on file I/O or visualization.

5. **Deterministic output** — same inputs → same results. Randomness discouraged unless seeded.

---

## 6. Output and Visualization

* **Preferred outputs:**

  * `.vtk` structured/unstructured grid (ASCII for small data)
  * `.csv` for quick analysis (x, y, Bx, By)

* **Visualization tools:**

  * Python scripts (`numpy`, `matplotlib`) in `python/`
  * ParaView for vector-field inspection

* Ensure all output writers include **units and field metadata** in headers.

---

## 7. Testing & Validation

### 7.1 Scope

* Keep at least one **smoke test** in `tests/` verifying build integrity.
* Add **unit tests** for math routines (e.g., Laplacian stencil, boundary handling).

### 7.2 Expected Physics Checks

* Single straight wire → circular field (`B ~ μ₀I / 2πr`)
* Two opposite wires → cancellation at midpoint
* Iron region → flux concentration

These can be visual or numerical comparisons to analytical expectations.

---

## 8. Commit & Branching Policy

* **Branch naming:**

  * Feature: `feat/<topic>`
  * Fix: `fix/<issue>`
  * Documentation: `docs/<section>`
  * Infra: `chore/<task>`

* **Commit structure:**

  ```
  <type>: <short description>

  [optional body explaining rationale]
  ```

  Example:

  ```
  feat: implement Gauss-Seidel solver core

  Adds first working version of iterative solver for Az(x,y).
  ```

* Always run `ctest` or smoke test before pushing.

---

## 9. Continuous Integration (CI)

* GitHub Actions workflow (`.github/workflows/ci.yml`) must:

  * Build with CMake in Debug and Release
  * Run smoke/unit tests
  * Enforce compilation warnings
  * Validate presence of `README.md` and licensing headers

---

## 10. Environment Setup & Reproducibility

* **WSL/Codespaces setup script:** `scripts/setup_env.sh` installs all dependencies.

* **Devcontainer**: defined in `.devcontainer/` for identical builds.

* **Build flow:**

  ```bash
  cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
  cmake --build build -j
  ./build/motor_sim
  ```

* Document new dependencies or environment changes in `README.md`.

---

## 11. Documentation & Learning Notes

* Maintain physics background, derivations, and theory references in `docs/`.
* Cite relevant equations and references (e.g., ANSYS Maxwell theory for magnetostatics).
* Add illustrative figures or derivations using Markdown or embedded LaTeX.

---

## 12. Use of AI Assistants

* AI-generated code must:

  * Adhere to this style guide.
  * Include math references or derivation notes for physics routines.
  * Pass local and CI smoke tests before merge.

When using agents (Codex, GitHub Copilot, etc.), ensure changes remain human-readable and educational.

---

## 13. Licensing & Attribution

* The project is open-source under the **MIT License**.
* Include the license header in all new source files.
* Cite external content (theory, equations, or code snippets) in comments or docs.

---

## 14. Extension & Integration Rules

* **Simulink Integration:** generated code must compile as standalone C (no MEX or MATLAB runtime).

  * Place under `external/simulink/` with separate CMake target.
  * Keep interfaces minimal (`Model_initialize()`, `Model_step()`).

* **Future modules (motion, eddy currents)**: isolate in new namespaces or subfolders.

  * Example: `motorsim::dynamic`, `motorsim::nonlinear`.

---

## 15. Quality Assurance Checklist

Before opening a Pull Request:

* [ ] Code compiles without warnings.
* [ ] Smoke test passes.
* [ ] README updated if build/run instructions changed.
* [ ] Added docstrings and physics comments.
* [ ] CI workflow passes on GitHub.
* [ ] Commit messages are clean and descriptive.

---

## 16. Philosophy Recap

> “Readable, reproducible, and physically accurate before fast.”

Every contribution should leave the repository:

* **More understandable** to a new learner,
* **More reproducible** for automation or CI,
* **Closer** to a complete understanding of how electromagnetic solvers work.

---

*Document version 1.1 — aligned with the Project Charter for the 2D Electromagnetic Motor Simulator (C/C++ Learning Project).*


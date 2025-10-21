# Scenario JSON

Scenario files define geometry, materials, sources, and outputs for the solver. This stub highlights the most frequently used fields; consult the schema in `python/scenario_api.py` for exhaustive options.

| Field | Type | Description |
| ----- | ---- | ----------- |
| `domain` | object | Simulation extents, grid resolution, and boundary conditions. |
| `materials` | array | List of material definitions referenced by regions. |
| `regions` | array | Geometry definitions (`polygon`, `halfspace`, etc.) with material assignments. |
| `sources` | array | Current-carrying conductors with winding metadata. |
| `magnets` | array | Permanent magnet regions with magnetisation vectors. |
| `timeline` | array | Optional per-frame overrides for currents, magnet orientation, and rotor angles. |
| `outputs` | object | Requested exports such as field maps, line probes, or VTK series. |
| `rotors` | array | Rotor assemblies that group regions, magnets, and wires for rigid motion. |

Additional helper fields (e.g., `phase_currents`, `transient`) are documented in the developer guide. Future revisions will expand this page with concrete JSON fragments and validation rules as the schema stabilises.

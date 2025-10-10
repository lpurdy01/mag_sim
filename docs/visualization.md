# Field Visualisation Options

The `python/visualize_scenario_field.py` helper renders solver field maps with a
number of switches that mirror the "Visualization Upgrades" milestone from the
roadmap.

## Command-line summary

```
python3 python/visualize_scenario_field.py \
  --scenario inputs/two_wire_cancel.json \
  --field-map outputs/two_wire_field_map.csv \
  --draw-boundaries \
  --streamlines \
  --overlay-analytic interface \
  --vector-mode log \
  --vector-log-floor 1e-6 \
  --color-scale log \
  --log-floor 1e-6
```

Key options:

- `--color-scale {linear,log}` toggles between a linear palette and logarithmic
  magnitude scaling. Pair with `--log-floor` to bound the minimum |B| value.
- `--draw-boundaries` outlines material and magnet regions extracted from the
  scenario JSON so discontinuities are easy to spot.
- `--streamlines` overlays streamlines traced from the field grid.
- `--vector-mode {linear,log,off}` controls the quiver arrows. The `log` mode
  compresses dynamic range using a log-normalised magnitude, while `off`
  disables the overlay entirely. `--vector-log-floor` mirrors `--log-floor` for
  the quiver transform.
- `--overlay-analytic interface` draws contour lines for the planar permeability
  interface analytic reference when available. Additional overlays can be added
  under the same flag in the future.
- `--save PATH` writes the rendered figure to disk; otherwise the viewer window
  opens interactively.

## CI integration

`scripts/run_ci_checks.sh` exercises these switches for the canonical demo
scenarios so GitHub Actions artifacts include the enhanced renders. The iron ring
and magnet strip plots use `--color-scale log` together with
`--vector-mode log` to improve readability when the dynamic range spans several
orders of magnitude.

# SOR vs CG Visual Comparison

These artefacts highlight the difference between the historical SOR solver and the modern preconditioned CG backend.

![CG vs SOR convergence overlay](../../assets/images/CG_vs_SOR_convergence_overlay.png)

CG reaches the target tolerance in roughly an order of magnitude fewer iterations on the iron-ring plus magnet scenario. Both runs use identical tolerances and grids.

## Field evolution snapshots

![CG progress snapshots](../../assets/images/progress_snapshots_cg.gif)

![SOR progress snapshots](../../assets/images/progress_snapshots_sor.gif)

The snapshots expose how CG removes low-frequency error earlier in the solve, yielding smooth fields while SOR continues to chew through residual ridges. Use these animations to explain solver choices to stakeholders or when preparing presentations.

## Warm-started frame

![Warm-started frame 2](../../assets/images/warm_start_frame_2.png)

Warm starts keep consecutive timeline frames aligned, especially when combined with `--parallel-frames`. The second frame in the three-phase stator demo converges substantially faster than a cold start when warm-starting is enabled.

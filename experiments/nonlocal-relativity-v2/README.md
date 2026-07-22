# Nonlocal Relativity v2 — reproducible experiment suite

Deterministic numerical experiments for the experimental
`scirust-nonlocal-relativity` layer (fractional-memory test-particle worldline
dynamics on fixed general-relativistic backgrounds). Every experiment is a
pure-Rust binary with **no RNG and no wall-clock dependence**, so identical
inputs produce byte-identical output. Each prints a `#`-prefixed metadata and
units header, then CSV rows, validates that every emitted number is finite, and
closes with a short interpretation.

**These are numerical experiments on a fixed phenomenological model. None of
them is a physical validation, and none of them establishes new physics.** See
[`docs/EXPERIMENTAL_NONLOCAL_RELATIVITY.md`](../../docs/EXPERIMENTAL_NONLOCAL_RELATIVITY.md)
for the scientific boundary.

## Running

```bash
# Optional provenance stamp (metadata only; does not affect numeric output):
export NLR_EXPERIMENT_COMMIT="$(git rev-parse HEAD)"

cargo run --release -p nonlocal-relativity-experiments --bin history_retention
```

Output goes to stdout as `#`-commented metadata plus CSV. Redirect it to a file
if you want to keep it; generated CSV and compiled binaries are intentionally
**not** committed (only the deterministic source is).

## Experiments

### `history_retention` — Phase 3 history-retention comparison

Compares the two persistent-history retention strategies of the step-doubling
adaptive integrator against an independent fine fixed-step reference:

- `EndpointOnly` (default): retain only each accepted step's endpoint.
- `RefinedAcceptedHistory`: additionally retain each accepted step's midpoint,
  at its true affine parameter.

Columns: `tolerance, strategy, accepted_steps, retained_samples,
op_count_proxy, endpoint_coord_err, endpoint_vel_err, memory_l2,
memory_force_l2, metric_norm_drift`. The reference is a fixed-step
semi-implicit-Euler run at `h = 5e-4` with the same non-uniform Caputo memory
law — a numerical reference, not an exact solution.

**Measured result and decision.** Across tolerances `1e-6 … 1e-9`, retaining
the midpoints leaves the accepted-step count unchanged and the endpoint
coordinate/velocity error identical to ~4 significant figures, while roughly
doubling the retained sample count and the operation-count proxy (the total
`O(N^2)` Caputo work). Representative rows (`commit`-independent):

| tolerance | strategy | accepted_steps | retained | op_proxy | coord_err |
|-----------|----------|----------------|----------|----------|-----------|
| 1e-8 | endpoint_only     | 40  | 41  | 861   | 1.197e-6 |
| 1e-8 | refined_accepted  | 40  | 81  | 1681  | 1.197e-6 |
| 1e-9 | endpoint_only     | 124 | 125 | 7875  | 3.387e-7 |
| 1e-9 | refined_accepted  | 124 | 249 | 15625 | 3.387e-7 |

The initial hypothesis — that a denser retained history would improve accuracy
for this hereditary equation — is **not supported** on this experiment: the
memory force is a small perturbation on the geodesic dynamics, and the endpoint
accuracy is dominated by the first-order integrator's truncation error, not by
the memory-quadrature density. Both strategies converge under tolerance
refinement (coord error ~`2.9e-6 → 3.4e-7` as tolerance tightens `1e-6 →
1e-9`).

**Decision: keep `EndpointOnly` as the default**; expose
`RefinedAcceptedHistory` only as an explicit research option
(`simulate_nonlocal_worldline_adaptive_with_stepper_policy_retention` with
`HistoryRetention::RefinedAcceptedHistory`). The structural invariants of both
strategies (strict parameter ordering, no duplicate parameters, true-midpoint
recording, exact retained-sample counts with no leakage from rejected trials)
are pinned by `scirust-nonlocal-relativity/tests/history_retention.rs`.

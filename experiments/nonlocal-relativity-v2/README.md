# Nonlocal Relativity v2 — reproducible experiment suite

Deterministic numerical experiments for the SciRust relativity stack. The
suite spans two scientific categories, kept strictly separate:

- **Experimental phenomenological layer** (`scirust-nonlocal-relativity`):
  fractional-memory test-particle worldline dynamics on fixed
  general-relativistic backgrounds (the [Experiments](#experiments) section).
- **Established general relativity** (`scirust-relativity` geometry core and
  the Layer 2 Covariant Gravity Workbench): textbook-GR primitives checked
  against exact analytic and closed-form oracles (the
  [Established general-relativity experiments](#established-general-relativity-experiments-geometry-core-and-layer-2)
  section).

Every experiment is a pure-Rust binary with **no RNG and no wall-clock
dependence**, so identical inputs produce byte-identical output. Each prints a
`#`-prefixed metadata and units header, then CSV rows, validates that every
emitted number is finite, and closes with a short interpretation.

**The phenomenological experiments are numerical experiments on a fixed
phenomenological model: none is a physical validation, and none establishes new
physics.** See
[`docs/EXPERIMENTAL_NONLOCAL_RELATIVITY.md`](../../docs/EXPERIMENTAL_NONLOCAL_RELATIVITY.md)
for the scientific boundary. **The established-GR experiments validate the
implementation against known results of textbook general relativity; they too
introduce no new physics.**

## Running

```bash
# Optional provenance stamp (metadata only; does not affect numeric output):
export NLR_EXPERIMENT_COMMIT="$(git rev-parse HEAD)"

cargo run --release -p nonlocal-relativity-experiments --bin history_retention
cargo run --release -p nonlocal-relativity-experiments --bin adaptive_convergence
cargo run --release -p nonlocal-relativity-experiments --bin complexity_scaling
cargo run --release -p nonlocal-relativity-experiments --bin bounded_memory_error
cargo run --release -p nonlocal-relativity-experiments --bin kerr_fd_sensitivity
cargo run --release -p nonlocal-relativity-experiments --bin modulation_sensitivity
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

### `adaptive_convergence` — Phase 9 adaptive-tolerance convergence

For both adaptive controllers (embedded Heun-Euler and step-doubling), reports
the accepted-step count and endpoint coordinate error against a fine fixed-step
reference of the matching method, at tolerances `1e-5 … 1e-9`, with the
error-reduction ratio between consecutive rows. The endpoint error decreases
monotonically as the tolerance tightens (embedded Heun-Euler reaches ~`3.5e-10`,
step-doubling ~`3.4e-7` at `1e-9`). This is numerical self-consistency toward a
fine-grid reference, **not** a validation of the model.

### `complexity_scaling` — Phase 10 empirical complexity (op-count proxy)

Measures the deterministic operation-count proxy (sum of retained sample counts
over all accepted evaluations — the total Caputo history-sample touches) for the
complete raw coordinate memory, the bounded short memory (`W = 16`), and the
discrete connection transport, at doubling fixed step counts `N ∈ {50, 100, 200,
400}`. The measured growth ratio `proxy(2N)/proxy(N)` matches the
implementation-derived complexity: `→ 4` (`O(N^2)`) for complete and discrete
transport, `→ 2` (`O(N*W)`) for bounded short memory. Discrete transport shares
complete memory's `O(N^2)` touch count but pays an extra `O(D^3)` Christoffel
contraction per touch. Scaling only — no wall-clock claim.

### `bounded_memory_error` — Phase 9 short-memory approximation error

Quantifies the accuracy cost of `BoundedShortMemoryHistory` (the `O(N*W)`
short-memory approximation) against the complete-history oracle (`O(N^2)`), on a
fixed-step Schwarzschild trajectory of 128 steps. For windows `W ∈ {4, 8, 16,
32, 64, 129}` it reports the endpoint coordinate/velocity error against the
oracle and the retained sample count. The error decreases monotonically as `W`
grows (`~3.7e-6` at `W=4` to `~2.5e-7` at `W=64`) and is **exactly zero** once
`W` covers every sample (the window then *is* the full history — a bit-for-bit
match, pinned by `tests/bounded_memory.rs`). This is the truncation cost the
user opts into by choosing a bounded backend, not a model validation.

### `kerr_fd_sensitivity` — Phase 9 Kerr finite-difference sensitivity

Unlike the analytic backgrounds, `Kerr`'s connection is evaluated by central
finite differences. At spin `a = 0` the Kerr metric reduces to Schwarzschild
exactly, so the finite-difference Christoffel symbols can be compared against
Schwarzschild's **exact analytic** ones. Sweeping the difference step exposes
the classic central-difference V-curve: the maximum Christoffel-component error
falls as `~O(h^2)` (`5.8e-8` at `h=1e-2` to `2.7e-11` near `h=1e-4`) then rises
as floating-point cancellation dominates (`3.2e-9` at `h=1e-6`). This quantifies
the disclosed truncation cost of the Kerr connection; every other background
uses exact analytic symbols.

### `modulation_sensitivity` — Phase 9 modulation sensitivity, β=0 baseline

Sweeps `beta` for `SchwarzschildKretschmannModulator` (weight `q = 1 + β·L⁴·K`)
and reports the endpoint deviation from the unmodulated baseline under adaptive
stepping. At `beta = 0` the deviation is **exactly zero and bit-identical** to
the baseline (the modulator's `beta = 0` bypass), and the deviation grows
monotonically (here ~linearly) with `beta` (`4.1e-11` at `β=0.1` to `8.3e-10` at
`β=2.0`). `beta` and the reference length are free, uncalibrated phenomenological
parameters — this quantifies a hook's sensitivity, not a physical effect. (The
`beta = 0` bit-identity is separately pinned by
`adaptive_beta_zero_modulator_matches_identity_modulator_bit_for_bit` in
`tests/adaptive.rs`.)

### Other scenarios (shipped as crate examples)

The remaining v2-plan scenarios already ship as deterministic, CSV-producing
examples of `scirust-nonlocal-relativity` and are not duplicated here:
coordinate-chart comparison (`coordinate_covariance`), exact flat-spacetime
transport convergence (`exact_transport_convergence`), exact Schwarzschild
circular-orbit transport convergence (`schwarzschild_orbit_transport`), and
proper-time vs affine memory (`proper_time_memory_comparison`). Run any with
`cargo run --release -p scirust-nonlocal-relativity --example <name>`.

## Established general-relativity experiments (geometry core and Layer 2)

The same binary crate also hosts the deterministic experiments that validate the
`scirust-relativity` geometry core and the Layer 2 Covariant Gravity Workbench
against **exact analytic and closed-form oracles of textbook general
relativity**. They share the reproducibility contract above (no RNG, no
wall-clock, finite-checked CSV), but belong to a different scientific category:
each is an implementation validation against established GR, not a study of the
phenomenological memory model. Run any with

```bash
cargo run --release -p nonlocal-relativity-experiments --bin <name>
```

Geometry-core primitives:

- `orthonormal_tetrad` — local orthonormal-frame (tetrad) construction:
  orthonormality, completeness, and the temporal/spatial split against the
  closed-form metric projection. Exact single-point construction, so residuals
  sit at the rounding floor.
- `parallel_transport` — parallel-transport metric-compatibility drift, flat
  closed-loop holonomy, and the holonomy/curvature identity around a small
  parallelogram (cross-checking the transport and curvature engines).
- `covariant_transport` — parallel transport of covectors and tensors under
  metric compatibility (`nabla g = 0`), including metric self-transport.
- `curvature_invariants` — Riemann/Ricci/Einstein/Kretschmann tensors against
  the exact Schwarzschild and (anti-)de Sitter scalar oracles, with a
  finite-difference step sweep exposing the central-difference trade-off.
- `coordinate_independence` — chart-independence of the Ricci and Kretschmann
  scalars (Cartesian vs spherical Minkowski; areal vs isotropic Schwarzschild,
  matched through the areal radius).
- `flrw_curvature` — spatially flat FLRW curvature against the exact Friedmann
  formulas (exponential/de Sitter and power-law scale factors).
- `geodesic_deviation` — Jacobi-field tidal focusing: exact linear growth in
  flat spacetime, (de)focusing in (anti-)de Sitter, and the Schwarzschild
  radial-vs-transverse tidal asymmetry.
- `exponential_map` — geodesic exponential/logarithm round-trip accuracy
  `|log_p(exp_p(v)) - v|` across backgrounds.
- `world_function` — Synge world function: flat exactness plus the
  convention-free symmetry, fundamental (`2 sigma = g sigma^mu sigma_mu`), and
  gradient-round-trip identities.
- `van_vleck_determinant` — van Vleck–Morette determinant from the
  exponential-map Jacobian: flat exactness, symmetry, and the near-coincidence
  expansion.

Layer 2 — Covariant Gravity Workbench:

- `linearized_gravity` — the weak-field Einstein equations and their four
  oracles (Newtonian Poisson limit, weak-field Schwarzschild vacuum, gauge
  invariance, and the `O(h^2)` nonlinear cross-check); see
  [`docs/LAYER_2_COVARIANT_GRAVITY.md`](../../docs/LAYER_2_COVARIANT_GRAVITY.md).
- `ppn_extraction` — asymptotic extraction of the Eddington–Robertson PPN
  parameters `gamma` and `beta` from static isotropic weak-field metrics, with
  exact and contaminated synthetic oracles, isotropic-Schwarzschild convergence
  to `gamma = beta = 1`, and areal-coordinate rejection; see
  [`docs/LAYER_2_PPN.md`](../../docs/LAYER_2_PPN.md).

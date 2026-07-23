# Relativistic-Computation Platform — Architecture Report

*Phase 0 "understand before modifying" report. This is an internal engineering
map of the relativistic-computation subgraph as it exists today: what is built,
what is validated, what the numerical assumptions are, where the technical debt
and duplication live, and which research capabilities are still missing. It is
descriptive, not aspirational — the forward-looking vision lives in
[`PLATFORM_ROADMAP.md`](PLATFORM_ROADMAP.md); this document records the current
state so that future work is additive and non-duplicating.*

## 1. Scope

The wider repository is a ~100-crate monorepo (a pure-Rust deep-learning
framework, a transpiler, industrial/robotics stacks, and more). This report
covers **only** the relativistic-computation dependency subgraph:

```
scirust-special      (Gamma, erf, Bessel, incomplete gamma/beta, zeta, Loader pmfs)
        ▲
        │ (gamma only)
scirust-fractional   (Caputo L1 uniform/non-uniform, Grünwald–Letnikov / RL)
        ▲
        │
scirust-nonlocal-relativity  ──────────────┐  (hereditary worldline dynamics)
        │                                   │
        ▼                                   ▼
scirust-relativity   (geometry core)   scirust-fractional
        │
        ▼
scirust-sim          (deterministic ODE integration engine)

experiments/nonlocal-relativity-v2   (7 deterministic experiment binaries)
```

Edges are `path` dependencies. `scirust-relativity` depends only on
`scirust-sim` (default features, i.e. core RK4/adaptive only).
`scirust-nonlocal-relativity` depends on `scirust-relativity` +
`scirust-fractional`. The experiments crate depends on all three relativity
crates.

Every crate in this subgraph is `#![forbid(unsafe_code)]` and (except the
experiment/bin crates) `#![deny(missing_docs)]`. None uses `unsafe`, and none
contains `panic!`/`todo!`/`unimplemented!`/`unreachable!` on a runtime path;
`unwrap`/`expect` appear only in `#[cfg(test)]` modules and doctests. There is
no RNG or wall-clock in any library or experiment path (the only PRNG in the
subgraph, `scirust_sim::SplitMix64`, is used exclusively by unrelated stochastic
domain models, never by the integrators).

## 2. Current capabilities, per crate

### 2.1 `scirust-special` — special functions

Single-file crate, zero runtime dependencies. Flat public API. Every function
is a documented numerical implementation of established mathematics, each with a
literature reference and test oracles against `scipy.special` / `mpmath`.

- **Gamma family:** `ln_gamma` (Lanczos g=7, Godfrey coefficients, Euler
  reflection), `gamma`, `digamma`, `beta`, `ln_beta`.
- **Error function:** `erf`, `erfc` (via regularized incomplete gamma),
  `erfinv` (Giles 2010 seed + Halley steps).
- **Incomplete gamma / χ²:** `regularized_gamma_p` / `_q` (series + modified
  Lentz continued fraction; Temme-1987 boundary-layer iteration cap).
- **Incomplete beta:** `regularized_incomplete_beta` (Lentz CF).
- **Riemann zeta:** `riemann_zeta`, `riemann_zeta_tail` (Euler–Maclaurin).
- **Loader (2000) saddle-point log-pmfs:** `stirling_error`, `binom_deviance`,
  `ln_poisson_pmf`, `ln_binomial_pmf`.
- **Bessel:** `bessel_j` (Miller downward recurrence), `bessel_y` (series /
  asymptotic / upward recurrence).

**Used by the platform:** only `gamma` (by `scirust-fractional`). The rest is
latent capability available to future inference / statistics layers (Layer 5).

### 2.2 `scirust-fractional` — fractional operators

Flat public API; `FractionalOrder` is a validated newtype over `f64`,
constructible only for the strict open interval `0 < α < 1`.

- `caputo_l1_uniform(samples, step, order)` — left Caputo derivative at the
  final sample, classical L1 scheme, uniform grid. Exact on piecewise-linear
  data.
- `caputo_l1_nonuniform(samples, sample_times, order)` — L1 Caputo with the
  kernel integrated exactly per sub-interval; correct on strictly increasing
  non-uniform grids.
- `grunwald_letnikov_weights(order, len)`, `riemann_liouville_gl_uniform(...)` —
  Grünwald–Letnikov coefficients (recursive) and the RL derivative.
- `FractionalError` typed enum (invalid order/step, empty/too-few samples,
  non-finite / non-monotonic times, mismatched lengths).

### 2.3 `scirust-sim` — deterministic ODE engine

The integration substrate. Core is dependency-free and a pure function of its
inputs (bit-identical trajectories for identical arguments).

- **`System` trait** (`y' = f(t, y)`, in-place `derivatives`), **`SecondOrderSystem`**
  (`q'' = a(t,q,v)`), **`FirstOrderForm`** adapter.
- **Integrators:** `simulate` (fixed-step classical RK4), `simulate_adaptive`
  (Dormand–Prince 5(4), FSAL, I-controller, Hairer automatic initial step),
  `simulate_second_order` (symplectic/semi-implicit Euler).
- **`Trajectory`** (sampled `(t, y)`), **`SimError`** typed enum
  (`BadInput`, `DimMismatch`, `NonFinite`, `StepUnderflow`), with a hard
  `MAX_STEPS` budget and an adaptive step-underflow floor.
- Optional `stiff` (Rosenbrock-W / Backward Euler via `scirust-stiff`) and `rl`
  (RL-environment bridge) features — **off** in the relativity build.

### 2.4 `scirust-relativity` — geometry core (Layer 1)

The established-GR geometry engine. Trait-based, const-generic over dimension.

- **Traits:** `Metric<D>` (`components`), `Connection<D>` (`christoffel`,
  indexed `[rho][mu][nu]`).
- **Operations:** `invert_metric` (deterministic Gauss–Jordan, partial
  pivoting), `determinant` (deterministic Gauss elimination), `metric_norm`,
  `numerical_christoffel` (central-difference Levi-Civita), and
  `CurvatureTensors<D>` — Riemann, Ricci, Ricci scalar,
  Einstein, and Kretschmann via central differences of the Christoffel symbols,
  with typed errors and a guarantee of no non-finite output.
- **Backgrounds:** `Minkowski`, `MinkowskiSpherical`, `Schwarzschild`,
  `IsotropicSchwarzschild`, `ReissnerNordstrom`, `Kerr`, `DeSitter`,
  `AntiDeSitter`, and `Flrw<S: ScaleFactor>` (spatially flat cosmology, generic
  over the scale factor — `ExponentialScaleFactor` = de Sitter,
  `PowerLawScaleFactor` = radiation/matter eras). Analytic connections
  everywhere except `Kerr` and `IsotropicSchwarzschild`, which use
  `numerical_christoffel` (a disclosed finite-difference connection). Static
  spherically symmetric `f(r)` metrics share the `static_spherical` lapse helper.
- **Parallel transport:** `transport_along_segment` / `transport_along_polyline`
  / `holonomy_defect` for vectors, plus `transport_covector_*` and
  `transport_covariant_tensor_*` for lower-index objects, integrating the
  transport ODE with the `scirust-sim` RK4 engine; validated by flat exactness,
  the holonomy/curvature identity, and metric compatibility (norm preservation,
  metric self-transport, index-lowering commutation, contraction preservation).
- **Geodesic deviation:** `integrate_geodesic_deviation` solves the covariant
  Jacobi system `(x, u, xi, w)` with the Riemann source, validated against the
  coordinate separation of two nearby geodesics (convention-free ground truth).
- **Exponential / logarithm maps:** `geodesic_exponential` (geodesic endpoint at
  unit affine parameter) and `geodesic_logarithm` (Newton-shooting inverse using
  the exponential's finite-difference Jacobian and the crate's matrix inverse),
  validated by flat exactness and the curved round-trip identity.
- **Orthonormal frames (tetrads):** `orthonormal_tetrad` builds a local
  orthonormal frame `{e_a}` (`g(e_a,e_b) = eta_ab`, `e_0` the normalized
  four-velocity) for a timelike observer by metric Gram-Schmidt, exposing the
  shared `OrthonormalTetrad<D>` type; validated by orthonormality, completeness
  (any vector reconstructs from its frame components), agreement with the
  closed-form metric split, and preservation of orthonormality under parallel
  transport of the legs. The experimental worldline observer tetrad delegates to
  this primitive rather than duplicating the Gram-Schmidt construction.
- **Synge's world function and bitensors:** `world_function` /
  `world_function_with_gradients` give the biscalar
  `sigma(x', x) = (1/2) g(x')(v, v)` (with `v = log_{x'}(x)`) and its
  first-derivative gradient bitensors `sigma^{mu'} = -v`, `sigma^mu = -log_x(x')`;
  `van_vleck_determinant` gives the second-derivative bitensor
  `Delta(x', x) = sqrt|g(x')| / (sqrt|g(x)| det J)` from the exponential-map
  Jacobian. All reuse the geodesic logarithm map; validated by flat exactness,
  base/field symmetry, the field-point fundamental identity
  `2 sigma = g(x) sigma^mu sigma^mu`, the gradient round trip, van Vleck
  flat/coincidence unity and `Delta(x', x) = Delta(x, x')` symmetry, and the
  known maximally-symmetric expansion `(Delta - 1)/sigma -> Lambda/3`.
- **Linearized gravity (Layer 2 opening):** `LinearizedField` computes the
  field equations to first order in a metric perturbation `h = g - eta` about
  Minkowski — the linearized Riemann/Ricci/Einstein tensors and the
  trace-reversed perturbation — by central differences of a perturbation sampler.
  Validated by weak-field-Schwarzschild vacuum, the Newtonian Poisson limit, gauge
  invariance of `R^(1)`, and an `O(h^2)` cross-check against the nonlinear
  `CurvatureTensors`. This is the first Layer 2 (Covariant Gravity Workbench)
  capability; it lives in this crate for now (a dedicated `scirust-gravity` crate
  is an option if the surface grows). See `docs/LAYER_2_COVARIANT_GRAVITY.md`.
- **PPN extraction (Layer 2):** the `ppn` module extracts the Eddington–Robertson
  `gamma`, `beta` from static isotropic weak-field metrics (`extract_ppn`), by
  deterministic polynomial extrapolation of the finite-radius effective estimators
  to zero compactness (a self-contained least-squares solver). It enforces an
  explicit isotropic-coordinate contract (`StaticIsotropicMetric` +
  `IsotropicChartAdapter`, which rejects areal charts via a conformal-flatness
  check), carries its own typed `PpnError`, and reports window/order/resolution/
  conditioning diagnostics as an *estimated* uncertainty (not a bound). Validated
  against exact and contaminated synthetic metrics and exact isotropic
  Schwarzschild (`gamma = beta = 1`). See `docs/LAYER_2_PPN.md`.
- **Einstein–Hilbert action variation (Layer 2):** the `action` module
  numerically varies `S = integral (R - 2 Lambda) sqrt(-g) d^4x`
  (`einstein_hilbert_action_variation`) for a static, axisymmetric background
  against a compact test perturbation, by a central difference in the amplitude
  of a Simpson-quadratured action, and compares it to the analytic-integrand
  prediction `-integral sqrt(-g) E^{ab} h_{ab}` from the Einstein tensor. The
  integrand uses a new metric-only nested-difference Ricci scalar
  (`ricci_scalar_from_metric`, a Layer 1 generalization reusing the curvature
  assembly); the static + axisymmetric symmetry reduces the 4D variation to a 2D
  `(r, theta)` integral and the compact bump kills the boundary term. Carries a
  typed `ActionError`; validated by vacuum stationarity (Schwarzschild,
  `Lambda`-matched de Sitter — residual ~`O(dx^4)`), a mismatched-`Lambda`
  nonzero cross-check, and grid convergence. A numerical approximation, never an
  exact variation. See `docs/LAYER_2_ACTION_VARIATION.md`.
- **3+1 (ADM) kinematics (Layer 2):** the `adm` module decomposes a 4-metric on
  the constant-time foliation into lapse, shift, spatial metric, and extrinsic
  curvature (`adm_decomposition`), and evaluates the Gauss–Codazzi Hamiltonian
  and momentum constraints (`adm_constraints`), which vanish for exact solutions.
  It reuses `ricci_scalar_from_metric` and `numerical_christoffel` at `D = 3` for
  the spatial curvature and connection, carries a typed `AdmError`, and is
  validated on Schwarzschild, static de Sitter, FLRW, and the horizon-penetrating
  `PainleveGullstrand` background (the non-zero-shift, spatially-varying-`K`
  oracle) plus the algebraic reconstruction identity. The bridge to Layer 3; it
  evolves nothing. See `docs/LAYER_2_ADM.md`.
- **Backgrounds:** Minkowski (Cartesian and spherical), Schwarzschild,
  isotropic Schwarzschild, Reissner–Nordström, Kerr, de Sitter, anti-de Sitter,
  spatially flat FLRW, and Painlevé–Gullstrand (a horizon-penetrating
  Schwarzschild foliation, `Metric` only — the ADM non-zero-shift oracle).
- **Dynamics:** `GeodesicSystem<C, D>` implements `scirust_sim::System` (state
  `[x, u]`, RHS the geodesic equation `−Γ^ρ_{μν} u^μ u^ν`).
- **Errors:** `RelativityError` (non-finite coordinate/metric/curvature/transport/
  world-function, singular metric, invalid difference/affine step, non-convergent
  logarithm map, and tetrad failures: invalid floor, non-timelike frame vector,
  non-finite leg, degenerate frame).
- **Tests:** 130 across eighteen integration-test files (curvature, geometry,
  kerr, reissner_nordstrom, schwarzschild, coordinate_independence,
  parallel_transport, covariant_transport, flrw, geodesic_deviation,
  exponential_map, tetrad, synge, van_vleck, linearized, ppn, action, adm).
- **Benchmarks:** `benches/geometry_core.rs`, `benches/ppn.rs`,
  `benches/action.rs`, and `benches/adm.rs` (`criterion`, `harness = false`) time
  the hot paths — Christoffel, `invert_metric`, the curvature engine, RK4
  transport, world-function / van Vleck shooting, PPN sampling / extrapolation /
  extraction, the metric-only Ricci scalar / action variation, and the ADM
  decomposition / constraints. Wall-clock, so machine-dependent — the library
  stays deterministic.

### 2.5 `scirust-nonlocal-relativity` — hereditary worldline dynamics (Layer 4, experimental)

An **experimental, phenomenological** layer: it integrates a single test
particle's worldline on a *fixed* GR background with a projected fractional
(Caputo) velocity-memory force added to the RHS. It does **not** modify the
Einstein equations, the Einstein tensor, or the stress-energy tensor, and claims
no empirical validity (scope stated verbatim in its README, `lib.rs`, and the
v2 docs).

Its strength is a clean **composable-policy architecture** — four orthogonal
trait axes plus a modulator axis, all consumed by one generic physics core:

- `HistoryBackend<D>` — `CompleteUniformHistory` (exact oracle),
  `BoundedShortMemoryHistory` (windowed approximation).
- `HistoryTransport<D>` — `IdentityHistoryTransport`,
  `DiscreteConnectionTransport` (Heun discrete parallel transport).
- `MemoryLaw<D>` — `CaputoCoordinateMemory`, `ModulatedCaputoCoordinateMemory`,
  and non-uniform variants.
- `WorldlineStepper<D>` — `SemiImplicitEulerStepper`, `HeunPeceStepper`.
- `HistoryModulator<D>` — `Identity`, `SchwarzschildKretschmann`,
  `ReissnerNordstromField` (phenomenological reweighting hooks).

Adaptive control is centralized in `adaptive_control` (one scaled-RMS error
norm, one accept/reject/step-control routine, shared by both the embedded
Heun–Euler and the step-doubling controllers). Exact validation oracles exist
for flat cylindrical-Minkowski transport (`charts`) and Schwarzschild
circular-orbit transport (`curved_transport`, a fixed-term matrix exponential).
Its observer-frame diagnostic (`build_orthonormal_tetrad`, `tetrad_state_error`)
delegates to the geometry core's `orthonormal_tetrad` and re-exports the shared
`OrthonormalTetrad`, so the Gram-Schmidt construction is not duplicated across
crates; the delegation maps each geometry-core error back to the crate's own
`NonlocalRelativityError` variant, keeping the public API and its tetrad tests
bit-for-bit unchanged. Determinism is enforced by `.to_bits()` bit-identity tests. ~45-variant typed
`NonlocalRelativityError`. 13 integration-test files, 12 examples.

### 2.6 `experiments/nonlocal-relativity-v2`

Twenty deterministic experiment binaries, each printing a `#`-prefixed
metadata header (units, determinism, provenance commit, scientific-category
label) then CSV, with finiteness validation and a non-overclaiming
interpretation. They split by scientific category: the **experimental,
phenomenological** worldline set (`adaptive_convergence`, `history_retention`,
`complexity_scaling`, `bounded_memory_error`, `kerr_fd_sensitivity`,
`modulation_sensitivity`) and the **established-GR** set — geometry core
(`curvature_invariants`, `coordinate_independence`, `parallel_transport`,
`covariant_transport`, `flrw_curvature`, `geodesic_deviation`, `exponential_map`,
`orthonormal_tetrad`, `world_function`, `van_vleck_determinant`) plus the Layer 2
`linearized_gravity`, `ppn_extraction`, `action_variation`, and `adm_kinematics`.

## 3. Validated mathematics (oracle inventory)

The platform's credibility rests on validating numerics against **exact
closed-form results**. Current oracle coverage:

| Domain | Oracle | Where |
|--------|--------|-------|
| Metric inversion | analytic inverse of diagonal / known metrics | geometry tests |
| Christoffel symbols | analytic symbols vs `numerical_christoffel` | schwarzschild / kerr tests |
| Curvature | Minkowski `R_{abcd}=0` exactly; Schwarzschild Ricci-flat, `K=48M²/r⁶`; dS/AdS `R_{μν}=Λg`, `R=4Λ`, `G=−Λg`, `K=8Λ²/3`; Riemann symmetries + first Bianchi | curvature tests |
| Coordinate independence | `R`, `K` agree across Cartesian/spherical Minkowski and areal/isotropic Schwarzschild charts | coordinate_independence tests |
| Caputo L1 | derivative of a constant = 0 (bit-exact); of a linear/power `t → t^{1−α}/Γ(2−α)` (exact); GL of `t²` = `Γ(3)/Γ(3−α)` | fractional tests |
| Flat transport | closed-form cylindrical-Minkowski parallel transport | nonlocal `exact_transport` |
| Curved transport | `exp(−λA)` circular-orbit transport | nonlocal `curved_transport` |
| Special functions | `scipy.special` / `mpmath` reference values | scirust-special inline tests |
| ODE integration | RK4 vs `e^{−t}` (order-4 convergence); DP5(4) tolerance behavior; symplectic energy bound | scirust-sim inline tests |

## 4. Numerical methods and assumptions

- **Curvature engine:** second-order central differences of the Christoffel
  symbols. For analytic-connection backgrounds this is a single FD layer
  (relative error ~`1e-6` near the optimal step, then a roundoff floor). For
  `Kerr` / `IsotropicSchwarzschild` the connection is itself an FD, so curvature
  is a *nested* FD with larger, disclosed error (~`1e-5`–`1e-6`).
- **Fractional operators:** `O(N)` per evaluation at the final point ⇒ `O(N²)`
  over a full trajectory; no history compression or short-memory truncation in
  `scirust-fractional` itself (the worldline crate adds a bounded-memory backend
  separately). Order restricted to `α ∈ (0,1)`.
- **Integration:** geodesics use fixed-step RK4 (the quadratic four-velocity
  norm is only approximately conserved — RK4 preserves only *linear* invariants
  exactly). The worldline crate adds its own first-order adaptive controllers.
- **Determinism:** fixed reduction/summation order everywhere; the one matrix
  exponential uses a fixed Taylor term count. Bit-for-bit reproducibility is a
  tested property, not an aspiration.

## 5. Architectural strengths

- **Trait-based, const-generic, composable.** Backgrounds are just
  `Metric + Connection`; the worldline layer's five orthogonal trait axes make
  new physics additive. Adding a background, memory law, transport, stepper, or
  modulator requires no change to existing engines.
- **One shared adaptive infrastructure** (`adaptive_control`) — the earlier
  duplication between the two controllers has already been consolidated into one
  error norm and one control routine.
- **Typed errors, no panics, no `unsafe`, deterministic** across the whole
  subgraph, enforced by crate attributes, CI, and bit-identity tests.
- **Disciplined scientific labeling.** Established GR (`scirust-relativity`),
  numerical methods, and phenomenological models (`scirust-nonlocal-relativity`)
  are separated in code, docs, and per-experiment headers.
- **Reproducibility as a first-class artifact:** a pinned-nightly CI gate and a
  deterministic experiment suite with provenance headers.

## 6. Technical debt and duplication

Concrete, actionable items (candidate follow-up increments — none blocks current
work):

- ~~**Duplicated fractional-memory helpers** in `scirust-nonlocal-relativity`:
  `nonuniform_caputo_velocity_memory` byte-for-byte identical in `adaptive.rs`
  and `nonuniform_memory.rs`; the transported-modulated variant likewise.~~
  **Resolved:** both callers now delegate to one shared `nonuniform_kernel`
  module (`pub(crate)` builders), moved verbatim so the crate's bit-identity
  golden tests still pass. This removed the highest-value duplication.
- **Two near-identical step-evaluation functions**: `evaluate_step_with_policy`
  (`lib.rs`, used by fixed-step + step-doubling) and `evaluate_adaptive_step`
  (`adaptive.rs`); they share the metric/Christoffel/gr/force/diagnostics body
  and differ mainly in memory-law dispatch. Consolidation is possible but
  higher-risk (touches both control paths).
- **Repeated Christoffel finiteness scans** against different error variants
  (`transport.rs`, `curved_transport.rs`, `lib.rs`) — a small shared validator
  would remove the churn.
- **Adaptive API asymmetry**: `AdaptiveSimulationPolicy` exposes a modulator but
  no memory law/stepper; `AdaptiveStepperPolicy` the opposite. Deliberate and
  documented, but a learning cost.
- **`BoundedShortMemoryHistory` uses `Vec::remove(0)`** (`O(N)` shift) instead
  of a ring buffer — minor inefficiency on the approximate path.
- **Documentation gaps:** no `README.md` for `scirust-relativity` or
  `scirust-fractional`; no standalone API / validation / numerical-methods
  handbooks (content is dispersed across the paper, the STATUS doc, and inline
  rustdoc). `scirust-fractional` cites its algorithms by name only (no
  author/year), unlike the well-cited `scirust-special`.
- **CI coverage gap:** the `nonlocal-relativity-experiments` crate is
  *path-triggered* and text-scanned for forbidden markers, but never compiled,
  tested, or run in CI (it is absent from every `-p` list, and the examples job
  only covers `scirust-nonlocal-relativity`). An experiment could break silently.
- **Documentation drift (largely resolved):** 20 experiment binaries on disk;
  the experiments README now itemises both the six phenomenological binaries and
  the fourteen established-GR ones (geometry core plus the Layer 2
  `linearized_gravity`, `ppn_extraction`, `action_variation`, and
  `adm_kinematics`). The remaining drift is the paper's reproduction section,
  which still lists only 3.
- ~~**No performance benchmarks anywhere** in the subgraph — no `benches/`, no
  `criterion`/`iai`/`divan`.~~ **Resolved:** `criterion` wall-clock benches now
  cover the geometry-core hot paths (`scirust-relativity/benches/geometry_core.rs`)
  and the `O(N^2)` Caputo history
  (`scirust-nonlocal-relativity/benches/caputo_history.rs`). Benchmark timings are
  machine-dependent and not bit-reproducible (inherent to timing); the
  deterministic, reproducible operation-count proxy in `complexity_scaling`
  remains the companion measure.

## 7. Extension points (designed-in, reuse-first)

- **New spacetime:** implement `Metric<D>` (+ `Connection<D>`, analytic or via
  `numerical_christoffel`). Proven background-agnostic (Kerr, dS/AdS, isotropic
  Schwarzschild all added this way).
- **New curvature diagnostic:** build on `CurvatureTensors` (as the
  coordinate-independence work does).
- **New memory physics / modulation:** implement `MemoryLaw<D>` and/or
  `HistoryModulator<D>`.
- **New history / transport strategy:** implement `HistoryBackend<D>` /
  `HistoryTransport<D>`.
- **New integrator:** implement `WorldlineStepper<D>` (fixed-step and
  step-doubling adaptive pick it up automatically) or, at the `scirust-sim`
  level, add a free function over `System` (RK45-Fehlberg, Verlet, Gauss–Legendre).
- **New fractional operator:** add to the flat `scirust-fractional` surface
  (higher orders `α>1`, Riesz, fast/compressed history are explicitly deferred
  contracts).

## 8. Missing research capabilities (mapped to the six-layer vision)

Relative to [`PLATFORM_ROADMAP.md`](PLATFORM_ROADMAP.md):

- **Layer 1 (Geometry Core) — partial.** Present: metrics, connections,
  curvature, geodesics, nine backgrounds (including spatially flat FLRW),
  coordinate-independence diagnostics, a reusable parallel-transport engine with
  first-class holonomy, covector/tensor transport, geodesic-deviation (Jacobi)
  fields, exponential / logarithm maps, local orthonormal frames (tetrads) in
  the geometry core (the worldline observer tetrad now delegates to this shared
  primitive), and Synge's world function with its gradient bitensors and the
  van Vleck–Morette determinant, with the first `criterion` performance benches
  over the geometry-core hot paths and the `O(N^2)` Caputo history. The
  differential-geometry surface and its near-term benchmark goal are now in place;
  Layers 2–6 are the next frontier.
- **Layer 2 (Covariant Gravity Workbench) — opening.** The design note
  (`docs/LAYER_2_COVARIANT_GRAVITY.md`) fixes the scope, category labels, and
  oracles. Four increments are delivered and validated: **linearized gravity**
  (`LinearizedField`: the weak-field Einstein equations to first order in
  `h = g - eta`), **PPN extraction** (the `ppn` module: Eddington–Robertson
  `gamma`, `beta` from static isotropic weak-field metrics, `docs/LAYER_2_PPN.md`),
  the **Einstein–Hilbert action variation** (the `action` module: a numerical
  `delta S / delta g` reproducing `G_{mu nu} + Lambda g_{mu nu} = 0` for static
  axisymmetric vacua, `docs/LAYER_2_ACTION_VARIATION.md`), and **3+1 (ADM)
  kinematics** (the `adm` module: lapse/shift/spatial-metric/extrinsic-curvature
  and the Gauss–Codazzi constraints, `docs/LAYER_2_ADM.md`). This completes the
  near-term Layer 2 sequence and bridges to Layer 3; full symbolic action
  machinery is deferred.
- **Layer 3 (Numerical Relativity) — absent.** No perturbation theory,
  self-force, or ADM/BSSN evolution. `scirust-sim` lacks the dense output,
  event detection, and constraint-preserving/projection integration such work
  needs.
- **Layer 4 (Gravitational Memory) — experimental prototype only.** The
  fractional-memory worldline layer exists and is clearly labelled
  phenomenological; standard/Christodoulou memory and detector response are
  absent.
- **Layer 5 (Astrophysical Inference) — absent.** No waveforms, noise models,
  likelihood, or samplers. (`scirust-special` provides much of the statistical
  substrate.)
- **Layer 6 (Relativistic Navigation) — absent.** No proper-time/Shapiro/redshift
  observables assembled into a navigation engine.

## 9. Recommended near-term sequence

Additive, each validated against an oracle, each one PR:

1. **Coordinate-independence diagnostics** (this increment) — invariants agree
   across charts; adds `MinkowskiSpherical` and `IsotropicSchwarzschild`.
2. **Consolidate the duplicated non-uniform Caputo memory helpers** — *done*
   (`nonuniform_kernel` module; bit-identity-preserving).
3. **Reusable parallel-transport engine** in the geometry core, with holonomy
   validated against the curvature tensor — *done*. Tetrads and Jacobi-field
   (geodesic-deviation) checks against maximally symmetric closed forms are the
   natural follow-on, layered on this same segment integrator.
4. **FLRW background** with its exact curvature oracle — *done* (generic over a
   scale factor; de Sitter and radiation/matter eras, coordinate-independence
   cross-check against the static de Sitter chart).
5. Geodesic-deviation (Jacobi) fields, exponential/logarithm maps,
   covector/tensor transport, and tetrads in the geometry core (generalising the
   worldline observer tetrad, not duplicating it) — *done*.
6. **Synge's world function and its gradient bitensors** — *done* (built on the
   geodesic logarithm map; validated by flat exactness, symmetry, the
   fundamental identity, and the gradient round trip).
7. **The van Vleck–Morette determinant** — *done* (from the exponential-map
   Jacobian on a new deterministic `determinant`; validated by flat/coincidence
   unity, `Delta(x', x) = Delta(x, x')` symmetry, and the known
   maximally-symmetric `(Delta - 1)/sigma -> Lambda/3` expansion).
8. **First performance benchmarks** — *done* (`criterion` wall-clock benches for
   the geometry-core hot paths and the `O(N²)` Caputo history; timings are
   machine-dependent, the deterministic op-count proxy remains the reproducible
   companion). This closes the near-term Layer 1 sequence.
9. **Layer 2 (Covariant Gravity Workbench)** — *near-term sequence complete*.
   Design note (`docs/LAYER_2_COVARIANT_GRAVITY.md`); **linearized gravity**
   (`LinearizedField`), **PPN extraction** (the `ppn` module,
   `docs/LAYER_2_PPN.md`), the **Einstein–Hilbert action variation** (the
   `action` module, `docs/LAYER_2_ACTION_VARIATION.md`), and **3+1 (ADM)
   kinematics** (the `adm` module, `docs/LAYER_2_ADM.md`) are all *done*. Next:
   **Layer 3 (Numerical Relativity)**, which evolves the ADM data in time.

Layers 2–6 open only after Layer 1 is broad and solid, each with a design note
fixing its oracles and category labels before code lands.

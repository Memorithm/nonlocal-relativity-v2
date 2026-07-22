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
  pivoting), `metric_norm`, `numerical_christoffel` (central-difference
  Levi-Civita), and `CurvatureTensors<D>` — Riemann, Ricci, Ricci scalar,
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
  / `holonomy_defect`, integrating the transport ODE with the `scirust-sim` RK4
  engine; validated by flat exactness, metric compatibility, and the
  holonomy/curvature identity against the Riemann tensor.
- **Geodesic deviation:** `integrate_geodesic_deviation` solves the covariant
  Jacobi system `(x, u, xi, w)` with the Riemann source, validated against the
  coordinate separation of two nearby geodesics (convention-free ground truth).
- **Exponential / logarithm maps:** `geodesic_exponential` (geodesic endpoint at
  unit affine parameter) and `geodesic_logarithm` (Newton-shooting inverse using
  the exponential's finite-difference Jacobian and the crate's matrix inverse),
  validated by flat exactness and the curved round-trip identity.
- **Dynamics:** `GeodesicSystem<C, D>` implements `scirust_sim::System` (state
  `[x, u]`, RHS the geodesic equation `−Γ^ρ_{μν} u^μ u^ν`).
- **Errors:** `RelativityError` (non-finite coordinate/metric/curvature,
  singular metric, invalid difference step).
- **Tests:** 45 across five files (curvature, geometry, kerr,
  reissner_nordstrom, schwarzschild) plus the coordinate-independence suite.

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
Determinism is enforced by `.to_bits()` bit-identity tests. ~45-variant typed
`NonlocalRelativityError`. 13 integration-test files, 12 examples.

### 2.6 `experiments/nonlocal-relativity-v2`

Seven deterministic experiment binaries, each printing a `#`-prefixed metadata
header (units, determinism, provenance commit, scientific-category label) then
CSV, with finiteness validation and a non-overclaiming interpretation:
`adaptive_convergence`, `history_retention`, `complexity_scaling`,
`bounded_memory_error`, `kerr_fd_sensitivity`, `modulation_sensitivity`,
`curvature_invariants` (established-GR category), plus the new
`coordinate_independence`.

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
- **Documentation drift:** 7 experiment binaries on disk; the experiments README
  documents 6 (omits `curvature_invariants`); the paper's reproduction section
  lists 3.
- **No performance benchmarks anywhere** in the subgraph — no `benches/`, no
  `criterion`/`iai`/`divan`. The only "performance" measure is a deterministic
  operation-count proxy in `complexity_scaling`/`history_retention`. The
  roadmap's performance-benchmark goals are entirely unmet.

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
  first-class holonomy, geodesic-deviation (Jacobi) fields, and exponential /
  logarithm maps. Missing: tetrads / orthonormal frames in the geometry core
  (one exists, observer-specialised, in the worldline crate — to be generalised,
  not duplicated), covector/tensor transport, and bitensors / Synge world
  function.
- **Layer 2 (Covariant Gravity Workbench) — absent.** No symbolic action,
  variational calculus, automatic field-equation derivation, or PPN/weak-field
  machinery.
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
5. Geodesic-deviation (Jacobi) fields and exponential/logarithm maps — *done*.
   Tetrads in the geometry core (generalising the worldline observer tetrad) and
   covector/tensor transport remain.
6. **First performance benchmarks** (curvature engine, Caputo `O(N²)` history)
   to close the empty-benchmarks gap.

Layers 2–6 open only after Layer 1 is broad and solid, each with a design note
fixing its oracles and category labels before code lands.

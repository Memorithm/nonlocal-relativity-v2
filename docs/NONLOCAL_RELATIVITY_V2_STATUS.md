# Nonlocal Relativity V2 — Status

This document is a single-page status snapshot of the `scirust-nonlocal-relativity`
experimental layer after Phases 1–4. It is a summary; the authoritative
technical description remains
[`docs/EXPERIMENTAL_NONLOCAL_RELATIVITY.md`](EXPERIMENTAL_NONLOCAL_RELATIVITY.md)
and [`scirust-nonlocal-relativity/README.md`](../scirust-nonlocal-relativity/README.md).

**This is an experimental research layer for fractional-memory test-particle
worldline dynamics on a fixed general-relativistic background. It does not
modify the Einstein field equations, the Einstein tensor, or the
stress-energy tensor; it computes no feedback of matter on curvature; it is
not a covariant field theory; and no empirical validation is claimed.**

## Delivered Features

### Phase 1/2 (prior work, unchanged)

- Fixed-background test-particle worldline integration on any
  `Metric<D> + Connection<D>` background.
- Ordinary state equation `du^rho/dlambda = a_GR^rho + F_memory^rho`, with a
  projected Caputo L1 velocity-memory force on the right-hand side.
- `HistoryBackend` (`CompleteUniformHistory` exact oracle,
  `BoundedShortMemoryHistory` explicit approximation), `HistoryTransport`
  (`IdentityHistoryTransport`), `MemoryLaw` (`CaputoCoordinateMemory`), and
  `WorldlineStepper` (`SemiImplicitEulerStepper`, `HeunPeceStepper`)
  separated as independent, composable components.
- Self-convergence studies and Schwarzschild exterior chart diagnostics.

### Phase 3 (this work)

- `HistoryEntry<D>`: typed accepted sample (coordinates, velocity,
  parameter), and additive `HistoryBackend::push_entry`/`entry` and
  `HistoryTransport::transport_segment` trait methods with compatible
  defaults.
- `DiscreteConnectionTransport`: deterministic Heun predictor-corrector
  discretization of `dV^mu/dlambda = -Gamma^mu_(alpha beta) u^alpha V^beta`,
  applied once per accepted segment to every retained history vector.
- `ParameterizationMode` (`AffineParameter` / `NormalizedTimelikeProperTime`)
  via `simulate_nonlocal_worldline_with_mode`, and `affine_trajectory_proper_time`
  diagnostics.
- `CylindricalMinkowski`, exact Jacobian coordinate/velocity transforms, and
  `examples/coordinate_covariance.rs`.

### Phase 4 (this work)

- `HistoryModulator<D>` (`IdentityHistoryModulator`,
  `SchwarzschildKretschmannModulator`) and `ModulatedCaputoCoordinateMemory<M>`,
  a `MemoryLaw` that applies a deterministic dimensionless scalar weight to
  each retained sample before the Caputo evaluation.
- `examples/curvature_modulated_memory.rs`.

### Follow-up: exact flat-spacetime transport oracle (this work)

- `exact_cylindrical_minkowski_transport`: a closed-form (non-discretized)
  parallel transport for flat spacetime in the Cartesian/cylindrical chart
  pair, exploiting the path-independence of transport under a curvature-free
  connection.
- `transport_vector_along_polyline`: exposes `DiscreteConnectionTransport`'s
  per-segment mechanism directly over an explicit waypoint list, independent
  of the full simulation/backend machinery.
- `examples/exact_transport_convergence.rs`: demonstrates
  `DiscreteConnectionTransport`'s numerical error converging to the exact
  oracle under path refinement (second-order, `~3.5e-5` to `~3.5e-8` from 4
  to 128 waypoints with the shipped parameters) — a direct validation against
  a known-exact answer, not just against another discretization.

### Follow-up: non-uniform Caputo operator and proper-time-based memory (this work)

- `scirust_fractional::caputo_l1_nonuniform`: the L1 Caputo scheme
  generalized to an explicitly non-uniform temporal grid (additive; the
  existing `caputo_l1_uniform` is unchanged). Validated independently:
  exactness for linear functions on a non-uniform grid, and numerical
  agreement with `caputo_l1_uniform` on a uniform grid.
- `proper_time_caputo_velocity_memory`: a pure post-hoc diagnostic that
  evaluates Caputo velocity memory against an already-computed trajectory's
  own non-uniform proper-time axis (from `affine_trajectory_proper_time`),
  resolving the "must never be passed to `caputo_l1_uniform`" gap that
  function's own documentation flagged. Does not touch the live integration
  loop.
- `examples/proper_time_memory_comparison.rs`: compares affine- and
  proper-time-based memory on a Schwarzschild trajectory across refinement
  levels.

### Follow-up: curved-background exact transport, adaptive stepping, and a second modulator (this work)

- `exact_schwarzschild_circular_orbit_transport`,
  `schwarzschild_circular_orbit_four_velocity`,
  `schwarzschild_circular_orbit_angular_velocity`: an exact, closed-form
  parallel transport along a circular equatorial geodesic orbit in a
  **curved** (Schwarzschild) background, via the matrix exponential of the
  constant transport generator that this special path family produces (a
  different mechanism than flat spacetime's path-independence). Validated
  against `DiscreteConnectionTransport` under refinement (second order,
  matching the flat-spacetime oracle's pattern) and against two exact,
  refinement-independent conservation laws of parallel transport under any
  metric-compatible connection (`g(V,V)` and `g(V,u)` constant).
- `AdaptiveNonlocalConfig`, `simulate_nonlocal_worldline_adaptive`: an
  integration loop that chooses its own non-uniform affine-parameter step,
  using the classical embedded Heun-Euler pair for error control (no extra
  acceleration evaluation beyond one ordinary Heun step) and
  `caputo_l1_nonuniform` for the memory force against the resulting
  non-uniform history. Coordinate memory only (no `HistoryTransport` or
  `HistoryModulator` composability yet). Cross-validated against a very fine
  fixed-step `HeunPeceStepper` run on the same trajectory (agreement to
  `1.0e-4` in the shipped example/test parameters); demonstrated reaching a
  comparable or better accuracy than an 800-step fixed run with as few as 23
  accepted steps at a loose tolerance in `examples/adaptive_worldline.rs`.
- `ReissnerNordstrom` (in `scirust-relativity`): a second fixed background
  (static, charged, spherically symmetric), with exact analytic Christoffel
  symbols (the same general `f(r)`-metric formula Schwarzschild already
  uses, with `f(r) = 1 - 2M/r + Q^2/r^2`), validated against
  `numerical_christoffel` and against exact agreement with `Schwarzschild`
  at `charge = 0`.
- `ReissnerNordstromFieldModulator`: a curvature-adjacent modulator built
  from the electromagnetic field invariant `F_(mu nu) F^(mu nu) = 2 Q^2 /
  r^4` of the Reissner-Nordström background's radial Coulomb field — not a
  curvature invariant, and not the Kretschmann scalar
  `SchwarzschildKretschmannModulator` uses. Same `beta = 0` bit-identical
  bypass pattern.
- `examples/schwarzschild_orbit_transport.rs`,
  `examples/adaptive_worldline.rs`,
  `examples/reissner_nordstrom_field_modulation.rs`.

### Follow-up: adaptive composability and a Kerr background (this work)

- `simulate_nonlocal_worldline_adaptive_with_policy`, `AdaptiveSimulationPolicy<H, T, M>`:
  composes the adaptive integrator with `HistoryTransport` and
  `HistoryModulator` — `DiscreteConnectionTransport`,
  `SchwarzschildKretschmannModulator`, and `ReissnerNordstromFieldModulator`
  now all work under adaptive stepping, individually and together, by
  reusing `HistoryBackend::push_entry` (the same segment-transport
  mechanism the fixed-step integrators already use). `simulate_nonlocal_worldline_adaptive`
  is now defined as the `IdentityHistoryTransport` +
  `IdentityHistoryModulator` + `CompleteUniformHistory` special case,
  verified to reproduce its pre-composition numbers bit-for-bit against a
  captured golden regression value.
- `Kerr` (in `scirust-relativity`): a rotating (stationary, axisymmetric)
  background in standard Boyer-Lindquist coordinates. Unlike every other
  background in this crate, its connection is evaluated by central finite
  differences (`numerical_christoffel`) rather than an exact analytic
  formula — an explicit, disclosed tradeoff given the far greater algebraic
  complexity of Kerr's Christoffel symbols. Validated at `spin = 0` against
  `Schwarzschild` (metric bit-for-bit, Christoffel symbols to
  finite-difference tolerance), against the known Lense-Thirring
  frame-dragging sign, and end-to-end via a stationary-observer worldline
  simulation showing spin-scaling frame dragging emerge from the ordinary
  geodesic equation with no Kerr-specific code beyond the metric and
  connection.
- `examples/adaptive_transported_modulated.rs`, `examples/kerr_worldline.rs`.

## Validations Performed

- `cargo fmt --all -- --check` clean.
- `cargo test --locked -p scirust-fractional -p scirust-relativity -p scirust-nonlocal-relativity`:
  all tests and doctests passing (exact counts in the phase commit messages
  and PR description).
- `cargo clippy --locked -p scirust-fractional -p scirust-relativity -p scirust-nonlocal-relativity --all-targets -- -D warnings`
  clean.
- All crate examples (`schwarzschild_memory`, `convergence_study`,
  `coordinate_covariance`, `curvature_modulated_memory`,
  `exact_transport_convergence`, `proper_time_memory_comparison`,
  `schwarzschild_orbit_transport`, `adaptive_worldline`,
  `reissner_nordstrom_field_modulation`, `adaptive_transported_modulated`,
  `kerr_worldline`) run to completion and produce deterministic CSV output.
- Bit-for-bit regression: every Phase 1/2 test file is unmodified and passes
  unchanged; Phase 3/4 additions include explicit bit-identity tests for the
  compatibility paths (`beta = 0`, identity transport, affine mode).
- A dedicated test mechanically scans this crate's own source for item
  declarations (`struct`/`enum`/`trait`/`fn`) whose name suggests a modified
  field equation, Einstein tensor, or stress-energy structure, and fails if
  one is found.
- `DiscreteConnectionTransport` validated directly against the exact
  flat-spacetime oracle under refinement, in addition to the cross-chart
  disagreement comparison: 7 dedicated tests plus a CSV-producing example,
  confirming second-order convergence to a known-exact answer.
- `DiscreteConnectionTransport` validated the same way against the
  circular-equatorial-orbit exact oracle in the **curved** Schwarzschild
  background (11 dedicated tests plus a CSV-producing example), including
  two exact conservation-law checks independent of the discrete method.
- `ReissnerNordstrom`'s exact analytic Christoffel symbols cross-checked
  against `numerical_christoffel` (central finite differences) and against
  bit-for-bit/machine-precision agreement with `Schwarzschild` at zero
  charge (10 dedicated tests).
- `simulate_nonlocal_worldline_adaptive` cross-validated against a very fine
  independent fixed-step `HeunPeceStepper` run on the same trajectory, in
  addition to bit-for-bit determinism, non-uniform-step, target-reaching,
  and typed-error-on-budget-exhaustion tests (16 dedicated tests including
  the transport/modulation composition follow-up below).
- `simulate_nonlocal_worldline_adaptive_with_policy` validated against a
  bit-for-bit golden regression captured from the pre-composition
  implementation, and tested composing `DiscreteConnectionTransport`,
  `SchwarzschildKretschmannModulator`, `ReissnerNordstromFieldModulator`,
  `BoundedShortMemoryHistory`, and combinations thereof, each checked for
  finiteness and (where applicable) a measurable, non-identical departure
  from the identity baseline.
- `Kerr` validated at `spin = 0` against `Schwarzschild`'s metric
  (bit-for-bit) and exact analytic Christoffel symbols (finite-difference
  tolerance), plus a frame-dragging sign check and a symmetric-metric check
  (11 dedicated tests).

## Complexities (as actually implemented)

| Component | Cost |
|---|---|
| Raw coordinate memory (`CaputoCoordinateMemory` + `IdentityHistoryTransport`) | `O(D * N^2)` over `N` fixed steps (unchanged from Phase 1/2) |
| Bounded short memory (`BoundedShortMemoryHistory`) | `O(D * N * W)` for window `W` (unchanged) |
| Discrete parallel transport (`DiscreteConnectionTransport`) | `O(D^3 * N^2)`: `O(N)` transported vectors per accepted step (`O(N^2)` total), each a Christoffel contraction (`O(D^3)`) |
| Curvature modulation (`ModulatedCaputoCoordinateMemory`) | Adds `O(1)` work per retained sample per evaluation on top of whichever transport/backend it wraps |
| Non-uniform Caputo (`caputo_l1_nonuniform`) | `O(N)` per evaluation, same order as `caputo_l1_uniform` |
| Proper-time memory (`proper_time_caputo_velocity_memory`) | `O(D * N)` for one post-hoc evaluation over an `N`-sample trajectory (builds the proper-time axis once, then one non-uniform Caputo evaluation per component) |
| Exact circular-orbit transport (`exact_schwarzschild_circular_orbit_transport`) | `O(1)`: one Christoffel evaluation plus a fixed-size (4x4) matrix exponential, independent of path length or refinement |
| Adaptive worldline integration (`simulate_nonlocal_worldline_adaptive`) | `O(D * N^2)` in the best case (no rejections) over `N` accepted steps, matching the coordinate-memory baseline's order; each accepted-step attempt costs one extra `O(D * N)` non-uniform Caputo evaluation at the trial point, and a rejected attempt is discarded and retried at a smaller step |
| Reissner-Nordström field modulation (`ReissnerNordstromFieldModulator`) | Adds `O(1)` work per retained sample per evaluation, identical in structure to `SchwarzschildKretschmannModulator` |
| Adaptive worldline with transport/modulation (`simulate_nonlocal_worldline_adaptive_with_policy`) | Same order as plain adaptive integration; `DiscreteConnectionTransport` composed in adds the same `O(D^3)` per-retained-vector transport cost the fixed-step path pays, applied once per accepted (and once per trial) segment |
| Kerr connection (`Kerr::christoffel`) | `O(1)` per evaluation: `numerical_christoffel` evaluates the metric `2D+1 = 9` times (central differences in each of 4 directions) and inverts one 4x4 metric, versus `O(1)` for Schwarzschild/Reissner-Nordström's single closed-form evaluation — asymptotically the same order, a larger constant factor |

## Assumptions

- The supplied background implements `Metric<D>` and `Connection<D>`
  consistently; the crate validates finiteness, not geometric consistency.
- The worldline is non-null; `|g(u,u)|` must exceed a configurable positive
  floor (affine mode) or stay within a configurable tolerance of `-1`
  (proper-time mode).
- `NormalizedTimelikeProperTime` assumes a `(-,+,+,+)` signature background;
  it detects an incompatible signature indirectly, through the same
  closeness-to-`-1` check used for drift.
- `CylindricalMinkowski` and `SchwarzschildKretschmannModulator` assume the
  crate's established 4D coordinate ordering (`(t, r, phi/theta, z/phi)` as
  appropriate to each background).
- `SchwarzschildKretschmannModulator` assumes evaluation points are in the
  Schwarzschild exterior (`r` strictly greater than the horizon radius).
- `exact_cylindrical_minkowski_transport` assumes flat spacetime and a path
  staying within a simply connected region of the chart (not winding around
  `r = 0`); path-independence of transport does not hold once curvature is
  non-zero.
- `proper_time_caputo_velocity_memory` assumes every sampled state is
  timelike; its proper-time axis is only as accurate as
  `affine_trajectory_proper_time`'s first-order quadrature.
- `exact_schwarzschild_circular_orbit_transport` assumes the transport path
  is exactly a circular equatorial geodesic orbit (`r` constant, `theta =
  pi/2`, the specific four-velocity `schwarzschild_circular_orbit_four_velocity`
  returns) at a radius strictly exceeding `3 M`; it is not valid for any
  other path, including non-circular, non-equatorial, or non-geodesic ones.
- `ReissnerNordstrom` assumes sub-extremal parameters (`charge^2 <
  mass^2`), guaranteeing two distinct, real horizons; extremal and
  super-extremal parameters are rejected at construction.
- `ReissnerNordstromFieldModulator` assumes evaluation points are in the
  Reissner-Nordström exterior (`r` strictly greater than the outer horizon
  radius), mirroring `SchwarzschildKretschmannModulator`'s Schwarzschild
  exterior assumption.
- `simulate_nonlocal_worldline_adaptive_with_policy` assumes its supplied
  `HistoryBackend` retains complete or explicitly bounded history
  consistently with the fixed-step architecture's own assumptions about
  that backend; adaptivity changes only the affine-parameter grid, not the
  history-retention contract.
- `Kerr` assumes sub-extremal parameters (`spin^2 < mass^2`); its
  `Christoffel` symbols are a finite-difference approximation with a fixed
  internal difference step, not an exact analytic result, unlike every
  other background in this crate.

## Limitations

- The memory kernel remains coordinate-dependent; transport and modulation
  change this quantitatively (measurably smaller cross-chart disagreement,
  shrinking under refinement) but do not eliminate it.
- `DiscreteConnectionTransport` is a discrete, segment-by-segment
  approximation, not an exact bitensor propagator, and its cost is
  asymptotically worse than the coordinate-memory baseline. It now has a
  known-exact reference for validation in the flat-spacetime case only;
  curved backgrounds (`Schwarzschild`) still have no exact reference in this
  crate.
- `NormalizedTimelikeProperTime` validates but does not adapt the step; drift
  beyond tolerance is a hard error, not a corrected trajectory.
- `SchwarzschildKretschmannModulator`'s `beta` and reference length are free,
  uncalibrated phenomenological parameters specific to the Schwarzschild
  exterior chart.
- Event handling, error estimation beyond the embedded Heun-Euler pair, and
  history compression are not included anywhere in the crate.
- `proper_time_caputo_velocity_memory` is a post-hoc diagnostic, not an
  adaptive integrator: it resamples an already uniformly-stepped trajectory
  onto an estimated proper-time axis after the fact. A genuinely adaptive
  live loop exists (`simulate_nonlocal_worldline_adaptive`), and it now
  composes with geometric transport and curvature/field modulation via
  `simulate_nonlocal_worldline_adaptive_with_policy`.
- `exact_schwarzschild_circular_orbit_transport` is exact only for the
  circular-equatorial-orbit special case; it has no known-exact reference
  for an eccentric, inclined, or otherwise general curved path, and neither
  it nor `DiscreteConnectionTransport` closes that gap.
- `simulate_nonlocal_worldline_adaptive`'s step-doubling-free embedded error
  estimate (Heun-Euler) is a standard but relatively simple adaptive scheme;
  it does not include event handling, dense output, or higher-order
  embedded pairs. It still does not compose with `WorldlineStepper` or
  `MemoryLaw` (as opposed to `HistoryTransport`/`HistoryModulator`, which it
  now does), since those two traits thread a single fixed `NonlocalConfig`
  step through their signatures.
- `ReissnerNordstromFieldModulator`'s `beta` and reference length are, like
  `SchwarzschildKretschmannModulator`'s, free, uncalibrated phenomenological
  parameters, specific to the Reissner-Nordström exterior chart.
- `Kerr`'s Christoffel symbols are evaluated by finite differences, not an
  exact analytic formula; they carry a small, difference-step-dependent
  truncation error that every other background in this crate avoids. No
  Kerr-specific circular-orbit, transport, or modulation construction is
  provided; the crate's Kerr example uses only a simple stationary initial
  state for exactly this reason.

## Future Work (not implemented here)

- General curved-path exact parallel transport (a full bitensor propagator),
  as a replacement for the discrete segment-by-segment transport in the
  general case. (The flat-spacetime case has an exact closed-form reference,
  `exact_cylindrical_minkowski_transport`, and the circular-equatorial-orbit
  case in a curved background now also has one,
  `exact_schwarzschild_circular_orbit_transport`; neither extends to a
  general curved path, where neither flatness's path-independence nor a
  circular orbit's constant-transport-generator argument applies.)
- Composing the adaptive-step integrator with `WorldlineStepper` or
  `MemoryLaw` themselves (as opposed to `HistoryTransport`/
  `HistoryModulator`, composed this round via
  `simulate_nonlocal_worldline_adaptive_with_policy`). Both traits still
  thread a single fixed `NonlocalConfig` step through their signatures
  (`StepperContext` for the former), which a variable step size cannot
  satisfy without changing those contracts.
- An exact analytic Kerr connection, as a replacement for the
  finite-difference Christoffel symbols delivered this round; Kerr-specific
  transport, modulation, or circular-orbit constructions (this round
  deliberately used only a simple stationary initial state, avoiding any
  Kerr orbital-mechanics formula).
- Curvature or field modulators for backgrounds other than Schwarzschild,
  Reissner-Nordström, or Kerr (all delivered), or built from invariants
  other than the Kretschmann scalar or the electromagnetic field invariant
  delivered so far — for example the Ricci-squared invariant that a
  non-vacuum background other than the traceless-stress-tensor
  Reissner-Nordström case would make nonzero.
- **Any investigation of modified field equations is not future work for
  this crate — it is permanently out of scope.** This is not a deferred
  item awaiting a future round; it is excluded by the crate's own
  non-negotiable scientific boundary (see the next section), and no future
  change to this crate should attempt it.

## Explicitly Forbidden Claims

This crate, in code, comments, documentation, and any future change, must
never claim or imply that it:

- modifies the Einstein field equations;
- implements fractional Einstein equations;
- modifies the Einstein tensor;
- modifies the stress-energy tensor;
- computes feedback of matter on curvature;
- constitutes a complete covariant theory;
- establishes new physics;
- carries experimental validation.

`DiscreteConnectionTransport` is a discrete numerical approximation, not an
exact bitensor propagator or a proof of covariance.
`SchwarzschildKretschmannModulator` is a phenomenological scalar reweighting,
not a consequence of general relativity, a quantum-gravity prediction, or an
experimentally derived law.
`exact_cylindrical_minkowski_transport` is exact only for flat spacetime; it
is not a general bitensor propagator and must never be described as valid
for curved backgrounds.
`exact_schwarzschild_circular_orbit_transport` is exact only for a circular
equatorial geodesic orbit in Schwarzschild; it is not a general bitensor
propagator either, and must never be described as valid for a general
curved path (eccentric, inclined, non-geodesic, or otherwise).
`simulate_nonlocal_worldline_adaptive` is a standard embedded-Runge-Kutta
adaptive step-size scheme applied to this crate's existing ordinary state
equation; it is not a new numerical method, does not change the state
equation, and does not extend to a claim about physical accuracy beyond
what the underlying discretization already provides.
`ReissnerNordstromFieldModulator` is a phenomenological scalar reweighting,
like `SchwarzschildKretschmannModulator`; it is not a consequence of general
relativity or electromagnetism, a quantum-field-theory prediction, or an
experimentally derived law.
`Kerr`'s Christoffel symbols are a finite-difference numerical
approximation, not an exact analytic result like `Schwarzschild`'s or
`ReissnerNordstrom`'s; this must always be disclosed alongside any use of
`Kerr`, and no claim of machine-precision exactness may be made for it.
`simulate_nonlocal_worldline_adaptive_with_policy` composes standard,
independently-established components (an embedded Runge-Kutta pair, the
existing `HistoryTransport`/`HistoryModulator` contracts); the composition
itself is not a new physical claim beyond what each component already
carries.

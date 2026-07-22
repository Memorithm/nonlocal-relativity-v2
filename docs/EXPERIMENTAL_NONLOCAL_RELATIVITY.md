# Experimental Nonlocal Relativity Layer

This document describes an **EXPERIMENTAL** SciRust research layer for
fractional-memory test-particle worldline dynamics on a fixed
general-relativistic background.

It is not a theory of fractional Einstein equations. It does not modify the
Einstein field equations, the Einstein tensor, the stress-energy tensor,
matter-generated curvature, or established general relativity. No empirical
validation is claimed.

## Scientific Boundary

A single-paragraph statement of exactly what this layer is and is not, so no
reader has to reconstruct it from the sections below.

- **What equation is solved.** A single ordinary first-order state equation
  for a test particle's coordinates and contravariant velocity,
  `du^rho/dlambda = a_GR^rho + F_memory^rho`, where `a_GR` is the ordinary
  geodesic acceleration and `F_memory^rho = -kappa P^rho_sigma m^sigma` is a
  projected Caputo velocity-memory force. It is **not** a Caputo fractional
  differential equation for the state itself, and **not** any field equation.
- **What is held fixed.** The background metric and connection are supplied
  externally and never change; nothing here solves the Einstein or Maxwell
  equations or computes any backreaction of the particle on curvature.
- **What is phenomenological.** The memory coupling `kappa`, the fractional
  order `alpha`, and every modulator (`SchwarzschildKretschmannModulator`,
  `ReissnerNordstromFieldModulator`) with its free `beta` and reference
  length. None is calibrated against data or derived from a field theory.
- **What is coordinate dependent.** The Caputo memory is evaluated
  componentwise in whatever chart the background supplies, so it is chart
  dependent by construction. The Phase 1 scaled error norm reduces the
  adaptive controllers' sensitivity to the chart and to component magnitudes
  but is itself componentwise and **not** a covariant measure; discrete
  transport reduces but does not remove the memory's chart dependence.
- **What has an exact oracle.** Two, and only two, transport families: flat
  spacetime (`exact_cylindrical_minkowski_transport`) and Schwarzschild
  circular equatorial geodesic orbits
  (`exact_schwarzschild_circular_orbit_transport`). No exact reference is
  currently implemented for a general curved path.
- **What has only self-convergence / fine-grid reference.** Everything else:
  `run_convergence_study` (self-convergence at `h`, `h/2`, `h/4`) and the
  adaptive/retention comparisons against a fine independent fixed-step run.
  These are numerical consistency diagnostics, not validations of the model.
- **What has not been empirically validated.** The physical model, in full.
  No result here is evidence of new physics, and none should be described as
  experimentally confirmed.

## Index Conventions

- Coordinates are `x^rho`.
- Contravariant coordinate velocity is `u^rho = dx^rho / d lambda`.
- Greek indices such as `rho`, `mu`, `nu`, and `sigma` range over the supplied
  coordinate dimension `D`.
- Repeated indices are summed in the equations.
- The metric is covariant, `g_(mu nu)`.
- Christoffel symbols are indexed as `Gamma^rho_(mu nu)`.
- All quantities use the chart and geometric units of the supplied background.

## Equations

The fixed background implements both `Metric<D>` and `Connection<D>`. The
ordinary geodesic acceleration is

```text
a_GR^rho = - Gamma^rho_(mu nu)(x) u^mu u^nu .
```

The complete uniformly sampled velocity history is used to define a coordinate
Caputo memory vector:

```text
m^rho(lambda_n) = CaputoDerivative_alpha[u^rho](lambda_n) .
```

The implementation delegates this operation to
`scirust_fractional::caputo_l1_uniform`. It does not duplicate the fractional
operator.

The current velocity is lowered with the metric:

```text
u_sigma = g_(sigma nu) u^nu .
```

The metric norm is

```text
s = g_(mu nu) u^mu u^nu .
```

For non-null worldlines, the memory vector is projected orthogonally to the
current velocity:

```text
P^rho_sigma = delta^rho_sigma - u^rho u_sigma / s .
```

The experimental memory force is

```text
F_memory^rho = - kappa P^rho_sigma m^sigma .
```

The total trajectory-level equation is

```text
du^rho / d lambda = a_GR^rho + F_memory^rho .
```

The derivative of the state remains an ordinary first derivative in affine
parameter. The fractional operator appears in the history-dependent force on
the right-hand side; this implementation must not be described as a Caputo
fractional differential equation for the state variables themselves.

The diagnostic residual

```text
u_rho F_memory^rho
```

is exposed so the projection can be audited numerically.

## Assumptions

- The spacetime geometry is fixed externally.
- The worldline is a test-particle trajectory and does not source curvature.
- The sampled worldline is non-null; `|s|` must exceed a configurable positive
  floor.
- The fractional order is in the first-release supported interval
  `0 < alpha < 1`.
- `kappa` is finite and non-negative.
- Positive `kappa` is a phenomenological damping-like coupling, not a new
  fundamental constant.
- The discretization is coordinate-dependent.
- The uniform affine-parameter step is finite and positive.
- Invalid non-finite numerical values are rejected rather than repaired.

## Numerical Architecture

The implementation separates the numerical responsibilities used by the
experimental worldline layer:

- `HistoryBackend<D>` stores accepted velocity samples and reports how many
  retained samples are used by memory evaluation.
- `HistoryTransport<D>` maps retained samples into the current coordinate
  frame before memory evaluation. The current production implementation is
  coordinate identity/no-transport.
- `MemoryLaw<D>` evaluates the coordinate memory vector from retained,
  transported samples. The current production law is Caputo L1 coordinate
  velocity memory.
- `WorldlineStepper<D>` advances the ordinary first-order state equation once
  the total acceleration has been evaluated.

Transport is abstracted separately from memory because future research may
study transported histories or chart-specific comparison maps without changing
the Caputo L1 stencil or the storage backend. The current identity transport
does not make the model covariant; it preserves the Phase 1 coordinate-memory
contract.

## Transported Memory and Discrete Parallel Transport (Phase 3)

The Phase 1/2 transport contract, `HistoryTransport::transport_velocity`,
receives a retained sample's bare velocity components and the current
worldline state, with no source point for the sample. That is enough for
coordinate-identity transport, but a transport cannot carry a vector between
two distinct tangent spaces using components alone. Phase 3 adds this
additively:

- `HistoryEntry<D>` is a typed accepted sample: coordinates, contravariant
  velocity, and the accepted parameter value.
- `HistoryBackend::push_entry` and `HistoryBackend::entry` are new trait
  methods built on `HistoryEntry`. Their default implementations fall back to
  `push_velocity`/`sample`-only behavior, so a backend that does not override
  them is honestly limited to coordinate-identity transport; it must not
  claim otherwise.
- `HistoryTransport::transport_segment` is a new trait method, called once
  per accepted (or provisional) segment for every currently retained vector,
  with the segment's `from`/`to` states and the fixed background connection.
  Its default is the identity, matching `transport_velocity`'s existing
  contract.
- `DiscreteConnectionTransport` implements `transport_segment` with a single
  Heun predict-evaluate-correct-evaluate step of the linear transport
  equation `dV^mu/dlambda = -Gamma^mu_(alpha beta) u^alpha V^beta` per
  accepted segment:

```text
1. evaluate the transport derivative at the segment start;
2. predict the vector at the segment end;
3. evaluate the connection and velocity at the segment end;
4. correct with the average of the two derivatives.
```

Because every retained vector is advanced by exactly one such step each time
a new segment is accepted, transport accumulates along the actual accepted
worldline polyline — this matters because parallel transport under a curved
connection is path-dependent, not just endpoint-dependent. `IdentityHistoryTransport`
is unchanged; its `transport_segment` uses the trait default (identity), so
the original coordinate-memory model is preserved exactly, including
bit-for-bit output for existing callers.

`DiscreteConnectionTransport` **is not**:

- an exact analytic bitensor propagator;
- an analytic parallel propagator;
- a proof of covariance for this model.

It is a discrete numerical approximation whose error accumulates with the
segment step and the number of transported segments. `examples/coordinate_covariance.rs`
shows, on a controlled experiment comparing Cartesian and cylindrical charts
of flat Minkowski spacetime, that this approximation measurably reduces
disagreement between charts compared to raw coordinate memory, and that the
disagreement shrinks under refinement — without claiming exact chart
agreement.

**Complexity.** Transporting every retained vector once per accepted segment
costs `O(N)` transport evaluations per step, `O(N^2)` over `N` accepted
steps; each evaluation is a Christoffel contraction, `O(D^3)`. The
discrete-transport pipeline therefore costs `O(D^3 * N^2)` overall, more
expensive than the `O(D * N^2)` raw coordinate-memory baseline. This is the
complexity actually implemented, not the baseline complexity restated.

## Parameterization Modes (Phase 3)

`ParameterizationMode` does not change the fixed-step numerical scheme —
both modes advance the same uniform `parameter_n = n * h` — it only fixes
what the parameter means and what is validated:

- `AffineParameter` is the default: no timelike assumption, bit-for-bit
  compatible with every pre-Phase-3 entry point.
- `NormalizedTimelikeProperTime { tolerance }` interprets the step as a
  proper-time step. The initial state must be timelike with metric norm
  within `tolerance` of `-1` under a `(-,+,+,+)`-signature background; every
  subsequently sampled step's metric norm must stay within `tolerance` of
  `-1`. Null and spacelike states, and states from an incompatible-signature
  background, are all rejected by the same check: their norm is not close to
  `-1`. No automatic four-velocity renormalization is ever performed; a
  violation is a typed `ProperTimeNormDrift` error.

`simulate_nonlocal_worldline_with_mode` wraps `simulate_nonlocal_worldline_with_policy`
with these checks before and after the run; it does not alter the stepper.

`affine_trajectory_proper_time` is a separate, purely diagnostic estimate of
proper time elapsed along an *affine*-parameter trajectory of timelike
states. Its quadrature is the left-endpoint rule
`Delta tau_n ~= h * sqrt(-g(u,u)_n)`, evaluated at the accepted state at the
start of each step — first-order accurate in `h`. It is not a resampling
onto a uniform proper-time grid: the returned increments are generally
non-uniform, and they must never be passed to
`scirust_fractional::caputo_l1_uniform`, which requires uniform sample
spacing.

## Proper-Time-Based Caputo Memory (Non-Uniform Grid, Follow-Up)

`scirust_fractional::caputo_l1_nonuniform` generalizes the L1 scheme to an
explicitly non-uniform temporal grid: on each subinterval `[t_k, t_(k+1)]`
the derivative is approximated by the finite difference
`(f_(k+1) - f_k) / (t_(k+1) - t_k)`, and the Caputo weight kernel is
integrated exactly over that subinterval using the *actual* sample times
rather than an assumed uniform step. It is validated independently
(exactness for linear functions on a non-uniform grid, and numerical
agreement with `caputo_l1_uniform` on a uniform one) and does not change
`caputo_l1_uniform` at all.

`proper_time_caputo_velocity_memory` uses this operator to evaluate the
Caputo velocity-memory vector of an already-computed affine-parameter
trajectory with respect to its own estimated proper-time axis (built from
`affine_trajectory_proper_time`'s cumulative estimates), resolving the gap
that function's own documentation flagged: "must never be passed to
`caputo_l1_uniform`". This is a **pure post-hoc diagnostic** over an
already-computed trajectory — it does not feed back into the live
integration loop, does not change any accepted state, and is unrelated to
`NormalizedTimelikeProperTime` mode (which advances the state equation with
a *uniform* proper-time step by construction, and never needs this
operator).

Because the memory-force law is built to be orthogonal to the four-velocity,
`g(u,u)` stays close to constant along an accepted trajectory (up to a small,
refinement-shrinking numerical drift), so proper time advances at an
approximately constant rate `c = sqrt(-g(u,u))` relative to the affine
parameter. A Caputo derivative computed against a linearly rescaled
parameter differs from the original by a factor of `c^(-alpha)` — a fact
about the Caputo operator's own scaling behavior, not a discretization
artifact — so the difference between proper-time-based and affine-based
memory is expected to stay roughly constant under refinement, not shrink to
zero. `examples/proper_time_memory_comparison.rs` shows exactly this.

## Coordinate-Chart Comparison (Phase 3)

`CylindricalMinkowski` is the same flat spacetime as `Minkowski`, expressed
in cylindrical coordinates `(t, r, phi, z)`, with connection

```text
Gamma^r_(phi phi)   = -r
Gamma^phi_(r phi)   = Gamma^phi_(phi r) = 1/r
```

which is not identically zero, unlike the Cartesian chart. Exact Jacobian
coordinate and velocity transforms (`cartesian_to_cylindrical_coordinates`,
`cartesian_to_cylindrical_velocity`, and their inverses) let the same
physical initial condition be expressed in either chart and the resulting
trajectories compared.

This comparison has real limits:

- the chart is only regular for `r > 0`; the example and any use of it must
  keep the trajectory away from the coordinate singularity at `r = 0`;
- because the memory model is inherently coordinate-dependent (componentwise
  Caputo differentiation), no finite refinement is expected to drive the
  disagreement between charts to exactly zero, even with transported memory;
  transport reduces the disagreement and improves its convergence trend, it
  does not eliminate the model's coordinate-dependence;
- the comparison is a controlled numerical demonstration on one worldline
  family, not a general covariance proof.

## Exact Flat-Spacetime Transport (Validation Oracle)

`DiscreteConnectionTransport`'s accumulated error can be validated against a
known-exact answer in one specific case: flat spacetime, in the
Cartesian/cylindrical chart pair already established above. Parallel
transport on a flat (curvature-free) manifold is path-independent —
transporting along two different paths between the same two points differs
only by transport around the closed loop formed by one path followed by the
reverse of the other, and a flat connection has trivial holonomy around every
contractible loop. `exact_cylindrical_minkowski_transport` uses this
directly, with no discretization:

1. convert the vector to the Cartesian chart at the source point;
2. leave it unchanged (the Cartesian connection vanishes identically, so
   transport along *any* path does nothing to the Cartesian components);
3. convert back to the cylindrical chart at the destination point.

This is exact **only** for this flat-spacetime chart pair, along paths that
stay within a simply connected region of the chart (in particular, paths
that do not wind around `r = 0`). It is not a general bitensor propagator; it
does not extend to `Schwarzschild` or any other curved background, where
`DiscreteConnectionTransport` remains a discrete approximation with no
closed-form exact counterpart in this crate.

`transport_vector_along_polyline` exposes the same per-segment mechanism
`HistoryBackend::push_entry` uses internally, applied directly to an explicit
waypoint list, independent of the full simulation/backend machinery.
`examples/exact_transport_convergence.rs` transports a fixed test vector
along a straight-line Cartesian path (converted to cylindrical waypoints) at
increasing waypoint density and compares the numerical result to the exact
oracle. With the shipped parameters, the error shrinks by a factor of
essentially `4` each time the waypoint count doubles — second-order
convergence, consistent with the Heun predictor-corrector scheme's local
truncation order — from `~3.5e-5` at 4 waypoints to `~3.5e-8` at 128. This is
a direct validation of `DiscreteConnectionTransport`'s numerical correctness
against a known-exact answer, strictly stronger than comparing two
discretizations to each other.

## Exact Curved-Background Transport (Circular Orbit, Follow-Up)

`exact_cylindrical_minkowski_transport` is exact because flat spacetime has
zero curvature, so transport is path-independent. Schwarzschild has nonzero
curvature everywhere outside its horizon, so that argument does not apply
there — but a *different* structural fact makes one more exact case
available: along a **circular equatorial geodesic orbit** (constant `r`,
`theta = pi/2`, the four-velocity `u = (u^t, 0, 0, u^phi)` returned by
`schwarzschild_circular_orbit_four_velocity`), Schwarzschild's Christoffel
symbols are constant, because they depend only on `r` and `theta`
(stationarity and axisymmetry), both fixed along this path. The parallel
transport equation

```text
dV^mu/dlambda = -Gamma^mu_(alpha beta) u^alpha V^beta
```

therefore reduces, along this one path family, to a **linear,
constant-coefficient** ODE `dV/dlambda = -A V` for the fixed generator
`A^mu_beta = Gamma^mu_(alpha beta) u^alpha`, with the exact closed-form
solution `V(lambda) = exp(-lambda A) V(0)`.
`exact_schwarzschild_circular_orbit_transport` evaluates this directly,
using the same already-validated `Schwarzschild::christoffel` this crate
uses everywhere else (no new Christoffel derivation), and a deterministic
4x4 matrix exponential (scaling-and-squaring with a fixed-length Taylor
series — a standard numerical linear algebra primitive, not a new numerical
method or physics construction).

`schwarzschild_circular_orbit_angular_velocity` returns
`sqrt(M / r^3)` (the general-relativistic form of Kepler's third law, exact
for circular equatorial orbits in these coordinates) and
`schwarzschild_circular_orbit_four_velocity` returns the corresponding
proper-time-normalized four-velocity (`g(u,u) = -1`); both require `r`
finite and strictly greater than `3 M`, the existence bound for a circular
equatorial timelike geodesic (a separate, larger bound at `6 M`, the
innermost *stable* circular orbit, is not enforced, since stability is
irrelevant to evaluating transport along a mathematically valid orbit).

This is validated three ways, none of which assumes the answer:

1. **Metric-compatibility conservation.** For *any* metric-compatible
   connection, parallel transport preserves inner products along any curve:
   `g(V, V)` and `g(V, u)` must stay exactly constant. These are checked
   directly and are true regardless of whether the specific closed-form
   solution above is correct — they test the general transport contract.
2. **Convergence to `DiscreteConnectionTransport`.** Exactly mirroring
   `examples/exact_transport_convergence.rs`'s flat-spacetime pattern,
   `examples/schwarzschild_orbit_transport.rs` shows the discrete scheme's
   numerical error against this new oracle shrinking at second order under
   path refinement (`error_ratio_to_previous` converging to `4.0`).
3. **Round-trip and determinism.** Forward transport followed by the
   reverse (`-delta_lambda`) recovers the original vector; repeated
   evaluation is bit-for-bit identical; `delta_lambda = 0` returns the input
   vector bit-for-bit unchanged.

**This is exact only for a circular equatorial geodesic orbit.** It is
**not** a general bitensor propagator, does **not** extend to eccentric,
inclined, or non-geodesic paths, and must never be described as valid for a
general curved trajectory. For a general curved path, `DiscreteConnectionTransport`
remains the only transport strategy, and — to state the exact-reference status
precisely — **no exact reference is currently implemented for a general curved
path; an exact special-case oracle exists only for the two families
implemented here: flat-spacetime transport (`exact_cylindrical_minkowski_transport`)
and Schwarzschild circular equatorial geodesic orbits
(`exact_schwarzschild_circular_orbit_transport`).**

## Curvature-Modulated Memory (Phase 4)

Phase 4 adds one more additive, opt-in `MemoryLaw`: a deterministic scalar
modulation of retained history vectors, applied componentwise before the
Caputo evaluation. `HistoryModulator<D>` transforms one finite
`HistoryEntry<D>` into a finite, dimensionless scalar weight:

- `IdentityHistoryModulator` always returns `1.0`.
- `SchwarzschildKretschmannModulator` is an explicitly experimental,
  phenomenological instance:

```text
K = 48 M^2 / r^6
q = 1 + beta * L^4 * K
```

`K` is the Schwarzschild Kretschmann scalar; `L` is a strictly positive
reference length that makes the modulation dimensionless; `beta` is a
finite, non-negative phenomenological coefficient, not a new fundamental
constant. Construction requires a strictly positive finite mass and
reference length and a finite non-negative `beta`; evaluation requires a
finite radius strictly outside the horizon and a finite positive resulting
weight.

**`beta == 0.0` is a full bypass.** Evaluation returns exactly `1.0` without
computing the Kretschmann scalar at all, so `ModulatedCaputoCoordinateMemory`
with `beta = 0` reproduces the unmodulated `CaputoCoordinateMemory` pipeline
bit-for-bit whenever the rest of the numerical path (backend, transport,
integrator, mode) is identical.

`ModulatedCaputoCoordinateMemory<M>` implements `MemoryLaw` by applying
`M::weight` to each retained sample's (possibly already-transported)
velocity, componentwise, immediately before the Caputo L1 stencil runs. The
resulting quantity is exactly a Caputo derivative of a dimensionless
*modulated* velocity history. It composes with every other Phase 2/3
component — either history backend, either transport, either fixed-step
integrator, either parameterization mode — because it only changes what
number the Caputo stencil consumes at each retained sample; it does not
change how that sample was stored or transported.

**This law must never be described as:**

- a unique consequence of general relativity;
- a quantum-gravity prediction;
- an experimentally derived law;
- a modification of the Einstein field equations.

No structure resembling a modified field equation, an Einstein tensor, or a
stress-energy tensor is introduced anywhere in this crate. This is
mechanically checked by a dedicated test that scans the crate's own source
for item declarations using such names.

## Reissner-Nordström Field Modulation (Follow-Up)

`SchwarzschildKretschmannModulator` is tied to one background and one
invariant. Schwarzschild is a vacuum solution, so its Ricci tensor vanishes
identically and the Kretschmann scalar is essentially the only nontrivial
polynomial curvature invariant available — there is no independent "other
invariant" to build a second Schwarzschild modulator from. This follow-up
instead adds a second **background**,
`scirust_relativity::ReissnerNordstrom` (a static, spherically symmetric,
electrically charged black hole), and a modulator built from a genuinely
different kind of invariant.

`ReissnerNordstrom`'s metric shares Schwarzschild's `(t, r, theta, phi)`
structure with lapse `f(r) = 1 - 2 M/r + Q^2/r^2` in place of `1 - 2M/r`. Its
Christoffel symbols use the same general formula that this lapse structure
always produces (verified, term by term, against the existing Schwarzschild
implementation while designing this follow-up), so they are exact and
analytic, not finite-differenced. At `charge = 0` the metric and Christoffel
symbols reduce to `Schwarzschild`'s exactly (bit-identical for the metric,
machine-precision-identical for the Christoffel symbols, which take a
different but mathematically equivalent computational path); this crate's
test suite checks this directly, alongside an independent cross-check
against `scirust_relativity::numerical_christoffel`. Construction requires a
strictly positive finite mass and a finite charge satisfying the
sub-extremal bound `charge^2 < mass^2`, guaranteeing two distinct, real
horizons.

`ReissnerNordstromFieldModulator` weights retained history samples by the
electromagnetic field invariant of the background's radial Coulomb field,
not a curvature invariant:

```text
F^2 = F_(mu nu) F^(mu nu) = 2 Q^2 / r^4
q = 1 + beta * L^2 * |F^2|
```

Note the reference-length power `L^2`, not `L^4`: in the geometric units
this crate uses (charge carries the same dimension as mass), `F^2` has
dimension `length^-2`, unlike the Kretschmann scalar's `length^-4`, so a
different power is needed to make the product dimensionless. This crate
uses the Reissner-Nordström metric exactly as it uses Schwarzschild's:
as fixed, externally specified background data. Nothing here solves the
Einstein or Maxwell equations, or computes the electromagnetic field's
backreaction on the metric — the metric formula is simply taken as known,
the same way Schwarzschild's is.

Construction and evaluation follow the same pattern as
`SchwarzschildKretschmannModulator`: a strictly positive finite reference
length and a finite non-negative `beta`; a finite radius strictly outside
the outer horizon; a finite positive resulting weight; and **`beta ==
0.0` is a full bypass** returning exactly `1.0` without computing the field
invariant, reproducing the unmodulated baseline bit-for-bit.

**This law must never be described as:**

- a unique consequence of general relativity or electromagnetism;
- a quantum-gravity or quantum-electrodynamics prediction;
- an experimentally derived law;
- a modification of the Einstein field equations or Maxwell's equations.

## Adaptive-Step Worldline Integration (Follow-Up)

Every integrator described so far advances a **fixed** affine-parameter
step `h` for a fixed number of steps (`NonlocalConfig::step`/`steps`), and
`proper_time_caputo_velocity_memory` only *resamples* an already
uniformly-stepped trajectory after the fact. `simulate_nonlocal_worldline_adaptive`
closes the remaining gap: the live integration loop chooses its own
non-uniform affine-parameter step, using
`AdaptiveNonlocalConfig` (an error tolerance, step bounds, and a target
affine parameter in place of a fixed step and step count).

The step-size controller is the classical **embedded Heun-Euler pair**, a
standard, well-established adaptive-Runge-Kutta technique: the same Euler
predictor and Heun corrector this crate's `HeunPeceStepper` already computes
serve as a first-order/second-order embedded pair, so the local error
estimate costs no acceleration evaluation beyond what one ordinary Heun step
already needs. The error estimate is the **componentwise scaled
root-mean-square norm** `scaled_local_error_norm` (shared with the
step-doubling controller): for each coordinate and velocity component,
`ratio_i = (high_i - low_i) / (abs_tol_i + rel_tol * max(|low_i|, |high_i|))`,
and the norm is `sqrt(mean_i ratio_i^2)`. The tolerances come from an
`AdaptiveTolerance` (a relative tolerance and separate coordinate/velocity
absolute tolerances); `AdaptiveNonlocalConfig::new` seeds a uniform
`AdaptiveTolerance` from the single scalar `error_tolerance` for
compatibility, and `with_tolerance` sets the three fields independently. A
step is accepted when the norm is `<= 1`, and the next step size is proposed
from the standard `safety * norm^(-1/2)` control law (exponent `1/2` for the
lower method order `p = 1`). This scaling improves robustness to the
coordinate chart and to differing component magnitudes; it is **not** a
geometrically invariant error measure and does not establish coordinate
covariance.

A rejected step shrinks and retries, and the retry budget
`max_rejections_per_step` is now actively enforced by **both** adaptive
controllers with identical semantics: exceeding the retry count is a typed
`AdaptiveRejectionBudgetExhausted`, while a proposed shrink that would cross
`min_step` is the **distinct** typed error `AdaptiveMinimumStepExhausted` —
never a silently-accepted out-of-tolerance result. Integration stops once the
accumulated affine parameter reaches `target_affine_parameter` (the final
step is clamped so it does not overshoot); exceeding `max_accepted_steps`
before reaching the target is a typed error (`AdaptiveStepBudgetExhausted`),
never a silently truncated trajectory.

Because the resulting history is non-uniform by construction, the memory
force is evaluated with `caputo_l1_nonuniform` directly against the
accumulated non-uniform affine-parameter axis — not
`caputo_l1_uniform`, and not a post-hoc resample. The returned
`NonlocalTrajectory` therefore samples a generally non-uniform axis: it must
**never** be passed to `affine_trajectory_proper_time`, whose `step`
argument assumes uniform spacing; read `diagnostics()[i].affine_parameter`
directly instead.

`simulate_nonlocal_worldline_adaptive` does not itself reuse `WorldlineStepper`
or `MemoryLaw` — both thread a single fixed `NonlocalConfig` step through
their signatures (`StepperContext` for the former), which a variable step
size cannot satisfy without changing those contracts. (A later follow-up,
"Composing Adaptive Stepping with `MemoryLaw` and `WorldlineStepper`" below,
closes this for `SemiImplicitEulerStepper` specifically, via a different
step-size-control mechanism. `HeunPeceStepper` is now itself sound under a
varying step — `StepperContext` carries the true accumulated affine parameter,
so its predictor no longer reconstructs the parameter from `step_index` — but
it is deliberately not offered through the step-doubling entry point because
that controller's error estimate is specialised to a first-order method, and
adaptive Heun-PECE already exists as this embedded Heun-Euler controller; see
that follow-up for the full reasoning.)
`examples/adaptive_worldline.rs` cross-validates the adaptive path against a
very fine, independent fixed-step `HeunPeceStepper` run (agreement to
`1.0e-4` or better with the shipped parameters) and demonstrates it reaching
comparable or better accuracy than an 800-step fixed run with as few as 23
accepted steps at a loose tolerance — a concrete illustration of *why*
adaptive stepping is useful here, not merely that it runs.

This is a standard numerical technique applied to this crate's existing
state equation. It does not change the state equation, does not claim
improved physical accuracy beyond what the underlying discretization
already provides, and is not a new numerical method.

### Composing Adaptive Stepping with Transport and Modulation (Follow-Up)

`WorldlineStepper` and `MemoryLaw` are the blockers for reusing the
fixed-step architecture directly, but `HistoryTransport` and
`HistoryModulator` were never coupled to a fixed step in the first place:
`HistoryTransport::transport_segment` takes its segment step as an explicit
argument, and `HistoryModulator::weight` takes only a `HistoryEntry` — so
both compose with a variable step size exactly as they are.
`simulate_nonlocal_worldline_adaptive_with_policy` takes an
`AdaptiveSimulationPolicy<H, T, M>` bundling a `HistoryBackend`, a
`HistoryTransport`, and a `HistoryModulator` (mirroring
`NonlocalSimulationPolicy`'s role for the fixed-step path, narrowed to the
three components adaptive stepping actually varies), and reuses
`HistoryBackend::push_entry` — the identical mechanism
`CompleteUniformHistory` and `BoundedShortMemoryHistory` already use to
transport every retained vector across each newly accepted segment for the
fixed-step integrators — to transport history across each adaptively-sized
accepted segment. Each retained sample is weighted by the modulator before
the non-uniform Caputo evaluation, exactly like
`ModulatedCaputoCoordinateMemory` does for the fixed-step path, just with
`caputo_l1_nonuniform` in place of `caputo_l1_uniform`.

`simulate_nonlocal_worldline_adaptive` is now defined as the special case
`AdaptiveSimulationPolicy::new(CompleteUniformHistory::new(),
IdentityHistoryTransport, IdentityHistoryModulator)`; a dedicated
bit-for-bit regression test (captured from the pre-composition
implementation before it was refactored) confirms this reproduces the
original numbers exactly, not merely approximately.
`DiscreteConnectionTransport`, `SchwarzschildKretschmannModulator`, and
`ReissnerNordstromFieldModulator` all compose with adaptive stepping —
individually and together — exactly as they compose with the fixed-step
integrators; `examples/adaptive_transported_modulated.rs` runs all four
combinations of identity/discrete transport and unmodulated/modulated
memory on the same trajectory and tolerance, showing each component
measurably changes the result while every combination remains finite and
well-behaved.

### Composing Adaptive Stepping with `MemoryLaw` and `WorldlineStepper` (Follow-Up)

`WorldlineStepper` and `MemoryLaw` remain blockers for
`simulate_nonlocal_worldline_adaptive_with_policy`'s embedded-pair
controller (the previous section) for the reason already given: both thread
a single fixed `NonlocalConfig` step through their signatures, and
`CaputoCoordinateMemory` specifically applies that one step to the *entire*
retained history via `caputo_l1_uniform`.
`simulate_nonlocal_worldline_adaptive_with_stepper_policy` closes this gap,
but only partially and only by a different mechanism — classical
step-doubling rather than an embedded pair — not by generalizing the
embedded-pair controller itself.

**The memory-law half of the gap is closed completely.** Two new types,
`NonuniformCaputoCoordinateMemory` and
`NonuniformModulatedCaputoCoordinateMemory<M>`, implement the *existing*
`MemoryLaw` trait unmodified: instead of applying one `NonlocalConfig::step`
value to the whole retained history, they read each retained sample's own
recorded `HistoryEntry::parameter` and evaluate `caputo_l1_nonuniform`
directly. Because `MemoryLaw::memory_vector` already receives
`history: &H where H: HistoryBackend<D>` — which already exposes
`HistoryBackend::entry`, since `ModulatedCaputoCoordinateMemory` already
needs it for `HistoryModulator::weight` — no trait signature changed. These
two types compose with the fixed-step architecture too
(`simulate_nonlocal_worldline_with_policy`), not only the adaptive one:
under uniform spacing they produce numerically close results to
`CaputoCoordinateMemory`/`ModulatedCaputoCoordinateMemory` (the two Caputo
evaluators are algebraically equivalent term-by-term under exactly uniform
spacing, though they reach that value by different floating-point paths, so
whether two runs agree to the bit or only closely is a property of the
specific input, not something either type guarantees).

**The stepper half is closed only for `SemiImplicitEulerStepper`, for a
numerical-analysis reason.** Both steppers are now sound under a varying step
size: `StepperContext` carries `current_parameter`, the true accumulated
affine parameter at the accepted state, and `HeunPeceStepper::advance`
computes its provisional predictor point at `current_parameter + config.step`
rather than reconstructing it as `step_index * config.step` (which was exact
only when every accepted step shared one size). The reason `HeunPeceStepper`
is not offered through this step-doubling entry point is therefore no longer a
parameter-formula bug — it is that this controller's error estimate is
*specialised to a first-order method* (the raw one-step/two-half-step
difference is the Richardson estimate only because the divisor `2^p - 1`
equals `1` for `p = 1`), whereas Heun-PECE is second order. Moreover,
step-doubling is the wrong adaptive scheme for Heun: a second-order method
already has a natural embedded first-order partner (its Euler predictor), and
the embedded Heun-Euler pair the *previous* controller uses **is** adaptive
Heun-PECE — it computes the same Heun corrector `HeunPeceStepper::advance`
computes, and additionally exposes the Euler predictor the error estimate
needs. Adding Heun-PECE here would be a strictly inferior duplicate, not a new
capability.

Since `SemiImplicitEulerStepper` alone has no natural embedded higher-order
partner (unlike the Euler-predictor/Heun-corrector pair the previous
section's controller reuses), error control instead uses classical
**step-doubling**: one full trial step of size `h` is compared against two
steps of size `h/2` (with a memory-law evaluation at the midpoint, against a
throwaway provisional history clone — the same pattern
`HeunPeceStepper::advance`'s predictor push already uses). Semi-implicit
Euler is a first-order method (local truncation error `O(h^2)`, same as the
embedded pair's lower method), so the same `1/(p+1) = 0.5` growth/shrink
exponent applies, and for a first-order method the Richardson error estimate
needs no rescaling: the raw one-step/two-half-step difference is already the
local error estimate, used directly.

`AdaptiveStepperPolicy<H, L, T>` bundles a `HistoryBackend`, a `MemoryLaw`,
and a `HistoryTransport` (mirroring `NonlocalSimulationPolicy`'s role for
the fixed-step path, narrowed further than `AdaptiveSimulationPolicy` above:
there is no stepper type parameter, since `SemiImplicitEulerStepper` is the
only sound choice). `simulate_nonlocal_worldline_adaptive_with_stepper` is
the `NonuniformCaputoCoordinateMemory` + `IdentityHistoryTransport` +
`CompleteUniformHistory` special case.
`examples/adaptive_worldline_stepper.rs` runs several `MemoryLaw`/transport
combinations on the same trajectory and tolerance, plus a sanity-anchor row
against the plain entry point.

Both adaptive entry points now enforce
`AdaptiveNonlocalConfig::max_rejections_per_step` identically, through one
shared control routine (`adaptive_control::control_step`): each keeps an
explicit per-accepted-step rejection counter, reset on every acceptance,
returns `AdaptiveRejectionBudgetExhausted` when the retry count is reached,
and returns the distinct `AdaptiveMinimumStepExhausted` when a proposed shrink
would cross `min_step` first. Earlier revisions of the embedded controller
never consulted the retry budget (it stopped only at `min_step`, and reported
that as a rejection-budget error); that inconsistency is fixed.

## Kerr Background (Follow-Up)

Every background used so far (`Minkowski`, `Schwarzschild`,
`ReissnerNordstrom`) is static: its metric does not depend on time, and none
has an off-diagonal `t`-`phi` term. `scirust_relativity::Kerr` (mass `M`,
spin `a = J/M`) is a **rotating** background — stationary and axisymmetric,
but not static, in standard Boyer-Lindquist coordinates:

```text
Sigma = r^2 + a^2 cos^2(theta)
Delta = r^2 - 2 M r + a^2

g_tt   = -(1 - 2 M r / Sigma)
g_tphi = -2 M a r sin^2(theta) / Sigma
g_rr   = Sigma / Delta
g_thetatheta = Sigma
g_phiphi = (r^2 + a^2 + 2 M a^2 r sin^2(theta) / Sigma) sin^2(theta)
```

At `a = 0`, `Sigma = r^2`, `Delta = r^2 - 2Mr`, and every component reduces
algebraically to `Schwarzschild`'s exactly (`g_tphi` vanishes).

**Unlike `Schwarzschild` and `ReissnerNordstrom`, `Kerr`'s connection is
evaluated by central finite differences** (`numerical_christoffel`), not an
exact analytic formula. Kerr's Christoffel symbols are algebraically far
more complex than either of those backgrounds': the metric depends on both
`r` and `theta` (not `r` alone), and the off-diagonal `t`-`phi` term
couples further components together, so many more symbols are nonzero and
mix all four coordinates. Hand-deriving them by the same process used for
`ReissnerNordstrom` (structurally verifying a general formula against
already-correct code, term by term) is not available here — there is no
comparably simple general formula this crate already implements correctly
that Kerr's Christoffels reduce to — so hand-derivation would carry a real
risk of a transcription error with no independent way to catch it in this
codebase. Using `numerical_christoffel`, itself already validated elsewhere
in this crate against every background with an exact analytic connection,
trades exact analytic Christoffels for a small, documented finite-difference
truncation error. This is an explicit, honestly disclosed engineering
choice, not an oversight, and is stated directly in `Kerr`'s own
documentation.

This choice is validated three ways: (1) at `a = 0`, the metric matches
`Schwarzschild`'s bit-for-bit and the finite-difference Christoffel symbols
match `Schwarzschild`'s exact analytic ones to the finite-difference
tolerance; (2) the zero-angular-momentum-observer (ZAMO) angular velocity
`-g_tphi / g_phiphi` is positive for positive spin outside the horizon,
matching the well-known Lense-Thirring frame-dragging direction (a weak-field
expansion of this same ratio reduces to the standard `2 M a / r^3` result);
(3) `examples/kerr_worldline.rs` runs this crate's ordinary worldline and
memory machinery — completely unmodified, with no Kerr-specific code beyond
`Kerr::components`/`Kerr::christoffel` themselves — on a stationary
observer at increasing spin, and frame dragging emerges as expected: the
`phi` coordinate stays exactly zero at `spin = 0` and picks up a positive,
spin-scaling drift at `spin > 0`, entirely from the geodesic equation and
the finite-difference Christoffel symbols working together.

The stationary-observer initial state deliberately avoids any Kerr-specific
circular-orbit formula. Unlike Schwarzschild's, a Kerr circular equatorial
orbit's four-velocity involves a prograde/retrograde asymmetry and an
ISCO shift that depend on `a` in a more complex way; this crate does not
derive, implement, or claim such a formula anywhere.

## Numerical Algorithms

The complete-history Caputo L1 evaluation is the numerical memory oracle for
both fixed-step integrators. The default compatibility policy uses complete
uniform history, coordinate identity transport, the Caputo coordinate-memory
law, and semi-implicit Euler. For the auditable compatibility implementation:

1. Validate configuration, initial coordinates, and initial velocity.
2. Retain the complete velocity history.
3. At sample zero, use a zero memory vector because the Caputo L1 stencil has
   insufficient history.
4. For each later sample, evaluate each component of the Caputo memory vector
   from the complete uniform velocity history.
5. Evaluate the metric, metric norm, Christoffel symbols, ordinary geodesic
   acceleration, projected memory force, and diagnostics in fixed loop order.
6. Reject non-finite metric components, Christoffel symbols, memory values,
   forces, accelerations, generated states, and diagnostics.
7. Reject `|g_(mu nu) u^mu u^nu|` below the configured floor.
8. Advance with semi-implicit Euler:

```text
u_(n+1) = u_n + h a_n
x_(n+1) = x_n + h u_(n+1)
```

The second integrator is named Heun PECE
(`predict_evaluate_correct_evaluate`). It is a second-order predictor-corrector
for the ordinary state equation with a fractional-history force, not a
fractional Adams method for the state equation:

```text
a_n      = a(x_n, u_n, accepted history)
u*       = u_n + h a_n
x*       = x_n + h u*
a*       = a(x*, u*, provisional history including u*)
u_(n+1)  = u_n + h/2 (a_n + a*)
x_(n+1)  = x_n + h/2 (u_n + u_(n+1))
```

The predicted velocity is inserted only into the provisional history used to
evaluate `a*`. The accepted complete history stores the corrected velocity.
All predicted and corrected coordinates and velocities are checked for
finiteness.

There is no RNG, no hidden global state, no parallel reduction, and no
automatic four-velocity renormalization. Metric-norm drift is measured and
reported instead.

The complete-history backend has `O(N)` memory use and `O(D * N^2)` history
cost over `N` fixed steps because each step recomputes direct Caputo histories
for all `D` velocity components. It is exact with respect to this discrete
complete-history contract and remains the oracle used by default.

The bounded short-memory backend retains only the most recent `W >= 2`
accepted velocity samples. It has `O(W)` memory use and `O(D * N * W)` history
cost over `N` fixed steps. It is an explicit approximation: it reports
`Approximate` history diagnostics, retained sample counts, and used sample
counts; it rejects windows smaller than two samples; and it is never selected
automatically from an exact configuration. Constant retained histories still
produce exactly zero Caputo memory because every retained first difference is
zero.

Semi-implicit Euler is a reference integrator for reproducible experiments, not
a precision integrator. Heun PECE is usually more accurate on smooth problems,
but it does not change the model's scientific status or make the force law
covariant.

## Convergence Methodology

The V2 convergence utility runs the same initial condition and final affine
parameter with steps `h`, `h/2`, and `h/4`. It reports:

- endpoint coordinate L2 differences between successive refinements;
- endpoint velocity L2 differences between successive refinements;
- observed self-convergence ratios when the denominator is non-zero;
- endpoint metric-norm drift;
- endpoint memory-force norm.

The `h/4` result is only an internal refinement reference. Self-convergence can
show that two discretizations are approaching each other for a chosen chart,
step sequence, and parameter set. It is not empirical validation, not a proof
of the continuous model, and not a substitute for an exact solution or an
independent numerical oracle.

For Schwarzschild standard exterior coordinates, the crate also exposes
chart-specific diagnostics:

```text
E   = -u_t
L_z = u_phi
s   = g_(mu nu) u^mu u^nu
```

These helpers are explicitly tied to the fixed Schwarzschild exterior
background and do not define generic invariants for every metric.

## Falsifiable Observables

The current API exposes quantities that can be compared across `kappa = 0`,
small positive `kappa`, and independent implementations:

- coordinate trajectory samples `x^rho(lambda_n)`;
- velocity samples `u^rho(lambda_n)`;
- metric-norm drift from the initial sample;
- coordinate L2 norm of the Caputo memory vector;
- coordinate L2 norm of the projected memory force;
- orthogonality residual `u_rho F_memory^rho`;
- coordinate L2 norm of the ordinary geodesic acceleration;
- deviations from an uncoupled geodesic baseline in a fixed chart;
- cross-chart position and velocity disagreement between a Cartesian and a
  cylindrical computation of the same physical worldline, for raw coordinate
  memory and for `DiscreteConnectionTransport` memory, under refinement;
- estimated proper-time increments along an affine-parameter trajectory, and
  metric-norm drift from `-1` under `NormalizedTimelikeProperTime` mode;
- the Schwarzschild-Kretschmann modulation weight and its effect on final
  radius, for `beta = 0` versus small positive `beta`, under refinement;
- the L2 norm difference between affine-parameter-based and
  proper-time-based Caputo velocity memory on the same trajectory, under
  refinement;
- `DiscreteConnectionTransport`'s numerical error against the exact
  circular-equatorial-orbit transport oracle in Schwarzschild, under path
  refinement;
- the Reissner-Nordström electromagnetic-field-invariant modulation weight
  and its effect on final radius, for `beta = 0` versus small positive
  `beta`;
- the adaptive integrator's accepted-step count and combined
  coordinate-and-velocity local error estimate, as a function of the
  configured error tolerance, and its final-state distance to an
  independent fine fixed-step reference;
- the adaptive integrator's final state and memory norm across identity and
  discrete transport, and unmodulated and curvature-modulated memory,
  composed individually and together, on the same trajectory and tolerance;
- the ZAMO angular velocity `-g_tphi / g_phiphi` in Kerr, as a function of
  spin, and the coordinate `phi` drift of a stationary-observer trajectory
  at increasing spin.

These are numerical observables of this discretized model. Agreement or
disagreement with physical data is not claimed.

## Known Limitations

- The memory kernel is applied componentwise in coordinates and is therefore
  coordinate-dependent.
- The current implementation is a trajectory-level constitutive experiment,
  not a covariant field theory.
- Complete-history direct evaluation has quadratic cost in the number of
  samples.
- Bounded short memory reduces history cost but deliberately changes the
  discrete memory model; it is useful only as an approximation to compare
  against the complete-history oracle.
- The Euler update is low order and intended for auditability, not accuracy;
  Heun PECE reduces time-discretization error on smooth tests but still uses a
  coordinate memory force.
- The background connection and metric are assumed to be supplied consistently
  by the caller; the crate validates finiteness but does not prove geometric
  compatibility.
- Null and nearly null worldlines are outside this first implementation.
- Event handling and history compression are not included anywhere in the
  crate.
- `DiscreteConnectionTransport` is a discrete, segment-by-segment
  approximation of parallel transport, not an exact analytic bitensor
  propagator; its discretization error grows with the segment step and the
  number of transported segments, and it costs `O(D^3 * N^2)`, more than the
  `O(D * N^2)` coordinate-memory baseline.
- Transported memory reduces, but does not eliminate, the chart-dependence of
  the underlying coordinate-memory model: componentwise Caputo
  differentiation is not covariant even when its inputs are transported.
- `NormalizedTimelikeProperTime` mode validates but does not adapt: the step
  is still fixed and uniform, drift beyond tolerance is a hard error, and
  there is no mechanism to shrink the step automatically when drift
  approaches the tolerance.
- `CylindricalMinkowski` is regular only for `r > 0`; it is a second chart of
  flat spacetime for comparison purposes, not a new physical background.
- `SchwarzschildKretschmannModulator` is a deliberately simple, hand-chosen
  phenomenological weight; `beta` and the reference length `L` are free
  parameters with no calibration against data, and the modulator is specific
  to the Schwarzschild exterior chart, not a generic curvature-modulation
  mechanism for arbitrary backgrounds.
- Modulation is applied to whatever velocity sample the transport pipeline
  already produced; it does not change the coordinate-dependence discussed
  above, and it composes with, but does not replace, the transported-memory
  discretization error already documented for `DiscreteConnectionTransport`.
- `exact_cylindrical_minkowski_transport` is exact only for flat spacetime in
  the Cartesian/cylindrical chart pair, and only along paths that stay within
  a simply connected region of the chart (paths winding around `r = 0` are
  out of scope). It provides no exact reference for `Schwarzschild` or any
  curved background.
- `proper_time_caputo_velocity_memory` is a post-hoc diagnostic over an
  already-computed trajectory; it does not make the live integration loop
  adaptive, and its non-uniform proper-time axis is itself only a
  first-order-accurate estimate from `affine_trajectory_proper_time`, not an
  independently resolved proper-time integration.
- `exact_schwarzschild_circular_orbit_transport` is exact only for a
  circular equatorial geodesic orbit at a radius strictly exceeding `3 M`
  in Schwarzschild; it provides no exact reference for any other path
  (eccentric, inclined, non-geodesic) or for any other curved background.
- `simulate_nonlocal_worldline_adaptive`'s embedded Heun-Euler error
  estimate is a relatively simple adaptive scheme: no dense output, no
  event handling, no higher-order embedded pair. It now composes with
  `HistoryTransport` and `HistoryModulator` via
  `simulate_nonlocal_worldline_adaptive_with_policy`, so transported and
  curvature/field-modulated memory are no longer exclusive to the
  fixed-step path.
- `ReissnerNordstrom` and `ReissnerNordstromFieldModulator` are, like
  `Schwarzschild` and `SchwarzschildKretschmannModulator`, a fixed
  background and a phenomenological reweighting specific to that
  background's exterior chart; `beta` and the reference length are free,
  uncalibrated parameters.
- `Kerr`'s connection is evaluated by finite differences
  (`numerical_christoffel`), not an exact analytic formula like every other
  background in this crate; it carries a small, step-size-dependent
  truncation error instead of being exact to machine precision. No
  Kerr-specific circular-orbit, transport, or modulation construction is
  provided; `examples/kerr_worldline.rs` uses only a simple stationary
  initial state for exactly this reason.

## Roadmap

### 1. Current Worldline-Memory Model

The current crate implements a fixed-background, test-particle,
Caputo-memory modification of the worldline equation, with two interchangeable
memory pipelines: raw coordinate memory (Phase 1/2, unchanged) and a discrete
parallel-transported memory (Phase 3), an explicit affine-vs-proper-time
parameterization choice (Phase 3), and an optional deterministic curvature
modulation of retained history vectors (Phase 4). It exposes deterministic
diagnostics and explicit failure modes.

### 2. Future Covariant Kernel Research

Phase 3's `DiscreteConnectionTransport` is an initial, explicitly discrete and
approximate step toward transported history — a segment-by-segment numerical
scheme, not an analytic construction. Two exact closed-form cases now exist
as validation oracles for it, neither a general bitensor propagator:
`exact_cylindrical_minkowski_transport` for flat spacetime (exploiting
path-independence of transport under a curvature-free connection), and
`exact_schwarzschild_circular_orbit_transport` for a circular equatorial
orbit in a **curved** background (exploiting the constancy of the transport
generator along that one special path family — a different mathematical
mechanism than path-independence, which curved spacetime does not have in
general). Future research may still study exact or analytic propagators for
a general curved path (where neither shortcut applies), or other covariant
constructions. These remain research directions only and are not
established physics in this crate; none of `DiscreteConnectionTransport`,
`exact_cylindrical_minkowski_transport`, or
`exact_schwarzschild_circular_orbit_transport` preempts or substitutes for
them, and none is presented as a general covariance proof.

`scirust_fractional::caputo_l1_nonuniform`, `proper_time_caputo_velocity_memory`,
and `simulate_nonlocal_worldline_adaptive` together resolve the "proper-time
history sampled at its own resolution" item: the Caputo evaluator no longer
requires a uniform grid, a genuinely non-uniform proper-time axis can be used
post hoc, and the live integration loop can now choose its own non-uniform
affine-parameter step directly, via a standard embedded Heun-Euler error
estimate — not merely resample a uniformly-stepped trajectory after the
fact. `simulate_nonlocal_worldline_adaptive_with_policy` further composes
that adaptive loop with `HistoryTransport` and `HistoryModulator` (see
"Composing Adaptive Stepping with Transport and Modulation" above); what
remains open there is composing adaptivity with `WorldlineStepper` or
`MemoryLaw` themselves, which still thread a single fixed `NonlocalConfig`
step through their signatures.

`scirust_relativity::Kerr` extends the background catalog to a rotating
spacetime, evaluated by finite-difference Christoffel symbols rather than
an exact analytic formula (see "Kerr Background" above) — a deliberate,
disclosed engineering tradeoff, not a claim of the same precision
`Schwarzschild` and `ReissnerNordstrom` provide. Future research may still
study an exact analytic Kerr connection, Kerr-specific transport or
modulation constructions, or additional backgrounds (other rotating or
non-vacuum spacetimes). These remain research directions only.

### 3. Hypothetical Field-Equation Work

This item is **not future work for this crate — it is permanently out of
scope**, excluded by the non-negotiable scientific boundary stated at the
top of this document and in the top-level `README.md`. Any field-equation
investigation would require a separate mathematical and numerical contract,
independent validation, and clear distinction from established general
relativity; this crate does not implement such work and will not, and no
future change to this crate should attempt it. Neither
`SchwarzschildKretschmannModulator` nor `ReissnerNordstromFieldModulator` is
a step toward this item: both are phenomenological scalar reweightings of a
trajectory-level history force, with no Einstein tensor, stress-energy
tensor, or field equation anywhere in their construction.

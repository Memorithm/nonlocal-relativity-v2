# Experimental Nonlocal Relativity Layer

This document describes an **EXPERIMENTAL** SciRust research layer for
fractional-memory test-particle worldline dynamics on a fixed
general-relativistic background.

It is not a theory of fractional Einstein equations. It does not modify the
Einstein field equations, the Einstein tensor, the stress-energy tensor,
matter-generated curvature, or established general relativity. No empirical
validation is claimed.

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
  radius, for `beta = 0` versus small positive `beta`, under refinement.

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
- No adaptive stepping, event handling, error estimation, or history
  compression is included.
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
scheme, not an analytic construction. Future research may still study kernels
defined with exact bitensors, analytic parallel propagators, proper-time
history sampled at its own adaptive resolution, or other covariant
constructions. These remain research directions only and are not established
physics in this crate; `DiscreteConnectionTransport` does not preempt or
substitute for them, and it is not presented as covariant.

### 3. Hypothetical Field-Equation Work

Any future field-equation investigation would require a separate mathematical
and numerical contract, independent validation, and clear distinction from
established general relativity. This crate does not implement such work, and
this roadmap item is not a claim that fractional field equations are
established physics. Phase 4's `SchwarzschildKretschmannModulator` is not a
step toward this item: it is a phenomenological scalar reweighting of a
trajectory-level history force, with no Einstein tensor, stress-energy
tensor, or field equation anywhere in its construction.

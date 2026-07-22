# A Deterministic Numerical Platform for Fractional-Memory Test-Particle Worldline Dynamics on Fixed Relativistic Backgrounds (v2)

**Status: experimental research software.** This document describes the
`scirust-nonlocal-relativity` layer after the v2 hardening round. It is a
numerical-methods report, not a physics result. Nothing here modifies general
relativity, and no claim in it is empirically validated.

---

## 1. Scope and scientific boundary

This platform integrates a single ordinary first-order state equation for one
test particle on a **fixed**, externally supplied general-relativistic
background, with a projected fractional (Caputo) velocity-memory term added to
the right-hand side. In one paragraph:

- **Equation solved:** `du^rho/dlambda = a_GR^rho + F_memory^rho`, an ordinary
  first-order system for coordinates `x^rho` and contravariant velocity
  `u^rho`. It is *not* a fractional differential equation for the state itself,
  and *not* any field equation.
- **Held fixed:** the metric `g_{mu nu}(x)` and connection
  `Gamma^rho_{mu nu}(x)`. Nothing solves the Einstein or Maxwell equations or
  computes backreaction.
- **Phenomenological:** the coupling `kappa`, the fractional order `alpha`, and
  the modulators; none is calibrated to data.
- **Coordinate dependent:** the memory is evaluated componentwise in the
  supplied chart.
- **Exact oracles:** only flat-spacetime transport and Schwarzschild circular
  equatorial geodesic-orbit transport.
- **Otherwise:** self-convergence and fine-grid references only.
- **Not empirically validated:** the physical model, in full.

The taxonomy of evidence used throughout this paper is fixed:

| Term | Meaning here |
|------|--------------|
| **exact analytic result** | closed-form, correct to machine rounding (e.g. `g(u,u) = -1` for a normalized four-velocity) |
| **numerical oracle** | a closed-form or independently-derived answer for a special case, used to validate a discretization (flat / circular-orbit transport) |
| **fine-grid reference** | a much finer independent run of the *same* discretized model; a consistency target, not truth |
| **self-convergence** | comparison of a scheme against refined copies of itself; a stability diagnostic |
| **regression test** | a pinned expected value guarding against unintended change |
| **physical validation** | comparison against experiment or observation — **absent here** |

## 2. Background and motivation

Fractional (Caputo) derivatives are a standard tool for modelling
history-dependent ("memory") effects in viscoelasticity, anomalous diffusion,
and control. The question this layer explores is narrow and purely
computational: *if* one augments a test particle's worldline equation on a
fixed curved background with a projected Caputo velocity-memory force, what does
it take to integrate the resulting hereditary system **deterministically,
reproducibly, and with honest error control**? The physics content of the model
is deliberately minimal and uncalibrated; the engineering content — history
storage, non-uniform Caputo evaluation, adaptive step control, and validation
against exact special cases — is the actual subject.

## 3. Fixed-background worldline dynamics

The ordinary geodesic acceleration in the supplied chart is

```text
a_GR^rho = - Gamma^rho_{mu nu}(x) u^mu u^nu.
```

The background provides `g_{mu nu}` (`Metric<D>`) and `Gamma^rho_{mu nu}`
(`Connection<D>`). Backgrounds implemented: `Minkowski` (Cartesian and, via
`CylindricalMinkowski`, cylindrical), `Schwarzschild`, `ReissnerNordstrom`
(exact analytic connections), and `Kerr` (finite-difference connection,
`numerical_christoffel`). The platform validates finiteness of the metric,
connection, metric norm, and every generated quantity, and returns a typed
error for any violation; it does **not** verify geometric self-consistency of
the supplied background.

The metric norm `s = g_{mu nu} u^mu u^nu` must satisfy `|s| >= metric_norm_floor`
(a positive floor) so the projector below is well defined; a null or nearly
null worldline is rejected, not silently handled.

## 4. Projected fractional-memory force

The memory vector is a componentwise Caputo derivative of the retained velocity
history,

```text
m^rho(lambda_n) = D^alpha_Caputo[ u^rho ](lambda_n),   0 < alpha < 1,
```

projected orthogonally to the four-velocity and scaled by `kappa`:

```text
F_memory^rho = - kappa P^rho_sigma m^sigma,
P^rho_sigma  = delta^rho_sigma - u^rho u_sigma / s,   u_sigma = g_{sigma nu} u^nu.
```

Projection keeps the force in the local rest space of the particle; the
diagnostic residual `u_rho F_memory^rho` (reported per step) is zero up to
rounding. `kappa >= 0` is a finite phenomenological coupling, not a fundamental
constant. The state equation remains ordinary in `lambda`; the fractional
operator appears only as a history-dependent force.

## 5. Uniform and non-uniform Caputo L1 discretizations

The memory uses the classical **L1 discretization** of the Caputo derivative.
On a uniform grid with spacing `h`, for a sampled sequence `f_0, …, f_N`,

```text
D^alpha f(t_N) ≈ h^{-alpha} / Gamma(2 - alpha) *
                sum_{k=0}^{N-1} b_k ( f_{N-k} - f_{N-k-1} ),
b_k = (k+1)^{1-alpha} - k^{1-alpha}.
```

`scirust_fractional::caputo_l1_uniform` implements this. On a **non-uniform**
grid `t_0 < t_1 < … < t_N`, `caputo_l1_nonuniform` uses the general L1 weights
built from the actual node spacings; it is validated independently in
`scirust-fractional` for exactness on linear functions on a non-uniform grid and
for agreement with the uniform scheme on a uniform grid. The two evaluators are
algebraically equal term-by-term under exactly uniform spacing but reach the
value by different floating-point paths, so bit-identity is input-specific, not
guaranteed.

The non-uniform operator is what makes adaptive stepping sound: an adaptive run
produces a genuinely non-uniform history, and applying a single `h` to it (as
`caputo_l1_uniform` would) is incorrect.

## 6. History storage

`HistoryBackend<D>` abstracts retained-sample storage:

- `CompleteUniformHistory` retains **every** accepted sample. It is the memory
  **oracle** for this model — the most faithful available discretization of the
  hereditary integral — at `O(N)` storage and `O(N)` per evaluation (`O(N^2)`
  over `N` steps).
- `BoundedShortMemoryHistory` retains only the most recent `W >= 2` samples. It
  is an explicit **approximation** (`O(W)` per evaluation) that must be selected
  deliberately.

Each retained sample is a typed `HistoryEntry` (coordinates, velocity, and its
true accepted parameter), so a geometric transport and a non-uniform memory law
can read the source point and parameter rather than only bare components.

## 7. Parallel transport of retained tangent vectors

Because the memory differentiates velocity *vectors* sampled at different points
of a curved manifold, one may transport each retained vector into the current
tangent space before differencing. `DiscreteConnectionTransport` advances every
retained vector by one Heun predictor-corrector step of the linear transport
equation

```text
dV^mu / dlambda = - Gamma^mu_{alpha beta} u^alpha V^beta
```

across each newly accepted segment. This accumulates along the actual accepted
polyline (transport is path-dependent under curvature). It is a **discrete
approximation**, not an exact bitensor propagator; its error grows with the
segment step and the number of segments, and it costs `O(D^3 N^2)`. The
production coordinate-memory path uses `IdentityHistoryTransport` (no transport),
which preserves the original chart-componentwise model exactly.

## 8. Proper-time and affine-parameter handling

`ParameterizationMode::AffineParameter` (default) advances an unconstrained
affine parameter. `NormalizedTimelikeProperTime { tolerance }` interprets the
step as proper time and *requires* `g(u,u)` to stay within `tolerance` of `-1`
under a `(-,+,+,+)` signature; a drift beyond tolerance is the typed error
`ProperTimeNormDrift`. **No automatic four-velocity renormalization is ever
performed** — drift is reported, never silently repaired. Separately,
`affine_trajectory_proper_time` estimates elapsed proper time along an affine
trajectory by first-order quadrature, and `proper_time_caputo_velocity_memory`
re-evaluates the memory of an already-computed trajectory against its own
non-uniform proper-time axis (a post-hoc diagnostic, never fed back into the
loop).

## 9. Curvature and field modulation as phenomenological hooks

`HistoryModulator<D>` multiplies each retained (and possibly transported)
velocity sample by a finite dimensionless scalar weight before the Caputo
stencil. Instances:

- `SchwarzschildKretschmannModulator`: `q = 1 + beta L^4 K`, with the
  Kretschmann scalar `K = 48 M^2 / r^6` and reference length `L`;
- `ReissnerNordstromFieldModulator`: `q = 1 + beta L^2 |F^2|`, with the
  electromagnetic field invariant `F_{mu nu} F^{mu nu} = 2 Q^2 / r^4`.

Both require valid exterior parameters and reject a non-finite or non-positive
weight. **When `beta = 0` the modulator returns exactly `1.0` and the pipeline
reproduces the unmodulated baseline bit-for-bit** (a bit-identity regression
test guards this). These are hand-chosen phenomenological reweightings — never
a consequence of general relativity, a quantum-gravity prediction, or an
experimentally derived law.

## 10. Fixed-step integration

`simulate_nonlocal_worldline_with_components` advances a fixed uniform step for
a fixed count, composing a `HistoryBackend`, `MemoryLaw`, `HistoryTransport`,
and `WorldlineStepper`. Two steppers:

- `SemiImplicitEulerStepper` (first order): `u_{n+1} = u_n + h a_n`,
  `x_{n+1} = x_n + h u_{n+1}`.
- `HeunPeceStepper` (second order, predict-evaluate-correct-evaluate): an Euler
  predictor, one acceleration evaluation at the predicted point against a
  provisional history, then the trapezoidal corrector.

The step index and the true accumulated parameter `n h` are both available to a
stepper via `StepperContext` (see §11 for why the true parameter matters).

## 11. Adaptive embedded Heun-Euler integration

`simulate_nonlocal_worldline_adaptive[_with_policy]` chooses its own non-uniform
step using the **embedded Heun-Euler pair**: the Euler predictor (order 1) and
Heun corrector (order 2) that one Heun step already computes form a first/second
order error pair at no extra acceleration evaluation. The memory force uses
`caputo_l1_nonuniform` against the accumulated non-uniform history. Pseudocode
(one accepted step):

```text
loop:
    (predicted, corrected) = heun_euler_step(state, accepted_accel, history, step)
    err = scaled_local_error_norm(predicted, corrected, tolerance)      # §13
    (decision, next_step) = control_step(err, step, &rejections, ...)   # §13
    match decision:
        Accept => break with (corrected, step, next_step)
        Retry  => step = next_step
push corrected into history at (current_parameter + step)
```

**Correct affine parameter (v2 fix).** `HeunPeceStepper` previously reconstructed
its provisional point's parameter as `step_index * config.step`, exact only
under uniform spacing. `StepperContext` now carries `current_parameter`, the
true accumulated affine parameter, and the provisional point is
`current_parameter + config.step`. This makes the Heun stepper sound under the
non-uniform spacing adaptive stepping produces, verified by a direct test: the
recorded provisional parameter equals `current_parameter + step` and is
independent of `step_index`.

## 12. Adaptive step-doubling integration

`simulate_nonlocal_worldline_adaptive_with_stepper[_policy]` drives
`SemiImplicitEulerStepper` (which has no natural embedded higher-order partner)
with classical **step-doubling**: one full step of size `h` versus two of size
`h/2`, the raw difference being the Richardson error estimate (the divisor
`2^p - 1 = 1` for the first-order method). It is the first adaptive path to
genuinely reuse `MemoryLaw`/`WorldlineStepper` via the non-uniform memory laws.

`HeunPeceStepper` is deliberately **not** offered here — not because of the old
parameter formula (now fixed), but because this controller's error estimate is
specialized to a first-order method, and the appropriate adaptive scheme for a
second-order method is the embedded pair of §11, which **is** adaptive
Heun-PECE. A step-doubling Heun variant would be a strictly inferior duplicate.

## 13. Scaled error control

Both controllers share one error norm and one control routine
(`adaptive_control`). The **componentwise scaled root-mean-square norm** is, for
each of the `D` coordinate and `D` velocity components,

```text
scale_i = abs_tol_i + rel_tol * max(|y_low_i|, |y_high_i|),
ratio_i = (y_high_i - y_low_i) / scale_i,
norm    = sqrt( (1 / 2D) * sum_i ratio_i^2 ).
```

with an `AdaptiveTolerance { relative, coordinate_absolute, velocity_absolute }`.
A step is accepted when `norm <= 1`; the next step is
`clamp(step * safety * norm^{-1/2}, min_step, max_step)` (growth capped, shrink
floored). This replaced an earlier unscaled `||Δx|| + ||Δu||` against one
absolute tolerance, which made acceptance depend excessively on the chart and on
component magnitudes. `control_step` enforces the retry budget
(`AdaptiveRejectionBudgetExhausted`) and the minimum step
(`AdaptiveMinimumStepExhausted`, a distinct error) identically for both
controllers. The reduction runs in a fixed component order, so the norm is
bit-for-bit reproducible.

This is the standard scaled RMS control used by production adaptive
Runge-Kutta codes (Hairer, Nørsett & Wanner, *Solving Ordinary Differential
Equations I*, §II.4). It improves scaling robustness. It is **not** a
geometrically invariant measure (see §14).

Configuration fields (`AdaptiveNonlocalConfig`):

| field | meaning |
|-------|---------|
| `alpha`, `coupling` | fractional order, memory coupling `kappa` |
| `initial_step`, `min_step`, `max_step` | step bounds |
| `error_tolerance` / `AdaptiveTolerance` | scalar (uniform) or three-field tolerance |
| `metric_norm_floor` | positive lower bound on `|g(u,u)|` |
| `target_affine_parameter` | integration endpoint |
| `max_accepted_steps`, `max_rejections_per_step` | typed-error budgets |

## 14. Coordinate dependence and limits of covariance

The memory is differentiated componentwise in the supplied chart, so it is
chart dependent by construction. Three v2 facts bound the covariance question
precisely:

1. The scaled error norm (§13) reduces the *controller's* sensitivity to the
   chart and to units, but it is itself componentwise and says nothing about the
   metric — it is not covariant.
2. `DiscreteConnectionTransport` reduces the memory's cross-chart disagreement
   (measurably, and shrinking under refinement in `coordinate_covariance`), but
   does not remove it: componentwise Caputo differentiation is not covariant
   even on transported inputs.
3. The metric-aware diagnostic of §15's companion, `timelike_state_error`,
   handles the indefinite signature correctly but is a post-hoc scalar
   temporal/spatial split at one chart point, not an invariant comparison of
   distant states.

No part of this layer is a covariant field theory.

**Metric-aware error diagnostic.** `timelike_state_error` decomposes a
coordinate error `delta` relative to a timelike observer `u` into
`temporal = |g(delta,u)| / sqrt(-g(u,u))` and
`spatial = sqrt(g(P delta, P delta))`, with `P` the §4 projector. `P delta` is
spacelike, so `spatial` is a genuine non-negative length; a non-timelike `u` is
rejected. In flat spacetime with a static observer this is exact: a purely
spatial `(0,3,4,0)` gives `temporal = 0`, `spatial = 5` bit-for-bit.

## 15. Exact flat-spacetime validation oracle

Parallel transport on a curvature-free manifold is path-independent, so a vector
transported in Cartesian Minkowski (where the connection vanishes) is unchanged.
`exact_cylindrical_minkowski_transport` uses this with **no discretization**:
convert to Cartesian, transport trivially, convert back at the destination.
`transport_vector_along_polyline` exposes `DiscreteConnectionTransport` over an
explicit waypoint list; `exact_transport_convergence` shows the discrete error
shrinking by ≈4 per waypoint doubling (second order), from `~3.5e-5` at 4 to
`~3.5e-8` at 128 waypoints — a **numerical oracle** validation against a
known-exact answer. Exact only for this flat chart pair, within a simply
connected region.

## 16. Exact Schwarzschild circular-orbit validation oracle

Along a circular equatorial geodesic orbit (fixed `r`, `theta = pi/2`, constant
four-velocity), Schwarzschild's Christoffel symbols are constant, so the
transport equation is a linear constant-coefficient ODE with solution
`V(lambda) = exp(-lambda A) V(0)` for a fixed generator `A`.
`exact_schwarzschild_circular_orbit_transport` evaluates this by a deterministic
4x4 matrix exponential, validated against two exact conservation laws (`g(V,V)`
and `g(V,u)` constant under any metric-compatible transport) and against
`DiscreteConnectionTransport`'s second-order convergence
(`schwarzschild_orbit_transport`). This is the **only** exact oracle for a
*curved* background here. **No exact reference is currently implemented for a
general curved path** (eccentric, inclined, or non-geodesic); the flat and
circular-orbit oracles are the two special cases, and neither extends.

## 17. Numerical convergence experiments

`experiments/nonlocal-relativity-v2/adaptive_convergence` compares both adaptive
controllers against a fine fixed-step reference of the matching method
(`h = 5e-4`), Schwarzschild `M = 1`, near-circular `r0 = 10`, `alpha = 0.55`,
`kappa = 0.02`, target `0.8`. Representative endpoint coordinate errors:

| tolerance | embedded Heun-Euler (steps, err) | step-doubling (steps, err) |
|-----------|----------------------------------|-----------------------------|
| 1e-6 | 17, 5.48e-9 | 17, 2.92e-6 |
| 1e-7 | 19, 4.68e-9 | 17, 2.92e-6 |
| 1e-8 | 56, 1.89e-9 | 40, 1.20e-6 |
| 1e-9 | 175, 3.52e-10 | 124, 3.39e-7 |

The error decreases monotonically as the tolerance tightens and the step count
rises. The embedded second-order pair reaches ~`3.5e-10`; the first-order
step-doubling path ~`3.4e-7`, consistent with their orders. These are
**fine-grid reference** comparisons (self-consistency), not model validation.
`run_convergence_study` additionally performs `h`, `h/2`, `h/4`
**self-convergence** for the fixed-step methods.

Reproduce:

```bash
cargo run --release -p nonlocal-relativity-experiments --bin adaptive_convergence
```

## 18. History-retention comparison

The step-doubling controller computes a midpoint and an endpoint per accepted
step but retains only the endpoint (`HistoryRetention::EndpointOnly`).
`RefinedAcceptedHistory` additionally retains the midpoint at its true affine
parameter. The experiment `history_retention` compares both against a fine
reference:

| tolerance | strategy | steps | retained | op-proxy | coord err |
|-----------|----------|-------|----------|----------|-----------|
| 1e-8 | endpoint_only | 40 | 41 | 861 | 1.197e-6 |
| 1e-8 | refined | 40 | 81 | 1681 | 1.197e-6 |
| 1e-9 | endpoint_only | 124 | 125 | 7875 | 3.387e-7 |
| 1e-9 | refined | 124 | 249 | 15625 | 3.387e-7 |

**Result (against the initial hypothesis).** Retaining midpoints leaves the
accepted-step count unchanged and the endpoint error identical to ~4 significant
figures, while roughly doubling the retained sample count and the op-count
proxy. Denser history does not improve endpoint accuracy here: the memory force
is a small perturbation and the endpoint error is dominated by the integrator's
truncation order, not the memory-quadrature density.

**Decision.** Keep `EndpointOnly` as the default; expose
`RefinedAcceptedHistory` only as an explicit research option. Structural
invariants (strict parameter ordering, no duplicates, true-midpoint recording,
exact retained counts with no leakage from rejected trials) are pinned by
`tests/history_retention.rs`.

## 19. Complexity and performance

`complexity_scaling` reports a **deterministic operation-count proxy** (sum of
retained sample counts over accepted evaluations — the exact Caputo
history-sample touch count) at doubling fixed step counts. Wall-clock time is
not used, so no asymptotic claim rests on timing.

| pipeline | proxy(50) | proxy(100) | proxy(200) | proxy(400) | ratio → |
|----------|-----------|-----------|-----------|-----------|---------|
| complete raw | 1326 | 5151 | 20301 | 80601 | 3.97 (`O(N^2)`) |
| bounded (W=16) | 696 | 1496 | 3096 | 6296 | 2.03 (`O(N*W)`) |
| discrete transport | 1326 | 5151 | 20301 | 80601 | 3.97 (`O(N^2)`) |

The measured ratios match the implementation-derived complexity: complete raw
memory is `O(N^2)` (ratio → 4), bounded short memory is `O(N*W)` (ratio → 2 once
`N > W`), and discrete transport has the same `O(N^2)` touch count as complete
memory with an additional `O(D^3)` Christoffel contraction per touch. The proxy
for complete memory equals the triangular number `(N+1)(N+2)/2` exactly.

**Validated micro-optimizations** in v2 are limited to sharing one error-norm
and control routine between the two adaptive controllers (removing duplicated
code, not changing arithmetic); the exact-history oracle path was deliberately
left unaltered, and no floating-point reduction order was changed on the default
paths. Incremental-Caputo recurrences and buffer reuse were considered but not
applied to the oracle path, to preserve its bit-for-bit determinism.

## 20. Limitations

- The memory is coordinate dependent; transport and modulation reduce but do
  not remove this.
- The exact oracles cover only flat spacetime and Schwarzschild circular
  equatorial geodesic orbits; there is no exact reference for a general curved
  path.
- `Kerr`'s connection is finite-difference, carrying a small step-dependent
  truncation error unlike the analytic backgrounds.
- The adaptive schemes are standard but simple (embedded Heun-Euler and
  first-order step-doubling); no dense output, event handling, or higher-order
  embedded pairs.
- The geometric error diagnostic is a scalar temporal/spatial split at one
  chart point, not a full tetrad projection or an invariant distant comparison.
- No result is a physical validation.

## 21. Future work

- A general curved-path exact bitensor propagator (neither special-case oracle
  extends).
- A full local-orthonormal-frame (tetrad) error projection.
- An exact analytic Kerr connection; Kerr-specific orbit/transport constructions.
- Higher-order or dense-output adaptive schemes.
- Modulators from other invariants or backgrounds.

**Modified field equations are permanently out of scope** — excluded by the
project's scientific boundary, not deferred.

## 22. Explicitly forbidden interpretations

This layer must never be described as, or implied to be, any of: a modification
of the Einstein field equations; fractional Einstein equations; a modification
of the Einstein tensor or stress-energy tensor; a computation of matter
backreaction on curvature; a complete covariant theory; evidence of new physics;
or an experimentally validated theory. `DiscreteConnectionTransport` is a
discrete approximation, not a proof of covariance. The modulators are
phenomenological reweightings, not physical laws. The exact transport oracles
are valid only for their stated special cases. The adaptive integrators are
standard numerical techniques applied to an ordinary state equation; they carry
no physical claim beyond what the discretization provides.

---

### References

- E. Hairer, S. P. Nørsett, G. Wanner, *Solving Ordinary Differential Equations
  I: Nonstiff Problems*, Springer Series in Computational Mathematics — the
  standard reference for embedded Runge-Kutta error estimation, scaled
  root-mean-square local-error control, and step-doubling / Richardson error
  estimation used in §§11–13.

(No other bibliographic references are asserted. The Caputo derivative and its
L1 discretization are used as standard, widely-documented constructions;
specific claims here are backed by the platform's own tests and experiments,
cited inline by file, not by external result.)

### Reproducing every quantitative claim

```bash
# tests (three crates) and doctests
cargo test --locked -p scirust-fractional -p scirust-relativity -p scirust-nonlocal-relativity
cargo test --locked -p scirust-nonlocal-relativity --doc

# experiments (deterministic CSV to stdout)
export NLR_EXPERIMENT_COMMIT="$(git rev-parse HEAD)"
cargo run --release -p nonlocal-relativity-experiments --bin history_retention
cargo run --release -p nonlocal-relativity-experiments --bin adaptive_convergence
cargo run --release -p nonlocal-relativity-experiments --bin complexity_scaling

# example scenarios (deterministic CSV to stdout)
cargo run --release -p scirust-nonlocal-relativity --example exact_transport_convergence
cargo run --release -p scirust-nonlocal-relativity --example schwarzschild_orbit_transport
cargo run --release -p scirust-nonlocal-relativity --example coordinate_covariance
```

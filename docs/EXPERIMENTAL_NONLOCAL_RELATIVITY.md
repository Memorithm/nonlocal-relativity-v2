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
- deviations from an uncoupled geodesic baseline in a fixed chart.

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

## Roadmap

### 1. Current Worldline-Memory Model

The current crate implements a fixed-background, test-particle, coordinate
Caputo-memory modification of the worldline equation. It exposes deterministic
diagnostics and explicit failure modes.

### 2. Future Covariant Kernel Research

Future research may study kernels defined with bitensors, parallel transport,
proper-time history, or other covariant constructions. These are research
directions only and are not established physics in this crate.

### 3. Hypothetical Field-Equation Work

Any future field-equation investigation would require a separate mathematical
and numerical contract, independent validation, and clear distinction from
established general relativity. This crate does not implement such work, and
this roadmap item is not a claim that fractional field equations are
established physics.

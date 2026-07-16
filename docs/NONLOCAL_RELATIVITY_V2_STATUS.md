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

## Validations Performed

- `cargo fmt --all -- --check` clean.
- `cargo test --locked -p scirust-fractional -p scirust-relativity -p scirust-nonlocal-relativity`:
  all tests and doctests passing (exact counts in the phase commit messages
  and PR description).
- `cargo clippy --locked -p scirust-fractional -p scirust-relativity -p scirust-nonlocal-relativity --all-targets -- -D warnings`
  clean.
- All crate examples (`schwarzschild_memory`, `convergence_study`,
  `coordinate_covariance`, `curvature_modulated_memory`) run to completion
  and produce deterministic CSV output.
- Bit-for-bit regression: every Phase 1/2 test file is unmodified and passes
  unchanged; Phase 3/4 additions include explicit bit-identity tests for the
  compatibility paths (`beta = 0`, identity transport, affine mode).
- A dedicated test mechanically scans this crate's own source for item
  declarations (`struct`/`enum`/`trait`/`fn`) whose name suggests a modified
  field equation, Einstein tensor, or stress-energy structure, and fails if
  one is found.

## Complexities (as actually implemented)

| Component | Cost |
|---|---|
| Raw coordinate memory (`CaputoCoordinateMemory` + `IdentityHistoryTransport`) | `O(D * N^2)` over `N` fixed steps (unchanged from Phase 1/2) |
| Bounded short memory (`BoundedShortMemoryHistory`) | `O(D * N * W)` for window `W` (unchanged) |
| Discrete parallel transport (`DiscreteConnectionTransport`) | `O(D^3 * N^2)`: `O(N)` transported vectors per accepted step (`O(N^2)` total), each a Christoffel contraction (`O(D^3)`) |
| Curvature modulation (`ModulatedCaputoCoordinateMemory`) | Adds `O(1)` work per retained sample per evaluation on top of whichever transport/backend it wraps |

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

## Limitations

- The memory kernel remains coordinate-dependent; transport and modulation
  change this quantitatively (measurably smaller cross-chart disagreement,
  shrinking under refinement) but do not eliminate it.
- `DiscreteConnectionTransport` is a discrete, segment-by-segment
  approximation, not an exact bitensor propagator, and its cost is
  asymptotically worse than the coordinate-memory baseline.
- `NormalizedTimelikeProperTime` validates but does not adapt the step; drift
  beyond tolerance is a hard error, not a corrected trajectory.
- `SchwarzschildKretschmannModulator`'s `beta` and reference length are free,
  uncalibrated phenomenological parameters specific to the Schwarzschild
  exterior chart.
- No adaptive stepping, event handling, error estimation, or history
  compression is included anywhere in the crate.

## Future Work (not implemented here)

- Exact analytic bitensor parallel propagators, as a replacement for the
  discrete segment-by-segment transport.
- Proper-time history sampled at its own adaptive resolution.
- Curvature modulators for backgrounds other than Schwarzschild, or built
  from invariants other than the Kretschmann scalar.
- Any investigation of modified field equations — explicitly out of scope
  for this crate; see the next section.

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

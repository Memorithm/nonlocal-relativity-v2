# SciRust Nonlocal Relativity — Platform Roadmap

This document states the long-horizon vision for turning the current
experimental worldline simulator into a broad, rigorous, open-source platform
for relativistic geometry and gravitational physics, and — just as importantly
— it draws a hard line between **what is built and validated today** and **what
is planned**. It is a roadmap, not a status report; the authoritative
per-feature status lives in
[`NONLOCAL_RELATIVITY_V2_STATUS.md`](NONLOCAL_RELATIVITY_V2_STATUS.md),
[`EXPERIMENTAL_NONLOCAL_RELATIVITY.md`](EXPERIMENTAL_NONLOCAL_RELATIVITY.md),
and the crate `README`s.

The objective is not to advance a new theory of gravity. It is to build the
most rigorous, deterministic, auditable, and reusable open-source foundation
for relativistic geometry and gravitational-physics computation that we can,
and to grow it one *validated* increment at a time. Scientific integrity takes
precedence over novelty. Negative results are kept, not discarded.

## Non-negotiable engineering and scientific principles

These constraints govern every increment and are not traded away for scope:

- 100% Rust, zero FFI, zero `unsafe` unless mathematically unavoidable (and
  none is used today).
- No `TODO`, `FIXME`, placeholder, or stub implementations in the hardened
  crates; no panic-based control flow; typed errors everywhere.
- Deterministic execution — no RNG and no wall-clock in library or experiment
  code paths — and bit-for-bit reproducibility wherever mathematically
  possible.
- CI stays green: `cargo fmt --check`, `cargo clippy -D warnings`, tests,
  doctests, and a forbidden-marker scan over the hardened sources.
- Every algorithm is referenced to the literature; every numerical claim is
  backed by a reproducible experiment or test against a known oracle.
- No claim of new physics without both a mathematical derivation and
  observational validation. Established general relativity, phenomenological
  models, and speculative extensions are labelled as such and never blurred.

## The classification discipline

Every result the platform produces belongs to exactly one of the following
categories, and each is labelled explicitly in code, docs, and experiment
headers. Blurring these categories is treated as a defect.

1. **Exact mathematical result** — a closed-form identity (e.g. the
   Kretschmann scalar `48 M^2 / r^6` of Schwarzschild, or `R = 4 Lambda` for a
   4D maximally symmetric spacetime).
2. **Numerical implementation of established physics** — a deterministic
   computation whose target is textbook general relativity, validated against
   an exact oracle to a stated tolerance (e.g. the finite-difference curvature
   engine).
3. **Numerical approximation** — a scheme with a disclosed, quantified error
   (e.g. finite-difference truncation, bounded-memory truncation).
4. **Phenomenological model** — a deterministic but uncalibrated construction
   with free parameters and no claim of physical law (e.g. the curvature-memory
   modulation hooks in the experimental worldline layer).
5. **Speculative model** — an extension explored for its structure, explicitly
   not asserted to describe nature.
6. **Empirically validated physics** — a claim checked against observation.
   The platform currently makes **none** of these.

## The six-layer vision

| Layer | Scope | Status |
|------|-------|--------|
| 1 — Geometry Core | manifolds, metrics, tetrads, connections, curvature tensors, geodesics, parallel transport, bitensors, world function, geometry diagnostics | **partially delivered** (see below) |
| 2 — Covariant Gravity Workbench | symbolic actions, variational calculus, automatic field-equation derivation, weak-field / PPN / cosmological limits, stability & ghost checks | planned |
| 3 — Numerical Relativity | linear perturbations, self-force, EMRI; then ADM/BSSN, constraint damping, AMR, wave extraction | planned |
| 4 — Gravitational Memory Lab | standard / Christodoulou / fractional memory, observer and detector response | partially explored in the experimental worldline layer (phenomenological) |
| 5 — Astrophysical Inference | waveform generation, noise models, likelihood, MCMC / nested sampling, matched filtering | planned |
| 6 — Relativistic Navigation | proper time, Shapiro delay, redshift, GNSS corrections, filtering across Earth/Moon/Mars/Sun/deep space | planned |

"Planned" means exactly that: not started, or present only as an experimental,
clearly-labelled prototype. The table is revised only when a layer's status
genuinely changes, and only downward-conservatively.

## Current status — what is actually built and validated

The following exists today, is tested, and ships in a green CI:

**Layer 1 (Geometry Core), delivered subset:**

- Trait-based backgrounds: `Metric<D>` and `Connection<D>` over const-generic
  dimension, with deterministic metric inversion (`invert_metric`) and typed
  errors (`RelativityError`).
- Backgrounds with analytic connections: Minkowski, Schwarzschild,
  Reissner–Nordström, de Sitter, and anti-de Sitter; Kerr with a
  finite-difference connection whose truncation cost is measured and disclosed.
- A shared lapse-metric helper (`static_spherical`) for static, spherically
  symmetric `f(r)` spacetimes, reused by the de Sitter / anti-de Sitter
  backgrounds.
- **Curvature engine** (`CurvatureTensors`): Riemann, Ricci, Ricci scalar,
  Einstein, and Kretschmann tensors from any `Metric + Connection` background
  by central finite differences of the Christoffel symbols, with typed errors
  and no non-finite results.
- **Exact-oracle validation** of the curvature engine (category 2 against
  category 1):
  - Minkowski: every curvature component is *exactly* zero.
  - Schwarzschild: Ricci-flat, and `K = 48 M^2 / r^6` (e.g. exactly
    `1.8310546875e-4` at `M = 1, r = 8`) to finite-difference tolerance.
  - de Sitter / anti-de Sitter: `R_(mu nu) = Lambda g_(mu nu)`,
    `R = 4 Lambda`, `G_(mu nu) = -Lambda g_(mu nu)`, `K = 8 Lambda^2 / 3`,
    and the maximally-symmetric Riemann form
    `R_(abcd) = (Lambda / 3)(g_(ac) g_(bd) - g_(ad) g_(bc))`.
  - Riemann index symmetries and the first Bianchi identity.
  - A `curvature_invariants` experiment that reports invariants against these
    oracles and sweeps the finite-difference step to expose the `O(h^2)`
    truncation-vs-roundoff trade-off.
- Geodesic integration (`GeodesicSystem`) compatible with `scirust-sim`.

**Layer 4 (memory), experimental and phenomenological only:** the
`scirust-nonlocal-relativity` crate integrates fractional-memory *test-particle*
worldlines on a fixed background. It does **not** modify the field equations
and claims **no** physical validity; see its status document. It is listed here
as an existing prototype, not a delivered platform layer.

Everything not in this section is not yet built.

## Increment plan

The platform grows by small, individually-validated pull requests. Each lands
only with tests or experiments against an oracle, green CI, and honest labels.
The near-term ordering within Layer 1:

1. **Curvature core + de Sitter / anti-de Sitter** — *this increment.* Riemann
   through Kretschmann, validated against exact maximally-symmetric and
   Schwarzschild oracles.
2. Coordinate-independence / invariant diagnostics: confirm scalar invariants
   agree across charts of the same geometry to tolerance.
3. Tetrads and parallel transport as a reusable engine, with holonomy and
   geodesic-deviation (Jacobi) checks against closed forms in maximally
   symmetric spacetimes.
4. Bitensors and Synge's world function on backgrounds with known expansions.
5. FLRW cosmological background with its exact curvature oracle.

Layers 2–6 begin only after Layer 1 is broad and solid. Each will open with a
design note fixing its oracles and its category labels before any code lands.

## What this platform will not do

- It will not assert modified-gravity or nonlocal-gravity physics as
  established. Such models remain in clearly-labelled experimental layers.
- It will not present a numerical approximation as an exact result, or a
  phenomenological hook as a physical law.
- It will not hide a limitation. Disclosed error (finite-difference truncation,
  bounded memory) is part of the result, reported alongside it.

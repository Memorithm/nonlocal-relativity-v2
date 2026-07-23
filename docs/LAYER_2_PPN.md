# Layer 2 — PPN parameter extraction (γ, β): design & conventions

This is the architecture-decision record for the second Layer 2 increment:
oracle-backed extraction of the Eddington–Robertson PPN parameters **γ** and
**β** from static, spherically symmetric, asymptotically flat, weak-field metrics
given in a PPN-compatible **isotropic** radial coordinate. It is written before
the code, per the platform rule.

It extends the Layer 2 Covariant Gravity Workbench
([`LAYER_2_COVARIANT_GRAVITY.md`](LAYER_2_COVARIANT_GRAVITY.md)) and **reuses**
the Layer 1 geometry core (`Metric`, `IsotropicSchwarzschild`), the experiment
framework, the determinism guarantees, and the scientific-category discipline.
It does **not** restart or redesign any validated component, and it introduces
**no** symbolic-algebra system.

## 1. Scientific target and conventions (signs stated explicitly)

Signature is `(-,+,+,+)` throughout (the crate convention: `Minkowski` gives
`g_00 = -1`). For a static, spherically symmetric metric in a PPN-compatible
isotropic radial coordinate `ρ`, with Newtonian compactness

```text
U(ρ) = G M / ρ ,        0 < U ≪ 1  (weak field),
```

the Eddington–Robertson weak-field expansions are

```text
g_00(ρ) = -1 + 2U - 2β U^2 + O(U^3),
g_ij(ρ) = A(ρ) δ_ij ,   A(ρ) = 1 + 2γ U + O(U^2).
```

`A(ρ)` is the **spatial conformal factor** (the isotropic spatial metric is
conformally flat: `g_ρρ = g_θθ/ρ^2 = g_φφ/(ρ^2 sin^2θ) = A`). For **General
Relativity in isotropic Schwarzschild coordinates**,

```text
g_00 = -((1 - M/2ρ)/(1 + M/2ρ))^2 = -1 + 2U - 2U^2 + O(U^3)  ⇒  β = 1,
A    = (1 + M/2ρ)^4               = 1 + 2U + (3/2)U^2 + O(U^3)  ⇒  γ = 1,
```

(hand-derived series, documented alongside the oracle). The extractor must
asymptotically recover `γ = β = 1`. Signs are never inferred silently: the
estimators below are written for exactly this signature and expansion.

## 2. Finite-radius effective estimators

Solving the two expansions for the leading coefficients gives the per-radius
effective estimators

```text
γ_eff(ρ) = (A(ρ) - 1) / (2U)                    = γ + O(U),
β_eff(ρ) = -(g_00(ρ) + 1 - 2U) / (2U^2)         = β + O(U).
```

These are **diagnostic data**, not the answer: each carries an `O(U)` (γ) or
`O(U)` (β, after the `U^2` division) contamination from the higher-order terms.
They are exposed per radius (`FiniteRadiusEstimate { radius, compactness,
value }`) and never returned as the asymptotic result.

## 3. Asymptotic extraction (deterministic polynomial least-squares)

The asymptotic value is the `U → 0` intercept of the effective estimator. We fit

```text
γ_eff(U) ≈ γ + c_1 U + c_2 U^2 + … + c_n U^n ,
β_eff(U) ≈ β + d_1 U + d_2 U^2 + … + d_n U^n ,
```

by ordinary least squares and read the intercept. Implementation choices, all to
keep the fit deterministic and auditable:

- **Normalized fit variable** `s = U / U_max ∈ (0, 1]` (the intercept `s → 0` is
  the same `U → 0` limit). Normalization keeps the Vandermonde/normal-equations
  system well scaled; without it the moments span `1 … U_max^{2n}` and become
  hopelessly ill-conditioned for `U ≪ 1`.
- A **small self-contained dense solver** for the `(n+1)×(n+1)` normal equations
  `(VᵀV) c = Vᵀy`: Gaussian elimination with partial pivoting, deterministic row
  order and accumulation order, no allocation-order dependence, no randomness, no
  hidden regularization. (`invert_metric`/`determinant` are const-generic and
  cannot take a runtime degree, so a dedicated small solver is used — no new
  dependency, matching the crate's self-contained linear-algebra style.)
- **Configurable degree** within a small validated range `1 ≤ n ≤ 4`.
- **Rejections** (typed `PpnError`): fewer samples than `n + 1`
  (underdetermined), non-finite inputs, a (near-)singular pivot (`SingularFit` /
  `IllConditionedFit`), an unsupported degree.
- A reported **conditioning indicator** (the scaled minimum pivot magnitude);
  ill-conditioned systems are rejected rather than silently trusted.

## 4. Radial domain

`PpnDomain { mass_scale, radius_min, radius_max, sample_count, sampling }` with a
deterministic `sampling` policy — uniform in compactness `U`, logarithmic in `ρ`,
or explicit radii. Validated: `mass_scale` finite and positive; `radius_min <
radius_max`; the largest compactness `U_max = M/radius_min` is below a weak-field
cap; every radius outside the metric's validity domain (caught by the metric's
own finiteness/validity checks); deterministic sample ordering; at least `n + 1`
samples for the requested degree. No single radius interval is hardcoded as
universally valid.

## 5. Coordinate-validity contract (mandatory)

PPN coefficients are coordinate-sensitive, so the extractor never accepts a raw
`Metric` and silently reads its radial coordinate as isotropic. Extraction is
defined over a dedicated trait

```text
StaticIsotropicMetric {
    mass_scale(&self) -> f64,
    g_tt(&self, radius) -> Result<f64, PpnError>,
    spatial_conformal_factor(&self, radius) -> Result<f64, PpnError>,
}
```

Two ways to obtain one:

- **Directly** — the synthetic oracles implement it from their defining
  coefficients (isotropic by construction).
- **Via an explicit adapter** over a spherical `Metric<4>` that **checks
  conformal flatness** at each radius: `g_ρρ`, `g_θθ/ρ^2`, `g_φφ/(ρ^2 sin^2θ)`
  must agree, else `PpnError::NonIsotropicCoordinates`. `IsotropicSchwarzschild`
  passes this check; **areal Schwarzschild fails it** (`g_rr = 1/f` while
  `g_θθ/r^2 = 1`) and is rejected — the mandatory negative oracle. It is thus
  impossible to *silently* apply the isotropic PPN formulas to areal coordinates.

## 6. Error and contamination characterization

The returned uncertainty is not decorative. Several deterministic diagnostics are
computed and combined conservatively into an **estimated numerical uncertainty**
(explicitly *not* a rigorous bound):

- fit residual norm of the selected model;
- **window sensitivity** — re-extract on the weaker-field sub-window and compare
  intercepts;
- **order sensitivity** — compare adjacent valid extrapolation degrees;
- **resolution sensitivity** — re-extract on a deterministic sub-sample;
- **oracle error** — absolute/relative error when an exact answer is known;
- **conditioning** — the solver's conditioning indicator.

Names use "estimated uncertainty" / "sensitivity" / "fit diagnostic", never
"proof" or "certified bound".

## 7. Oracle hierarchy

1. **Exact synthetic PPN metric** — `g_00 = -1 + 2U - 2β⋆U^2`,
   `A = 1 + 2γ⋆U`, exactly. For these the effective estimators are already
   *constant* (`γ_eff ≡ γ⋆`, `β_eff ≡ β⋆`), so the extractor recovers the injected
   pair to machine precision. Pairs include `(1, 1)` and non-GR pairs. Isolates
   extraction correctness from contamination.
2. **Contaminated synthetic** — add known `a_3 U^3 + a_4 U^4` to `g_00` and
   `b_2 U^2 + b_3 U^3` to `A`. Finite-radius estimators are biased; the
   extrapolation converges to the injected values; a weaker-field window improves
   the result; the reported sensitivity tracks the actual error.
3. **Isotropic Schwarzschild** — the exact `IsotropicSchwarzschild` metric
   functions (not the truncated series the extractor uses) → `γ = β = 1`.
4. **Invalid coordinates** — areal-coordinate `Schwarzschild` through the adapter
   is rejected (`NonIsotropicCoordinates`), never silently trusted. Mandatory.

## 8. Scope and explicit exclusions

**In scope:** static · spherically symmetric · asymptotically flat · isotropic
(or explicitly PPN-compatible) coordinates · γ and β · finite-radius estimates ·
asymptotic extrapolation · convergence/contamination diagnostics · the synthetic
and isotropic-Schwarzschild oracles.

**Excluded for now** (documented so the boundary is explicit): preferred-frame
and nonconservative PPN parameters; the full ten-parameter formalism; rotating,
time-dependent, non-spherical/multipolar, or cosmological metrics; automated
arbitrary-coordinate conversion; observational likelihoods.

## 9. Category labels and honesty

- The estimators and the extrapolation are a **numerical approximation**; the
  extrapolation is never labelled exact.
- The synthetic-metric recoveries and the isotropic-Schwarzschild convergence are
  **validation oracles** for the *implementation* — passing the GR oracle
  validates the code, not an alternative theory.
- The exact synthetic coefficients and the isotropic-Schwarzschild series are
  **exact analytical** inputs.
- Estimated uncertainty is **not** automatically a rigorous bound; observational
  consistency would require a separate data-analysis layer (Layer 5), out of
  scope here.

## 10. Placement

The `ppn` module lives in `scirust-relativity` (as linearized gravity does; a
dedicated `scirust-gravity` crate remains an option if the Layer 2 surface
grows), with a dedicated `PpnError`, a `ppn`-scoped module tree, a
`linearized_gravity`-style experiment, and `criterion` benches — all reusing the
existing conventions.

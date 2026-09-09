# BSSN grid — current implementation contract

This note records the **current code contract** of `scirust-relativity/src/bssn_grid.rs` so that future work, audits, and autonomous agents do not reconstruct the solver from older prose that predates the live-gauge increments.

It is documentation of implemented numerical-relativity machinery, not a new physics claim.

## Flat state layout: 24 scalars per grid point

`COMPONENTS_PER_POINT` is currently `24`.

The structure-of-arrays slots are:

| slots | field | count |
|---|---|---:|
| `0` | conformal factor `phi` | 1 |
| `1..7` | six independent components of `gammatilde_ij` | 6 |
| `7` | mean curvature `K` | 1 |
| `8..14` | six independent components of `Atilde_ij` | 6 |
| `14..17` | conformal connection `Gammatilde^i` | 3 |
| `17` | lapse `alpha` | 1 |
| `18..21` | shift `beta^i` | 3 |
| `21..24` | Gamma-driver auxiliary `B^i` | 3 |

Total: `1 + 6 + 1 + 6 + 3 + 1 + 3 + 3 = 24` scalars per grid point.

Symmetric tensors store only their six independent components, so symmetry is represented structurally rather than repaired after a round trip.

## Slicing actually implemented

`BssnSlicing` currently exposes two modes:

- `Prescribed` — `d_t alpha = 0`; this remains the default.
- `OnePlusLog` — live 1+log slicing. With non-zero shift the implementation also includes lapse advection.

The lapse is a stored field in both modes. Selecting `Prescribed` freezes its right-hand side; it does not remove the lapse from the state layout.

## Shift conditions actually implemented

`BssnShiftCondition` currently exposes:

- `Prescribed` — freezes both `beta^i` and `B^i`; this remains the default.
- `GammaDriver { eta }` — evolves the hyperbolic Gamma-driver auxiliary system through `gamma_driver_rhs`.

The driver is fed the **full** conformal-connection rate used by the grid right-hand side. A live shift is therefore an explicit opt-in capability, not an absent future feature.

## Spatial dimensionality

The grid implementation is generic over `UniformGrid<D>` and has validation coverage for `D = 1`, `D = 2`, and `D = 3` in the existing experiment suite.

Every point still carries full three-dimensional BSSN tensor/vector data. When `D < 3`, axes outside the varying grid dimensions have zero spatial derivatives by construction.

Current restrictions remain important:

- periodic boxes;
- equal spacing on all active axes (`AnisotropicGrid` is rejected);
- weak, smooth fields in the validated scenarios;
- no adaptive mesh refinement;
- no puncture/excision machinery;
- no radiative outer boundary;
- no waveform-extraction pipeline;
- no claim that measured stability constitutes an analytic proof of strong hyperbolicity.

## Genuine BSSN principal part

The grid path uses `conformal_ricci_from_derivatives` with the **evolved** `Gammatilde^i`, rather than substituting the generic Ricci tensor of `gammatilde_ij` into BSSN variables.

This distinction is numerically load-bearing: the repository history records that the generic-Ricci path retained an ADM-like principal part and was unstable under refinement, while the genuine BSSN form restored the intended principal structure and the measured refinement behavior.

## Contract for future changes

Any future change to the BSSN grid should update this note when it changes one of the following:

1. `COMPONENTS_PER_POINT` or any slot constant;
2. `BssnSlicing` variants or equations;
3. `BssnShiftCondition` variants or equations;
4. supported grid dimensionality or spacing constraints;
5. the source of the conformal Ricci principal part;
6. the list of explicitly validated strong-field/boundary capabilities.

A future documentation hardening should also reconcile older module-level comments in `bssn_grid.rs` that still describe the earlier 17-component, prescribed-gauge stage. This note is the current authoritative inventory until that source-level prose is rewritten together with a focused code review of the full file.

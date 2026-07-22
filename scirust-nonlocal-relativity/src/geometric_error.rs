//! Optional research diagnostic: a Lorentzian-correct decomposition of a
//! coordinate state-error vector relative to a timelike observer.
//!
//! [`crate::scaled_local_error_norm`] (the adaptive step controller's norm) is
//! componentwise and chart-dependent: it treats the state as a point in a
//! Euclidean `R^{2D}` and says nothing about the indefinite spacetime metric.
//! That is the right, cheap tool for *step-size control*, but it is not a
//! geometrically meaningful error measure, and a Lorentzian-signature metric
//! must never be presented as an ordinary positive Euclidean norm.
//!
//! This module provides the geometric alternative for a coordinate
//! displacement `delta` (a small difference of two nearby states' coordinates,
//! treated as an approximate tangent vector at the reference point), decomposed
//! relative to a timelike four-velocity `u`:
//!
//! - the **temporal** part `|g(delta, u)| / sqrt(-g(u,u))` — the size of the
//!   component of `delta` along the observer's proper-time direction;
//! - the **spatial** part `sqrt(g(P delta, P delta))`, where
//!   `P^mu_nu = delta^mu_nu - u^mu u_nu / g(u,u)` is the crate's existing
//!   projector onto the subspace orthogonal to `u`. Because `P delta` is
//!   orthogonal to the timelike `u`, it is spacelike, so `g(P delta, P delta)`
//!   is non-negative and its square root is a genuine positive length.
//!
//! This handles the indefinite metric correctly (the temporal and spatial
//! parts are separated *by the metric*, not lumped into one Euclidean sum), and
//! it is exact in closed form — no tetrad construction is needed for the scalar
//! temporal/spatial magnitudes. It remains chart-dependent in the weaker sense
//! that `delta` is a coordinate difference and `u`, the metric are evaluated at
//! one chart point; it is **not** a full invariant comparison of two distant
//! states, and it is deliberately a post-hoc diagnostic, **not** wired into the
//! adaptive step controllers.
//!
//! A general local-orthonormal-frame (tetrad) projection, which would express
//! `delta` in an observer's own orthonormal basis component by component, is a
//! larger construction and is left as future work; the scalar temporal/spatial
//! split delivered here is the part that is exactly verifiable now (see the
//! flat-spacetime tests in `tests/geometric_error.rs`).

use crate::{
    NonlocalRelativityError, NonlocalResult, lower_index, metric_contraction, validate_scalar,
};

/// Metric-aware decomposition of a coordinate state-error vector relative to a
/// timelike observer. See the module documentation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TimelikeStateError {
    /// `|g(delta, u)| / sqrt(-g(u,u))`: the magnitude of the error along the
    /// observer's proper-time direction. Non-negative.
    pub temporal: f64,
    /// `sqrt(g(P delta, P delta))`: the metric length of the error projected
    /// orthogonal to `u`. Non-negative (the projected vector is spacelike).
    pub spatial: f64,
    /// `g(P delta, u)`: the metric inner product of the projected error with
    /// `u`. Zero in exact arithmetic (the projection is orthogonal to `u`);
    /// exposed as a rounding-level self-consistency residual.
    pub orthogonality_residual: f64,
}

/// Decompose a coordinate state-error vector `delta` at a point with metric
/// `metric` and timelike four-velocity `four_velocity` into metric-aware
/// temporal and spatial magnitudes.
///
/// `timelike_floor` must be finite and strictly positive; the observer is
/// required to be timelike with `g(u,u) < -timelike_floor` under a
/// `(-,+,+,+)`-signature convention. A non-timelike observer (null, spacelike,
/// or from an incompatible signature) is rejected with
/// [`NonlocalRelativityError::NonTimelikeMetricNorm`] rather than silently
/// producing a meaningless "spatial" length from an indefinite contraction.
///
/// This is a diagnostic, not a controller input, and it is chart-dependent in
/// the sense described in the module documentation; it must never be described
/// as a covariant error between two distant states or as establishing
/// coordinate covariance.
pub fn timelike_state_error<const D: usize>(
    metric: &[[f64; D]; D],
    four_velocity: &[f64; D],
    delta: &[f64; D],
    timelike_floor: f64,
) -> NonlocalResult<TimelikeStateError> {
    if !timelike_floor.is_finite() || timelike_floor <= 0.0
    {
        return Err(NonlocalRelativityError::InvalidMetricNormFloor(
            timelike_floor,
        ));
    }

    let metric_norm = metric_contraction(metric, four_velocity, four_velocity);
    if !metric_norm.is_finite()
    {
        return Err(NonlocalRelativityError::NonFiniteMetricNorm {
            step: 0,
            value: metric_norm,
        });
    }
    if metric_norm > -timelike_floor
    {
        // Not sufficiently timelike under the (-,+,+,+) convention.
        return Err(NonlocalRelativityError::NonTimelikeMetricNorm { metric_norm });
    }

    let lowered = lower_index(metric, four_velocity);
    let inner_delta_u = lowered
        .iter()
        .zip(delta)
        .fold(0.0, |sum, (lowered_component, delta_component)| {
            sum + *lowered_component * delta_component
        });
    validate_scalar("timelike_state_error_inner", inner_delta_u, 0)?;

    let temporal = inner_delta_u.abs() / (-metric_norm).sqrt();
    validate_scalar("timelike_state_error_temporal", temporal, 0)?;

    // Projection orthogonal to u: (P delta)^mu = delta^mu - u^mu * g(delta,u)/g(u,u).
    let projection_scale = inner_delta_u / metric_norm;
    let mut projected = [0.0_f64; D];
    for component in 0..D
    {
        projected[component] = delta[component] - four_velocity[component] * projection_scale;
        validate_scalar(
            "timelike_state_error_projected",
            projected[component],
            component,
        )?;
    }

    let spatial_squared = metric_contraction(metric, &projected, &projected);
    validate_scalar("timelike_state_error_spatial_squared", spatial_squared, 0)?;
    // The projected vector is spacelike, so this is non-negative up to
    // rounding; clamp a tiny negative rounding excursion to zero.
    let spatial = spatial_squared.max(0.0).sqrt();
    validate_scalar("timelike_state_error_spatial", spatial, 0)?;

    let orthogonality_residual =
        lowered
            .iter()
            .zip(projected)
            .fold(0.0, |sum, (lowered_component, projected_component)| {
                sum + *lowered_component * projected_component
            });
    validate_scalar(
        "timelike_state_error_orthogonality",
        orthogonality_residual,
        0,
    )?;

    Ok(TimelikeStateError {
        temporal,
        spatial,
        orthogonality_residual,
    })
}

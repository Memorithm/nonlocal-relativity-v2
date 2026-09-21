//! Reusable geometric comparison of two nearby candidate worldline states.
//!
//! This module implements the comparison proposed in research issue #42
//! without changing the adaptive acceptance controller. The low candidate is
//! compared at the high candidate point by combining a local logarithm map,
//! parallel transport of the low velocity, and observer-frame tetrad
//! projections. The result exposes each geometric component separately.
//!
//! This is a numerical diagnostic. Using geometric primitives does not by
//! itself establish coordinate invariance or covariance; cross-chart evidence
//! remains required before making either claim.

use crate::{
    NonlocalRelativityError, TetradStateError, WorldlineState, metric_contraction,
    tetrad_state_error, validate_scalar,
};
use scirust_relativity::{
    Connection, Metric, RelativityError, geodesic_logarithm, transport_along_segment,
};
use std::error::Error;
use std::fmt;

/// Numerical settings used to compare two nearby candidate states.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GeometricComparisonConfig {
    log_step: f64,
    jacobian_step: f64,
    log_tolerance: f64,
    log_max_iterations: usize,
    transport_substeps: usize,
    timelike_floor: f64,
}

impl GeometricComparisonConfig {
    /// Validate and construct geometric-comparison settings.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        log_step: f64,
        jacobian_step: f64,
        log_tolerance: f64,
        log_max_iterations: usize,
        transport_substeps: usize,
        timelike_floor: f64,
    ) -> Result<Self, GeometricComparisonError> {
        validate_positive_setting("log_step", log_step)?;
        validate_positive_setting("jacobian_step", jacobian_step)?;
        validate_positive_setting("log_tolerance", log_tolerance)?;
        validate_positive_setting("timelike_floor", timelike_floor)?;
        if log_max_iterations == 0
        {
            return Err(GeometricComparisonError::InvalidCount {
                parameter: "log_max_iterations",
                value: log_max_iterations,
            });
        }
        if transport_substeps == 0
        {
            return Err(GeometricComparisonError::InvalidCount {
                parameter: "transport_substeps",
                value: transport_substeps,
            });
        }

        Ok(Self {
            log_step,
            jacobian_step,
            log_tolerance,
            log_max_iterations,
            transport_substeps,
            timelike_floor,
        })
    }

    /// Geodesic RK4 step used by the logarithm map.
    #[must_use]
    pub const fn log_step(self) -> f64 {
        self.log_step
    }

    /// Central-difference step used by the logarithm-map Jacobian.
    #[must_use]
    pub const fn jacobian_step(self) -> f64 {
        self.jacobian_step
    }

    /// Euclidean endpoint tolerance used by logarithm-map Newton iteration.
    #[must_use]
    pub const fn log_tolerance(self) -> f64 {
        self.log_tolerance
    }

    /// Maximum logarithm-map Newton iterations.
    #[must_use]
    pub const fn log_max_iterations(self) -> usize {
        self.log_max_iterations
    }

    /// RK4 substeps used to transport the low candidate velocity.
    #[must_use]
    pub const fn transport_substeps(self) -> usize {
        self.transport_substeps
    }

    /// Timelike metric-norm floor used for observer-frame construction.
    #[must_use]
    pub const fn timelike_floor(self) -> f64 {
        self.timelike_floor
    }
}

/// Failure to construct a geometric comparison.
///
/// Failures are explicit: invalid numerical settings, an unavailable local
/// logarithm/transport operation, or a failed observer-frame/non-finite
/// diagnostic never falls back to a coordinate-component difference.
#[derive(Debug, Clone, PartialEq)]
pub enum GeometricComparisonError {
    /// A scalar numerical setting must be finite and strictly positive.
    InvalidPositiveSetting {
        /// Stable setting name.
        parameter: &'static str,
        /// Rejected value.
        value: f64,
    },
    /// An integer numerical setting must be non-zero.
    InvalidCount {
        /// Stable setting name.
        parameter: &'static str,
        /// Rejected value.
        value: usize,
    },
    /// A geometry-core logarithm-map or parallel-transport operation failed.
    Relativity {
        /// Stable operation name.
        operation: &'static str,
        /// Typed geometry-core failure.
        source: RelativityError,
    },
    /// Observer-frame or finite diagnostic validation failed.
    Nonlocal {
        /// Stable operation name.
        operation: &'static str,
        /// Typed nonlocal diagnostic failure.
        source: NonlocalRelativityError,
    },
}

impl fmt::Display for GeometricComparisonError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self
        {
            Self::InvalidPositiveSetting { parameter, value } => write!(
                formatter,
                "{parameter} must be finite and strictly positive; got {value}"
            ),
            Self::InvalidCount { parameter, value } =>
            {
                write!(formatter, "{parameter} must be non-zero; got {value}")
            },
            Self::Relativity { operation, source } =>
            {
                write!(formatter, "{operation} failed: {source}")
            },
            Self::Nonlocal { operation, source } =>
            {
                write!(formatter, "{operation} failed: {source}")
            },
        }
    }
}

impl Error for GeometricComparisonError {}

/// Geometric discrepancy between low- and high-order candidate states.
///
/// Position and velocity discrepancies live in the tangent space at the high
/// candidate. The low velocity is first parallel-transported to that point.
/// Frame fields are tetrad projections relative to the high candidate velocity.
/// Metric norms are exposed independently so transport drift remains visible.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GeometricTrialDiscrepancy<const D: usize> {
    /// Tangent log_(x_high)(x_low) at the high candidate.
    pub position_tangent: [f64; D],
    /// Low candidate velocity parallel-transported to x_high.
    pub transported_low_velocity: [f64; D],
    /// transported_low_velocity - high.velocity in the common tangent space.
    pub velocity_tangent: [f64; D],
    /// Tetrad components/magnitudes of the position tangent at x_high.
    pub position_frame: TetradStateError<D>,
    /// Tetrad components/magnitudes of the velocity tangent at x_high.
    pub velocity_frame: TetradStateError<D>,
    /// g(u_low, u_low) at x_low.
    pub low_velocity_metric_norm: f64,
    /// g(u_high, u_high) at x_high.
    pub high_velocity_metric_norm: f64,
    /// g(PT(u_low), PT(u_low)) at x_high.
    pub transported_low_velocity_metric_norm: f64,
    /// Absolute metric-norm drift introduced by numerical transport.
    pub transport_metric_norm_drift: f64,
}

/// Compare two nearby candidate worldline states in a common tangent space.
///
/// The logarithm is based at the high candidate and targets the low candidate.
/// The low velocity is transported along the geometry-core straight coordinate
/// segment from low to high before subtraction. Both tangent discrepancies are
/// then projected into a local orthonormal frame defined by high.velocity.
///
/// This routine deliberately returns no accept/reject decision. It is the
/// reusable diagnostic substrate needed before issue #42 can compare controller
/// decisions across charts.
///
/// # Errors
///
/// Returns GeometricComparisonError and fails closed when the logarithm map,
/// transport, tetrad construction, or finite metric diagnostics cannot be
/// constructed with the supplied settings.
pub fn geometric_trial_discrepancy<B, const D: usize>(
    background: &B,
    low: &WorldlineState<D>,
    high: &WorldlineState<D>,
    config: GeometricComparisonConfig,
) -> Result<GeometricTrialDiscrepancy<D>, GeometricComparisonError>
where
    B: Connection<D> + Metric<D> + Copy,
{
    let position_tangent = geodesic_logarithm(
        background,
        &high.coordinates,
        &low.coordinates,
        config.log_step,
        config.jacobian_step,
        config.log_tolerance,
        config.log_max_iterations,
    )
    .map_err(|source| GeometricComparisonError::Relativity {
        operation: "geodesic_logarithm",
        source,
    })?;

    let transported_low_velocity = transport_along_segment(
        background,
        &low.coordinates,
        &high.coordinates,
        &low.velocity,
        config.transport_substeps,
    )
    .map_err(|source| GeometricComparisonError::Relativity {
        operation: "transport_along_segment",
        source,
    })?;

    let mut velocity_tangent = [0.0_f64; D];
    for (component, slot) in velocity_tangent.iter_mut().enumerate()
    {
        *slot = transported_low_velocity[component] - high.velocity[component];
        validate_scalar("geometric_velocity_tangent", *slot, component).map_err(|source| {
            GeometricComparisonError::Nonlocal {
                operation: "velocity_tangent",
                source,
            }
        })?;
    }

    let high_metric = background.components(&high.coordinates);
    let position_frame = tetrad_state_error(
        &high_metric,
        &high.velocity,
        &position_tangent,
        config.timelike_floor,
    )
    .map_err(|source| GeometricComparisonError::Nonlocal {
        operation: "position_tetrad",
        source,
    })?;
    let velocity_frame = tetrad_state_error(
        &high_metric,
        &high.velocity,
        &velocity_tangent,
        config.timelike_floor,
    )
    .map_err(|source| GeometricComparisonError::Nonlocal {
        operation: "velocity_tetrad",
        source,
    })?;

    let low_metric = background.components(&low.coordinates);
    let low_velocity_metric_norm =
        metric_contraction(&low_metric, &low.velocity, &low.velocity);
    let high_velocity_metric_norm =
        metric_contraction(&high_metric, &high.velocity, &high.velocity);
    let transported_low_velocity_metric_norm = metric_contraction(
        &high_metric,
        &transported_low_velocity,
        &transported_low_velocity,
    );
    let transport_metric_norm_drift =
        (transported_low_velocity_metric_norm - low_velocity_metric_norm).abs();

    for (quantity, value) in [
        ("geometric_low_velocity_metric_norm", low_velocity_metric_norm),
        ("geometric_high_velocity_metric_norm", high_velocity_metric_norm),
        (
            "geometric_transported_low_velocity_metric_norm",
            transported_low_velocity_metric_norm,
        ),
        (
            "geometric_transport_metric_norm_drift",
            transport_metric_norm_drift,
        ),
    ]
    {
        validate_scalar(quantity, value, 0).map_err(|source| {
            GeometricComparisonError::Nonlocal {
                operation: quantity,
                source,
            }
        })?;
    }

    Ok(GeometricTrialDiscrepancy {
        position_tangent,
        transported_low_velocity,
        velocity_tangent,
        position_frame,
        velocity_frame,
        low_velocity_metric_norm,
        high_velocity_metric_norm,
        transported_low_velocity_metric_norm,
        transport_metric_norm_drift,
    })
}

fn validate_positive_setting(
    parameter: &'static str,
    value: f64,
) -> Result<(), GeometricComparisonError> {
    if value.is_finite() && value > 0.0
    {
        Ok(())
    }
    else
    {
        Err(GeometricComparisonError::InvalidPositiveSetting { parameter, value })
    }
}

#[cfg(test)]
mod tests {
    use super::{
        GeometricComparisonConfig, GeometricComparisonError, geometric_trial_discrepancy,
    };
    use crate::WorldlineState;
    use scirust_relativity::{Minkowski, RelativityError};

    fn config() -> GeometricComparisonConfig {
        GeometricComparisonConfig::try_new(0.05, 1.0e-6, 1.0e-12, 16, 8, 1.0e-12).unwrap()
    }

    #[test]
    fn flat_comparison_uses_common_high_tangent_space() {
        let low = WorldlineState::new(
            [0.1, 1.0, -0.5, 0.2],
            [1.0, 0.10, -0.02, 0.03],
        );
        let high = WorldlineState::new(
            [0.2, 1.2, -0.4, 0.25],
            [1.0, 0.12, -0.01, 0.02],
        );

        let comparison =
            geometric_trial_discrepancy(&Minkowski, &low, &high, config()).unwrap();

        let expected_position = [-0.1, -0.2, -0.1, -0.05];
        let expected_velocity = [0.0, -0.02, -0.01, 0.01];
        for component in 0..4
        {
            assert!((comparison.position_tangent[component] - expected_position[component]).abs() < 2.0e-12);
            assert_eq!(
                comparison.transported_low_velocity[component].to_bits(),
                low.velocity[component].to_bits()
            );
            assert!((comparison.velocity_tangent[component] - expected_velocity[component]).abs() < 2.0e-15);
        }
        assert!(comparison.transport_metric_norm_drift < 2.0e-15);
    }

    #[test]
    fn invalid_configuration_fails_before_geometry() {
        assert!(matches!(
            GeometricComparisonConfig::try_new(0.0, 1.0e-6, 1.0e-12, 16, 8, 1.0e-12),
            Err(GeometricComparisonError::InvalidPositiveSetting {
                parameter: "log_step",
                ..
            })
        ));
        assert!(matches!(
            GeometricComparisonConfig::try_new(0.05, 1.0e-6, 1.0e-12, 0, 8, 1.0e-12),
            Err(GeometricComparisonError::InvalidCount {
                parameter: "log_max_iterations",
                value: 0
            })
        ));
    }

    #[test]
    fn unavailable_logarithm_fails_closed() {
        let low = WorldlineState::new(
            [f64::NAN, 0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0, 0.0],
        );
        let high = WorldlineState::new(
            [0.0, 0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0, 0.0],
        );

        assert!(matches!(
            geometric_trial_discrepancy(&Minkowski, &low, &high, config()),
            Err(GeometricComparisonError::Relativity {
                operation: "geodesic_logarithm",
                source: RelativityError::NonFiniteCoordinate(_)
            })
        ));
    }
}

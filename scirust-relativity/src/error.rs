//! Errors reported by differential-geometry operations.

use std::error::Error;
use std::fmt;

/// Errors produced while evaluating metrics and connections.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum RelativityError {
    /// A coordinate component is not finite.
    NonFiniteCoordinate(usize),

    /// A metric component is not finite.
    NonFiniteMetricComponent {
        /// Row index.
        row: usize,
        /// Column index.
        column: usize,
    },

    /// The metric matrix is singular or numerically non-invertible.
    SingularMetric,

    /// The finite-difference step is non-finite or non-positive.
    InvalidDifferenceStep(f64),

    /// A curvature-tensor component evaluated to a non-finite value.
    NonFiniteCurvatureComponent {
        /// Short name of the offending quantity (for example `"riemann"`,
        /// `"ricci"`, `"ricci_scalar"`, `"einstein"`, `"kretschmann"`).
        quantity: &'static str,
    },
}

impl fmt::Display for RelativityError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self
        {
            Self::NonFiniteCoordinate(index) =>
            {
                write!(formatter, "coordinate at index {index} is not finite")
            },
            Self::NonFiniteMetricComponent { row, column } =>
            {
                write!(
                    formatter,
                    "metric component ({row}, {column}) is not finite"
                )
            },
            Self::SingularMetric =>
            {
                write!(
                    formatter,
                    "metric is singular or numerically non-invertible"
                )
            },
            Self::InvalidDifferenceStep(step) =>
            {
                write!(
                    formatter,
                    "finite-difference step must be finite and positive; got {step}"
                )
            },
            Self::NonFiniteCurvatureComponent { quantity } =>
            {
                write!(formatter, "curvature quantity '{quantity}' is not finite")
            },
        }
    }
}

impl Error for RelativityError {}

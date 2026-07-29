//! Uniform periodic one-dimensional grids and centred finite differences.
//!
//! This module is deliberately free of any relativity: it is the discretisation
//! substrate that [`crate::bssn_grid`] builds on, and nothing here knows what a
//! metric is.
//!
//! # Convention
//!
//! The domain is **half-open**. For `N` points on `[x_min, x_max)`:
//!
//! ```text
//! x_n = x_min + n * dx,    n = 0, ..., N-1,    dx = (x_max - x_min) / N
//! ```
//!
//! The upper endpoint is **not** stored — it is the periodic image of the
//! lower one. Storing both would duplicate a degree of freedom.
//!
//! # Periodic indexing
//!
//! Neighbour lookup wraps with Euclidean remainder, never with unsigned integer
//! wraparound, so the left neighbour of index `0` is `N-1`, the right neighbour
//! of `N-1` is `0`, and no offset of any sign can underflow.
//!
//! # Finite differences
//!
//! Second-order centred:
//!
//! ```text
//! D1 f_i = ( f_{i+1} - f_{i-1} ) / ( 2 dx )
//! D2 f_i = ( f_{i+1} - 2 f_i + f_{i-1} ) / ( dx^2 )
//! ```
//!
//! Both are validated against `sin(kx)` at several wave numbers and several
//! resolutions; the observed order is measured, not assumed from the formula.
//!
//! ```
//! use scirust_relativity::grid1d::{UniformGrid1d, periodic_first_derivative};
//!
//! let grid = UniformGrid1d::new(64, 0.0, 1.0).unwrap();
//! let k = 2.0 * std::f64::consts::PI;
//! let samples: Vec<f64> = (0..grid.points()).map(|i| (k * grid.coordinate(i)).sin()).collect();
//! let numerical = periodic_first_derivative(&samples, &grid, 8).unwrap();
//! let exact = k * (k * grid.coordinate(8)).cos();
//! assert!((numerical - exact).abs() < 1.0e-2);
//! ```

use core::fmt;

/// The minimum number of points a second-order centred stencil needs.
///
/// Three points is the smallest grid on which `f_{i-1}`, `f_i`, and `f_{i+1}`
/// are not all the same sample under periodic wrapping.
pub const MINIMUM_STENCIL_POINTS: usize = 3;

/// A failure constructing or using a grid.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Grid1dError {
    /// Fewer points than the centred stencil requires.
    InsufficientPoints {
        /// The requested point count.
        requested: usize,
        /// The minimum the stencil needs.
        minimum: usize,
    },
    /// A domain bound is not finite.
    NonFiniteBound {
        /// Which bound: `"lower"` or `"upper"`.
        bound: &'static str,
    },
    /// The domain length is not strictly positive.
    NonPositiveLength,
    /// The derived spacing is not finite and strictly positive.
    InvalidSpacing,
    /// A sample array's length does not match the grid.
    LengthMismatch {
        /// The grid's point count.
        expected: usize,
        /// The supplied array's length.
        actual: usize,
    },
}

impl fmt::Display for Grid1dError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self
        {
            Self::InsufficientPoints { requested, minimum } => write!(
                f,
                "grid needs at least {minimum} points for a centred stencil, got {requested}"
            ),
            Self::NonFiniteBound { bound } => write!(f, "{bound} domain bound is not finite"),
            Self::NonPositiveLength => write!(f, "domain length must be strictly positive"),
            Self::InvalidSpacing => write!(f, "grid spacing must be finite and strictly positive"),
            Self::LengthMismatch { expected, actual } =>
            {
                write!(f, "expected {expected} samples for this grid, got {actual}")
            },
        }
    }
}

impl std::error::Error for Grid1dError {}

/// A uniform, half-open, periodic one-dimensional grid.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct UniformGrid1d {
    points: usize,
    lower: f64,
    spacing: f64,
}

impl UniformGrid1d {
    /// Build a grid of `points` samples spanning the half-open domain
    /// `[lower, upper)`.
    ///
    /// The upper bound is *not* a stored point: it is the periodic image of
    /// `lower`, so the spacing is `(upper - lower) / points`, not
    /// `(upper - lower) / (points - 1)`.
    pub fn new(points: usize, lower: f64, upper: f64) -> Result<Self, Grid1dError> {
        if points < MINIMUM_STENCIL_POINTS
        {
            return Err(Grid1dError::InsufficientPoints {
                requested: points,
                minimum: MINIMUM_STENCIL_POINTS,
            });
        }
        if !lower.is_finite()
        {
            return Err(Grid1dError::NonFiniteBound { bound: "lower" });
        }
        if !upper.is_finite()
        {
            return Err(Grid1dError::NonFiniteBound { bound: "upper" });
        }
        // Both bounds are finite by the checks above, so the difference is
        // finite and a plain comparison is exact here.
        let length = upper - lower;
        if length <= 0.0
        {
            return Err(Grid1dError::NonPositiveLength);
        }
        let spacing = length / points as f64;
        if !spacing.is_finite() || spacing <= 0.0
        {
            return Err(Grid1dError::InvalidSpacing);
        }
        Ok(Self {
            points,
            lower,
            spacing,
        })
    }

    /// The number of stored points `N`.
    #[must_use]
    pub const fn points(&self) -> usize {
        self.points
    }

    /// The spacing `dx`.
    #[must_use]
    pub const fn spacing(&self) -> f64 {
        self.spacing
    }

    /// The lower domain bound `x_min` (a stored point).
    #[must_use]
    pub const fn lower(&self) -> f64 {
        self.lower
    }

    /// The upper domain bound `x_max` (**not** a stored point).
    #[must_use]
    pub fn upper(&self) -> f64 {
        self.lower + self.length()
    }

    /// The periodic domain length `x_max - x_min`.
    #[must_use]
    pub fn length(&self) -> f64 {
        self.spacing * self.points as f64
    }

    /// The coordinate of stored point `index`, which is wrapped periodically.
    #[must_use]
    pub fn coordinate(&self, index: usize) -> f64 {
        self.lower + self.wrap_usize(index) as f64 * self.spacing
    }

    /// Wrap a `usize` index into `[0, N)`.
    #[must_use]
    pub const fn wrap_usize(&self, index: usize) -> usize {
        index % self.points
    }

    /// The index `offset` steps from `index`, wrapped periodically.
    ///
    /// `offset` may be negative and of any magnitude; the Euclidean remainder
    /// makes underflow impossible.
    #[must_use]
    pub fn offset(&self, index: usize, offset: isize) -> usize {
        let n = self.points as isize;
        let raw = (index % self.points) as isize + offset;
        raw.rem_euclid(n) as usize
    }

    /// The stored index whose coordinate is nearest `x`, wrapped periodically.
    ///
    /// This is how a grid-backed field answers a query at an arbitrary
    /// coordinate. Rounding — rather than truncation — is what makes a sample
    /// taken at exactly `x_i + dx` land on index `i+1` despite floating-point
    /// representation error.
    #[must_use]
    pub fn nearest_index(&self, x: f64) -> usize {
        if !x.is_finite()
        {
            return 0;
        }
        let raw = ((x - self.lower) / self.spacing).round();
        // `raw` is finite here because `x`, `lower`, and `spacing` all are and
        // `spacing > 0`; the modulus keeps the cast in range.
        let n = self.points as f64;
        let wrapped = raw - (raw / n).floor() * n;
        let index = wrapped as usize;
        // Guard the boundary case where rounding lands exactly on `n`.
        self.wrap_usize(index)
    }
}

fn check_length(samples: &[f64], grid: &UniformGrid1d) -> Result<(), Grid1dError> {
    if samples.len() != grid.points()
    {
        return Err(Grid1dError::LengthMismatch {
            expected: grid.points(),
            actual: samples.len(),
        });
    }
    Ok(())
}

/// The second-order centred periodic first derivative at `index`.
///
/// ```text
/// D1 f_i = ( f_{i+1} - f_{i-1} ) / ( 2 dx )
/// ```
pub fn periodic_first_derivative(
    samples: &[f64],
    grid: &UniformGrid1d,
    index: usize,
) -> Result<f64, Grid1dError> {
    check_length(samples, grid)?;
    let left = samples[grid.offset(index, -1)];
    let right = samples[grid.offset(index, 1)];
    Ok((right - left) / (2.0 * grid.spacing()))
}

/// The second-order centred periodic second derivative at `index`.
///
/// ```text
/// D2 f_i = ( f_{i+1} - 2 f_i + f_{i-1} ) / ( dx^2 )
/// ```
pub fn periodic_second_derivative(
    samples: &[f64],
    grid: &UniformGrid1d,
    index: usize,
) -> Result<f64, Grid1dError> {
    check_length(samples, grid)?;
    let centre = samples[grid.wrap_usize(index)];
    let left = samples[grid.offset(index, -1)];
    let right = samples[grid.offset(index, 1)];
    Ok((right - 2.0 * centre + left) / (grid.spacing() * grid.spacing()))
}

/// Fill `out` with the centred periodic first derivative of `samples`.
///
/// Deterministic ascending point order.
pub fn periodic_first_derivative_all(
    samples: &[f64],
    grid: &UniformGrid1d,
    out: &mut [f64],
) -> Result<(), Grid1dError> {
    check_length(samples, grid)?;
    check_length(out, grid)?;
    let scale = 1.0 / (2.0 * grid.spacing());
    for index in 0..grid.points()
    {
        out[index] = (samples[grid.offset(index, 1)] - samples[grid.offset(index, -1)]) * scale;
    }
    Ok(())
}

/// Fill `out` with the centred periodic second derivative of `samples`.
///
/// Deterministic ascending point order.
pub fn periodic_second_derivative_all(
    samples: &[f64],
    grid: &UniformGrid1d,
    out: &mut [f64],
) -> Result<(), Grid1dError> {
    check_length(samples, grid)?;
    check_length(out, grid)?;
    let scale = 1.0 / (grid.spacing() * grid.spacing());
    for index in 0..grid.points()
    {
        out[index] = (samples[grid.offset(index, 1)] - 2.0 * samples[index]
            + samples[grid.offset(index, -1)])
            * scale;
    }
    Ok(())
}

/// A deterministic reduction of a residual field over the grid.
///
/// Every component is reported; nothing is blended into a single score.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GridReduction {
    /// The signed arithmetic mean.
    pub signed_mean: f64,
    /// The discrete `L1` norm `(1/N) sum |f_i|`.
    pub l1: f64,
    /// The discrete `L2` norm `sqrt( (1/N) sum f_i^2 )`.
    pub l2: f64,
    /// The maximum absolute value.
    pub max_abs: f64,
    /// The lowest index attaining `max_abs` (ties resolved deterministically).
    pub max_index: usize,
    /// How many samples were not finite.
    pub non_finite: usize,
}

impl GridReduction {
    /// Reduce `samples` in ascending index order.
    ///
    /// Non-finite samples are counted rather than propagated into the norms, so
    /// a single `NaN` reports itself instead of erasing every other diagnostic.
    #[must_use]
    pub fn of(samples: &[f64]) -> Self {
        let mut sum = 0.0_f64;
        let mut absolute = 0.0_f64;
        let mut squared = 0.0_f64;
        let mut max_abs = 0.0_f64;
        let mut max_index = 0_usize;
        let mut non_finite = 0_usize;
        let mut counted = 0_usize;

        for (index, &value) in samples.iter().enumerate()
        {
            if !value.is_finite()
            {
                non_finite += 1;
                continue;
            }
            counted += 1;
            sum += value;
            absolute += value.abs();
            squared += value * value;
            if value.abs() > max_abs
            {
                max_abs = value.abs();
                max_index = index;
            }
        }

        let denominator = if counted == 0 { 1.0 } else { counted as f64 };
        Self {
            signed_mean: sum / denominator,
            l1: absolute / denominator,
            l2: (squared / denominator).sqrt(),
            max_abs,
            max_index,
            non_finite,
        }
    }
}

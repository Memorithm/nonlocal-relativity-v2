//! Synge's world function and its first-derivative bitensors.
//!
//! Synge's world function `sigma(x', x)` is the biscalar equal to one-half the
//! signed squared geodesic distance between a base point `x'` and a field point
//! `x`: it is negative for timelike separation, positive for spacelike, and zero
//! for null (under the `(-,+,+,+)` convention). With the unique geodesic linking
//! the two points parametrized affinely on `[0, 1]`,
//!
//! ```text
//! sigma(x', x) = (1/2) g_{mu nu}(x') v^mu v^nu,   v = log_{x'}(x),
//! ```
//!
//! where `v` is the tangent at `x'` of that geodesic — exactly the geodesic
//! logarithm map ([`crate::geodesic_logarithm`]). Because `g(v, v)` is conserved
//! along an affine geodesic, evaluating it at `x'` gives the whole world
//! function; no separate distance integral is needed.
//!
//! Its first covariant derivatives are bitensors (each transforms as a vector at
//! one of the two points):
//!
//! - at the base point, `sigma^{mu'} = -v` (minus the tangent at `x'`);
//! - at the field point, `sigma^{mu} = -log_x(x')` (the tangent at `x` of the
//!   geodesic from `x'` to `x`, directed away from `x'`).
//!
//! Both obey the fundamental identity `2 sigma = g_{mu nu} sigma^mu sigma^nu`
//! (contracted with the metric at the respective point). The base-point identity
//! holds by construction; the field-point identity, evaluated with an
//! independent logarithm shooting from `x`, is a genuine convention-free check
//! that the two ends agree on the same geodesic distance — as is the symmetry
//! `sigma(x', x) = sigma(x, x')`.
//!
//! In flat spacetime everything is exact: `sigma = (1/2) eta_{mu nu} Δ^mu Δ^nu`
//! with `Δ = x - x'`, `sigma^{mu'} = -Δ`, and `sigma^mu = Δ`. This is established
//! general relativity, built entirely on the existing geodesic exponential /
//! logarithm maps; no new geodesic solver is introduced.

use crate::{Connection, Metric, RelativityError, geodesic_logarithm, metric_norm};

/// Numerical controls for the geodesic shooting used to evaluate the world
/// function, forwarded to [`crate::geodesic_logarithm`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WorldFunctionSettings {
    /// RK4 affine step for the underlying geodesic integration.
    pub step: f64,
    /// Central-difference step for the exponential map's Jacobian in the Newton
    /// shooting.
    pub jacobian_step: f64,
    /// Euclidean convergence tolerance on the shooting residual.
    pub tolerance: f64,
    /// Maximum number of Newton iterations.
    pub max_iterations: usize,
}

impl Default for WorldFunctionSettings {
    /// The settings used by the crate's world-function tests: `step = 1e-2`,
    /// `jacobian_step = 1e-5`, `tolerance = 1e-10`, `max_iterations = 50`.
    fn default() -> Self {
        Self {
            step: 1.0e-2,
            jacobian_step: 1.0e-5,
            tolerance: 1.0e-10,
            max_iterations: 50,
        }
    }
}

/// Synge's world function `sigma(base, field) = (1/2) g_{mu nu}(base) v^mu v^nu`,
/// with `v = log_base(field)` the tangent at `base` of the geodesic reaching
/// `field`.
///
/// Negative for timelike separation, positive for spacelike, zero for null
/// (under the `(-,+,+,+)` convention). Returns a typed [`RelativityError`] if the
/// connecting geodesic cannot be found (propagated from
/// [`crate::geodesic_logarithm`]) or the result is non-finite; it never panics.
///
/// # Example
///
/// In flat spacetime the world function is one-half the squared coordinate
/// separation. For the purely spatial separation `Δ = (0, 3, 4, 0)`,
/// `sigma = (1/2)(3^2 + 4^2) = 12.5`.
///
/// ```
/// use scirust_relativity::{Minkowski, WorldFunctionSettings, world_function};
///
/// let base = [0.0, 0.0, 0.0, 0.0];
/// let field = [0.0, 3.0, 4.0, 0.0];
/// let sigma = world_function(&Minkowski, &base, &field, &WorldFunctionSettings::default())
///     .expect("flat world function");
/// assert!((sigma - 12.5).abs() < 1.0e-12);
/// ```
pub fn world_function<B, const D: usize>(
    background: &B,
    base: &[f64; D],
    field: &[f64; D],
    settings: &WorldFunctionSettings,
) -> Result<f64, RelativityError>
where
    B: Metric<D> + Connection<D> + Copy,
{
    let tangent = geodesic_logarithm(
        background,
        base,
        field,
        settings.step,
        settings.jacobian_step,
        settings.tolerance,
        settings.max_iterations,
    )?;
    let metric_base = background.components(base);
    let sigma = 0.5 * metric_norm(&metric_base, &tangent);
    if !sigma.is_finite()
    {
        return Err(RelativityError::NonFiniteWorldFunction);
    }
    Ok(sigma)
}

/// Synge's world function together with its first-derivative gradient bitensors
/// at both the base and field points.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WorldFunction<const D: usize> {
    /// `sigma(base, field)`; see [`world_function`].
    pub sigma: f64,
    /// The contravariant gradient at the base point,
    /// `sigma^{mu'} = -[log_base(field)]^mu` — minus the tangent at `base` of the
    /// connecting geodesic. Satisfies `g_{mu nu}(base) sigma^{mu'} sigma^{nu'} =
    /// 2 sigma` by construction.
    pub gradient_base: [f64; D],
    /// The contravariant gradient at the field point,
    /// `sigma^mu = -[log_field(base)]^mu` — the tangent at `field` of the
    /// connecting geodesic, directed away from `base`. Satisfies
    /// `g_{mu nu}(field) sigma^mu sigma^nu = 2 sigma` (a convention-free check,
    /// since it uses an independent shooting from `field`).
    pub gradient_field: [f64; D],
}

/// Evaluate Synge's world function and both gradient bitensors for the base
/// point `base` and field point `field`.
///
/// This performs two geodesic logarithm shootings — one from each endpoint — so
/// the field-point gradient (and hence the fundamental identity `2 sigma =
/// g(field) sigma^mu sigma^mu`) is an independent cross-check on the world
/// function computed from the base point. Returns a typed [`RelativityError`] on
/// a shooting failure or a non-finite result; it never panics.
pub fn world_function_with_gradients<B, const D: usize>(
    background: &B,
    base: &[f64; D],
    field: &[f64; D],
    settings: &WorldFunctionSettings,
) -> Result<WorldFunction<D>, RelativityError>
where
    B: Metric<D> + Connection<D> + Copy,
{
    // Tangent at the base point of the geodesic reaching the field point.
    let tangent_base = geodesic_logarithm(
        background,
        base,
        field,
        settings.step,
        settings.jacobian_step,
        settings.tolerance,
        settings.max_iterations,
    )?;
    let metric_base = background.components(base);
    let sigma = 0.5 * metric_norm(&metric_base, &tangent_base);

    // Tangent at the field point of the geodesic reaching the base point.
    let tangent_field = geodesic_logarithm(
        background,
        field,
        base,
        settings.step,
        settings.jacobian_step,
        settings.tolerance,
        settings.max_iterations,
    )?;

    let gradient_base = tangent_base.map(|component| -component);
    let gradient_field = tangent_field.map(|component| -component);

    if !sigma.is_finite()
        || gradient_base
            .iter()
            .chain(gradient_field.iter())
            .any(|component| !component.is_finite())
    {
        return Err(RelativityError::NonFiniteWorldFunction);
    }

    Ok(WorldFunction {
        sigma,
        gradient_base,
        gradient_field,
    })
}

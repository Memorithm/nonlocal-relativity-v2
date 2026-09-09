//! Research experiment: geometric local-state discrepancy across two charts.
//!
//! This is the first executable slice of issue #42. It deliberately does not
//! replace the adaptive controller. Instead it composes geometry primitives
//! that already exist in SciRust and checks them on the same physical pair of
//! nearby states represented in Cartesian and cylindrical Minkowski charts:
//!
//! 1. `log_{x_high}(x_low)` for the position discrepancy;
//! 2. parallel transport of `u_low` to `x_high`;
//! 3. velocity discrepancy in the common tangent space at `x_high`;
//! 4. local observer-frame (tetrad) decomposition of both discrepancies;
//! 5. metric-norm comparisons as an independent structural diagnostic.
//!
//! The cylindrical transport uses the generic numerical segment transporter
//! and is cross-checked against the exact flat Cartesian/cylindrical transport
//! oracle. Passing this experiment is evidence that the proposed diagnostic is
//! consistent for this flat, local test case. It is not a proof of coordinate
//! invariance or covariance and it is not yet an accept/reject controller.

#![forbid(unsafe_code)]

use nonlocal_relativity_experiments::{print_common_header, require_finite};
use scirust_nonlocal_relativity::{
    CylindricalMinkowski, cartesian_to_cylindrical_coordinates,
    cartesian_to_cylindrical_velocity, cylindrical_to_cartesian_velocity,
    exact_cylindrical_minkowski_transport, tetrad_state_error,
};
use scirust_relativity::{Metric, Minkowski, geodesic_logarithm, transport_along_segment};

const LOG_STEP: f64 = 0.05;
const JACOBIAN_STEP: f64 = 1.0e-6;
const LOG_TOLERANCE: f64 = 1.0e-11;
const LOG_MAX_ITERATIONS: usize = 16;
const TRANSPORT_SUBSTEPS: usize = 64;
const TIMELIKE_FLOOR: f64 = 1.0e-12;

const POSITION_VECTOR_TOLERANCE: f64 = 2.0e-8;
const VELOCITY_VECTOR_TOLERANCE: f64 = 2.0e-9;
const FRAME_MAGNITUDE_TOLERANCE: f64 = 2.0e-8;
const TRANSPORT_ORACLE_TOLERANCE: f64 = 2.0e-10;
const METRIC_NORM_TOLERANCE: f64 = 2.0e-13;

fn main() -> Result<(), String> {
    // Two nearby timelike states, expressed first in Cartesian Minkowski
    // coordinates. They are intentionally off-axis so the cylindrical
    // connection and Jacobians are genuinely exercised.
    let low_coordinates = [0.200, 5.000, 1.000, 0.300];
    let high_coordinates = [0.205, 5.004, 1.003, 0.298];
    let low_velocity = [1.030, 0.120, -0.040, 0.020];
    let high_velocity = [1.031, 0.119, -0.039, 0.021];

    let cartesian = Minkowski;
    let cartesian_metric = cartesian.components(&high_coordinates);
    let cartesian_position_tangent = geodesic_logarithm(
        &cartesian,
        &high_coordinates,
        &low_coordinates,
        LOG_STEP,
        JACOBIAN_STEP,
        LOG_TOLERANCE,
        LOG_MAX_ITERATIONS,
    )
    .map_err(stringify)?;
    let cartesian_low_velocity_at_high = transport_along_segment(
        &cartesian,
        &low_coordinates,
        &high_coordinates,
        &low_velocity,
        TRANSPORT_SUBSTEPS,
    )
    .map_err(stringify)?;
    let cartesian_velocity_difference = subtract(cartesian_low_velocity_at_high, high_velocity);
    let cartesian_position_frame = tetrad_state_error(
        &cartesian_metric,
        &high_velocity,
        &cartesian_position_tangent,
        TIMELIKE_FLOOR,
    )
    .map_err(stringify)?;
    let cartesian_velocity_frame = tetrad_state_error(
        &cartesian_metric,
        &high_velocity,
        &cartesian_velocity_difference,
        TIMELIKE_FLOOR,
    )
    .map_err(stringify)?;

    let cylindrical = CylindricalMinkowski;
    let low_cylindrical =
        cartesian_to_cylindrical_coordinates(low_coordinates).map_err(stringify)?;
    let high_cylindrical =
        cartesian_to_cylindrical_coordinates(high_coordinates).map_err(stringify)?;
    let low_velocity_cylindrical =
        cartesian_to_cylindrical_velocity(low_coordinates, low_velocity).map_err(stringify)?;
    let high_velocity_cylindrical =
        cartesian_to_cylindrical_velocity(high_coordinates, high_velocity).map_err(stringify)?;
    let cylindrical_metric = cylindrical.components(&high_cylindrical);

    let cylindrical_position_tangent = geodesic_logarithm(
        &cylindrical,
        &high_cylindrical,
        &low_cylindrical,
        LOG_STEP,
        JACOBIAN_STEP,
        LOG_TOLERANCE,
        LOG_MAX_ITERATIONS,
    )
    .map_err(stringify)?;
    let cylindrical_low_velocity_at_high = transport_along_segment(
        &cylindrical,
        &low_cylindrical,
        &high_cylindrical,
        &low_velocity_cylindrical,
        TRANSPORT_SUBSTEPS,
    )
    .map_err(stringify)?;
    let exact_cylindrical_low_velocity_at_high = exact_cylindrical_minkowski_transport(
        low_cylindrical,
        low_velocity_cylindrical,
        high_cylindrical,
    )
    .map_err(stringify)?;
    let cylindrical_transport_error = max_abs_difference(
        &cylindrical_low_velocity_at_high,
        &exact_cylindrical_low_velocity_at_high,
    );

    let cylindrical_velocity_difference =
        subtract(cylindrical_low_velocity_at_high, high_velocity_cylindrical);
    let cylindrical_position_frame = tetrad_state_error(
        &cylindrical_metric,
        &high_velocity_cylindrical,
        &cylindrical_position_tangent,
        TIMELIKE_FLOOR,
    )
    .map_err(stringify)?;
    let cylindrical_velocity_frame = tetrad_state_error(
        &cylindrical_metric,
        &high_velocity_cylindrical,
        &cylindrical_velocity_difference,
        TIMELIKE_FLOOR,
    )
    .map_err(stringify)?;

    // Bring the cylindrical tangent-space discrepancies back to Cartesian
    // components at the common high point. This is an independent chart-level
    // comparison in addition to the tetrad scalar magnitudes below.
    let cylindrical_position_tangent_cartesian = cylindrical_to_cartesian_velocity(
        high_cylindrical,
        cylindrical_position_tangent,
    )
    .map_err(stringify)?;
    let cylindrical_velocity_difference_cartesian = cylindrical_to_cartesian_velocity(
        high_cylindrical,
        cylindrical_velocity_difference,
    )
    .map_err(stringify)?;

    let position_vector_disagreement = max_abs_difference(
        &cartesian_position_tangent,
        &cylindrical_position_tangent_cartesian,
    );
    let velocity_vector_disagreement = max_abs_difference(
        &cartesian_velocity_difference,
        &cylindrical_velocity_difference_cartesian,
    );
    let position_temporal_disagreement =
        (cartesian_position_frame.temporal - cylindrical_position_frame.temporal).abs();
    let position_spatial_disagreement =
        (cartesian_position_frame.spatial - cylindrical_position_frame.spatial).abs();
    let velocity_temporal_disagreement =
        (cartesian_velocity_frame.temporal - cylindrical_velocity_frame.temporal).abs();
    let velocity_spatial_disagreement =
        (cartesian_velocity_frame.spatial - cylindrical_velocity_frame.spatial).abs();

    let low_norm_cartesian = metric_contraction(
        &cartesian.components(&low_coordinates),
        &low_velocity,
        &low_velocity,
    );
    let low_norm_cylindrical = metric_contraction(
        &cylindrical.components(&low_cylindrical),
        &low_velocity_cylindrical,
        &low_velocity_cylindrical,
    );
    let high_norm_cartesian = metric_contraction(
        &cartesian_metric,
        &high_velocity,
        &high_velocity,
    );
    let high_norm_cylindrical = metric_contraction(
        &cylindrical_metric,
        &high_velocity_cylindrical,
        &high_velocity_cylindrical,
    );
    let transported_low_norm = metric_contraction(
        &cylindrical_metric,
        &cylindrical_low_velocity_at_high,
        &cylindrical_low_velocity_at_high,
    );
    let metric_norm_chart_disagreement = (low_norm_cartesian - low_norm_cylindrical)
        .abs()
        .max((high_norm_cartesian - high_norm_cylindrical).abs());
    let transport_metric_norm_drift = (transported_low_norm - low_norm_cylindrical).abs();

    require_finite(&[
        ("cylindrical_transport_error", cylindrical_transport_error),
        ("position_vector_disagreement", position_vector_disagreement),
        ("velocity_vector_disagreement", velocity_vector_disagreement),
        (
            "position_temporal_disagreement",
            position_temporal_disagreement,
        ),
        ("position_spatial_disagreement", position_spatial_disagreement),
        (
            "velocity_temporal_disagreement",
            velocity_temporal_disagreement,
        ),
        ("velocity_spatial_disagreement", velocity_spatial_disagreement),
        ("metric_norm_chart_disagreement", metric_norm_chart_disagreement),
        ("transport_metric_norm_drift", transport_metric_norm_drift),
    ])?;

    require_at_most(
        "cylindrical numerical transport versus exact flat oracle",
        cylindrical_transport_error,
        TRANSPORT_ORACLE_TOLERANCE,
    )?;
    require_at_most(
        "position tangent Cartesian/cylindrical disagreement",
        position_vector_disagreement,
        POSITION_VECTOR_TOLERANCE,
    )?;
    require_at_most(
        "velocity tangent Cartesian/cylindrical disagreement",
        velocity_vector_disagreement,
        VELOCITY_VECTOR_TOLERANCE,
    )?;
    for (name, value) in [
        ("position temporal frame magnitude", position_temporal_disagreement),
        ("position spatial frame magnitude", position_spatial_disagreement),
        ("velocity temporal frame magnitude", velocity_temporal_disagreement),
        ("velocity spatial frame magnitude", velocity_spatial_disagreement),
    ]
    {
        require_at_most(name, value, FRAME_MAGNITUDE_TOLERANCE)?;
    }
    require_at_most(
        "metric norm chart identity",
        metric_norm_chart_disagreement,
        METRIC_NORM_TOLERANCE,
    )?;
    require_at_most(
        "metric norm preservation under cylindrical transport",
        transport_metric_norm_drift,
        TRANSPORT_ORACLE_TOLERANCE,
    )?;

    print_common_header("geometric local-state diagnostic across flat-space charts");
    println!("# evidence: numerical diagnostic cross-checked against exact flat transport");
    println!("# claim boundary: NOT a proof of covariance/invariance; NOT wired into adaptive acceptance");
    println!(
        "# log settings: step={LOG_STEP}, jacobian_step={JACOBIAN_STEP}, tolerance={LOG_TOLERANCE}, max_iterations={LOG_MAX_ITERATIONS}"
    );
    println!("# transport_substeps={TRANSPORT_SUBSTEPS}, timelike_floor={TIMELIKE_FLOOR}");
    println!("quantity,cartesian,cylindrical,absolute_disagreement");
    println!(
        "position_temporal,{:.12e},{:.12e},{position_temporal_disagreement:.12e}",
        cartesian_position_frame.temporal, cylindrical_position_frame.temporal
    );
    println!(
        "position_spatial,{:.12e},{:.12e},{position_spatial_disagreement:.12e}",
        cartesian_position_frame.spatial, cylindrical_position_frame.spatial
    );
    println!(
        "velocity_temporal,{:.12e},{:.12e},{velocity_temporal_disagreement:.12e}",
        cartesian_velocity_frame.temporal, cylindrical_velocity_frame.temporal
    );
    println!(
        "velocity_spatial,{:.12e},{:.12e},{velocity_spatial_disagreement:.12e}",
        cartesian_velocity_frame.spatial, cylindrical_velocity_frame.spatial
    );
    println!(
        "metric_norm_low,{low_norm_cartesian:.12e},{low_norm_cylindrical:.12e},{:.12e}",
        (low_norm_cartesian - low_norm_cylindrical).abs()
    );
    println!(
        "metric_norm_high,{high_norm_cartesian:.12e},{high_norm_cylindrical:.12e},{:.12e}",
        (high_norm_cartesian - high_norm_cylindrical).abs()
    );
    println!("# position_vector_disagreement={position_vector_disagreement:.12e}");
    println!("# velocity_vector_disagreement={velocity_vector_disagreement:.12e}");
    println!("# cylindrical_transport_exact_oracle_error={cylindrical_transport_error:.12e}");
    println!("# transported_low_metric_norm_drift={transport_metric_norm_drift:.12e}");
    println!("# oracle_check: geometric_flat_cross_chart_diagnostic=true");
    Ok(())
}

fn subtract<const D: usize>(left: [f64; D], right: [f64; D]) -> [f64; D] {
    let mut difference = [0.0_f64; D];
    for (slot, (left_value, right_value)) in difference.iter_mut().zip(left.into_iter().zip(right))
    {
        *slot = left_value - right_value;
    }
    difference
}

fn max_abs_difference<const D: usize>(left: &[f64; D], right: &[f64; D]) -> f64 {
    left.iter()
        .zip(right.iter())
        .fold(0.0_f64, |maximum, (left_value, right_value)| {
            maximum.max((left_value - right_value).abs())
        })
}

fn metric_contraction<const D: usize>(
    metric: &[[f64; D]; D],
    left: &[f64; D],
    right: &[f64; D],
) -> f64 {
    let mut value = 0.0_f64;
    for mu in 0..D
    {
        for nu in 0..D
        {
            value += metric[mu][nu] * left[mu] * right[nu];
        }
    }
    value
}

fn require_at_most(name: &str, value: f64, tolerance: f64) -> Result<(), String> {
    if value <= tolerance
    {
        Ok(())
    }
    else
    {
        Err(format!(
            "{name} exceeded tolerance: value={value:.12e}, tolerance={tolerance:.12e}"
        ))
    }
}

fn stringify<E: std::fmt::Display>(error: E) -> String {
    error.to_string()
}

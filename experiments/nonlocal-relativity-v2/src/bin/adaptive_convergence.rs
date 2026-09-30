//! Phase 9 experiment: adaptive-tolerance convergence tables for both adaptive
//! controllers (embedded Heun-Euler and step-doubling), measured against an
//! independent fine fixed-step reference.
//!
//! For each error tolerance it reports the accepted-step count and the endpoint
//! coordinate error against the reference, plus the error-reduction ratio
//! between consecutive tolerances (a self-consistency trend, not a proof of a
//! convergence order). The reference is a fine fixed-step run of the matching
//! method family — a numerical reference, not an exact solution.

#![forbid(unsafe_code)]

use nonlocal_relativity_experiments::{
    circular_schwarzschild_state, euclidean_distance,
    oracle::{EvidenceClass, strictly_decreasing},
    print_common_header, require_finite,
};
use scirust_nonlocal_relativity::{
    AdaptiveNonlocalConfig, CaputoCoordinateMemory, CompleteUniformHistory, HeunPeceStepper,
    IdentityHistoryTransport, NonlocalConfig, NonlocalSimulationPolicy,
    NonuniformCaputoCoordinateMemory, SemiImplicitEulerStepper, WorldlineState,
    simulate_nonlocal_worldline_adaptive, simulate_nonlocal_worldline_adaptive_with_stepper,
    simulate_nonlocal_worldline_with_policy,
};
use scirust_relativity::Schwarzschild;

const MASS: f64 = 1.0;
const RADIUS: f64 = 10.0;
const ALPHA: f64 = 0.55;
const COUPLING: f64 = 0.02;
const TARGET: f64 = 0.8;
const TOLERANCES: [f64; 5] = [1.0e-5, 1.0e-6, 1.0e-7, 1.0e-8, 1.0e-9];

fn main() -> Result<(), String> {
    let background =
        Schwarzschild::try_new(MASS).ok_or_else(|| "invalid Schwarzschild mass".to_string())?;
    let mut initial = circular_schwarzschild_state(MASS, RADIUS);
    initial.velocity[1] = -0.01;

    let heun_reference = fine_reference_heun(&background, initial)?;
    let euler_reference = fine_reference_euler(&background, initial)?;

    print_common_header("adaptive-tolerance convergence (both controllers)");
    println!("# background: Schwarzschild M = {MASS}, exterior near-circular orbit r0 = {RADIUS}");
    println!("# alpha = {ALPHA}, kappa = {COUPLING}, target affine parameter = {TARGET}");
    println!(
        "# embedded reference: fixed-step HeunPeceStepper h = 5e-4; stepdoubling reference: fixed-step SemiImplicitEuler h = 5e-4"
    );
    println!(
        "# error_ratio = previous_row_error / this_row_error (trend only; not a proof of convergence order)"
    );
    println!("controller,tolerance,accepted_steps,endpoint_coord_err,error_ratio");

    run_controller(
        "embedded_heun_euler",
        &background,
        initial,
        &heun_reference,
        |bg, init, cfg| {
            simulate_nonlocal_worldline_adaptive(bg, init, cfg)
                .map(|t| (t.len() - 1, *t.final_state().unwrap()))
                .map_err(|e| e.to_string())
        },
    )?;

    run_controller(
        "step_doubling",
        &background,
        initial,
        &euler_reference,
        |bg, init, cfg| {
            simulate_nonlocal_worldline_adaptive_with_stepper(bg, init, cfg)
                .map(|t| (t.len() - 1, *t.final_state().unwrap()))
                .map_err(|e| e.to_string())
        },
    )?;

    println!("# interpretation: the tested sweep improves the endpoint error overall while");
    println!("# allowing explicitly counted, bit-identical endpoint/count plateaus. The strict");
    println!("# trend is reported separately; a plateau is not evidence of strict convergence.");
    println!("# This is a bounded fine-grid self-consistency check, not a convergence proof or");
    println!("# validation of the underlying phenomenological model.");
    Ok(())
}

fn run_controller<F>(
    label: &str,
    background: &Schwarzschild,
    initial: WorldlineState<4>,
    reference_final: &WorldlineState<4>,
    run: F,
) -> Result<(), String>
where
    F: Fn(
        &Schwarzschild,
        WorldlineState<4>,
        AdaptiveNonlocalConfig,
    ) -> Result<(usize, WorldlineState<4>), String>,
{
    let mut previous_error: Option<f64> = None;
    let mut errors = Vec::with_capacity(TOLERANCES.len());
    let mut step_counts = Vec::with_capacity(TOLERANCES.len());
    let mut endpoints = Vec::with_capacity(TOLERANCES.len());
    for tolerance in TOLERANCES
    {
        let config = AdaptiveNonlocalConfig::new(
            ALPHA, COUPLING, 0.02, 0.000005, 0.05, tolerance, 1.0e-8, TARGET, 200_000, 60,
        )
        .map_err(|error| error.to_string())?;
        let (accepted_steps, final_state) = run(background, initial, config)?;
        let error = euclidean_distance(&final_state.coordinates, &reference_final.coordinates);
        require_finite(&[("endpoint_coord_err", error)])?;
        let ratio = previous_error.map_or(f64::NAN, |previous| previous / error);
        if ratio.is_finite()
        {
            println!("{label},{tolerance:.0e},{accepted_steps},{error:.6e},{ratio:.3}");
        }
        else
        {
            println!("{label},{tolerance:.0e},{accepted_steps},{error:.6e},NA");
        }
        previous_error = Some(error);
        errors.push(error);
        step_counts.push(accepted_steps);
        endpoints.push(final_state);
    }

    let strict = strictly_decreasing(&errors).map_err(|error| error.to_string())?;
    // Preserve the failed strict-trend hypothesis in the output, even when the
    // distinct, explicitly versioned plateau-aware regression contract passes.
    println!(
        "# oracle_check: {label}_endpoint_error_strictly_decreasing={strict},evidence={}",
        EvidenceClass::FineGridReference.as_str()
    );
    let plateaus = verify_refinement_rows(&errors, &step_counts, &endpoints)
        .map_err(|error| format!("{label}: {error}"))?;
    println!(
        "# oracle_check: contract=adaptive-refinement/v2,controller={label},nonincreasing_error=true,overall_error_improved=true,nondecreasing_steps=true,overall_step_count_increased=true,identical_endpoint_count_plateaus={plateaus},evidence={}",
        EvidenceClass::FineGridReference.as_str()
    );
    Ok(())
}

fn fine_reference_heun(
    background: &Schwarzschild,
    initial: WorldlineState<4>,
) -> Result<WorldlineState<4>, String> {
    let step = 0.0005;
    let steps = (TARGET / step).round() as usize;
    let config = NonlocalConfig::new(ALPHA, COUPLING, step, steps, 1.0e-8)
        .map_err(|error| error.to_string())?;
    let trajectory = simulate_nonlocal_worldline_with_policy(
        background,
        initial,
        config,
        NonlocalSimulationPolicy::new(
            CompleteUniformHistory::<4>::with_capacity(steps + 1),
            CaputoCoordinateMemory,
            IdentityHistoryTransport,
            HeunPeceStepper,
        ),
    )
    .map_err(|error| error.to_string())?;
    trajectory
        .final_state()
        .copied()
        .ok_or_else(|| "empty reference".to_string())
}

fn fine_reference_euler(
    background: &Schwarzschild,
    initial: WorldlineState<4>,
) -> Result<WorldlineState<4>, String> {
    let step = 0.0005;
    let steps = (TARGET / step).round() as usize;
    let config = NonlocalConfig::new(ALPHA, COUPLING, step, steps, 1.0e-8)
        .map_err(|error| error.to_string())?;
    let trajectory = simulate_nonlocal_worldline_with_policy(
        background,
        initial,
        config,
        NonlocalSimulationPolicy::new(
            CompleteUniformHistory::<4>::with_capacity(steps + 1),
            NonuniformCaputoCoordinateMemory,
            IdentityHistoryTransport,
            SemiImplicitEulerStepper,
        ),
    )
    .map_err(|error| error.to_string())?;
    trajectory
        .final_state()
        .copied()
        .ok_or_else(|| "empty reference".to_string())
}

/// Validate a bounded refinement sweep without turning repeated output into
/// a strict convergence claim. No tolerance is applied to an error increase.
/// Equal errors are admitted only with the same accepted-step count and
/// bit-identical complete endpoints; that does not assert identical paths.
fn verify_refinement_rows(
    errors: &[f64],
    step_counts: &[usize],
    endpoints: &[WorldlineState<4>],
) -> Result<usize, String> {
    if errors.len() < 2 || errors.len() != step_counts.len() || errors.len() != endpoints.len()
    {
        return Err("refinement evidence requires at least two aligned rows".to_string());
    }
    for index in 0..errors.len()
    {
        let state = &endpoints[index];
        if !errors[index].is_finite()
            || errors[index] < 0.0
            || step_counts[index] == 0
            || state
                .coordinates
                .iter()
                .chain(state.velocity.iter())
                .any(|value| !value.is_finite())
        {
            return Err(format!("invalid refinement evidence at row {index}"));
        }
    }
    let mut plateaus = 0;
    for index in 1..errors.len()
    {
        if errors[index] > errors[index - 1]
        {
            return Err(format!("endpoint error increased at row {index}"));
        }
        if step_counts[index] < step_counts[index - 1]
        {
            return Err(format!("accepted-step count decreased at row {index}"));
        }
        if errors[index] == errors[index - 1]
        {
            let previous = &endpoints[index - 1];
            let current = &endpoints[index];
            if step_counts[index] != step_counts[index - 1]
                || previous.coordinates.map(f64::to_bits) != current.coordinates.map(f64::to_bits)
                || previous.velocity.map(f64::to_bits) != current.velocity.map(f64::to_bits)
            {
                return Err(format!(
                    "equal errors lack identical endpoint/count evidence at row {index}"
                ));
            }
            plateaus += 1;
        }
    }
    let last = errors.len() - 1;
    if errors[last] >= errors[0] || step_counts[last] <= step_counts[0]
    {
        return Err(
            "the sweep supplies no combined overall error improvement and step-count increase"
                .to_string(),
        );
    }
    Ok(plateaus)
}

#[cfg(test)]
mod refinement_evidence_tests {
    use super::{WorldlineState, strictly_decreasing, verify_refinement_rows};

    fn endpoint(value: f64) -> WorldlineState<4> {
        WorldlineState::new([value, 10.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0])
    }

    #[test]
    fn repeated_endpoint_is_counted_without_claiming_a_strict_trend() {
        let errors = [3.0, 3.0, 1.0];
        let states = [endpoint(3.0), endpoint(3.0), endpoint(1.0)];
        assert_eq!(
            verify_refinement_rows(&errors, &[17, 17, 40], &states),
            Ok(1)
        );
        assert!(!strictly_decreasing(&errors).unwrap());
        assert_eq!(
            verify_refinement_rows(
                &[3.0, 2.0, 1.0],
                &[17, 20, 40],
                &[endpoint(3.0), endpoint(2.0), endpoint(1.0)]
            ),
            Ok(0)
        );
    }

    #[test]
    fn even_one_ulp_of_error_regression_is_rejected() {
        let increased = f64::from_bits(3.0_f64.to_bits() + 1);
        assert!(
            verify_refinement_rows(
                &[3.0, increased, 1.0],
                &[17, 18, 40],
                &[endpoint(3.0), endpoint(increased), endpoint(1.0)]
            )
            .is_err()
        );
    }

    #[test]
    fn equal_errors_require_both_identical_endpoint_and_count() {
        let errors = [3.0, 3.0, 1.0];
        let states = [endpoint(3.0), endpoint(3.0), endpoint(1.0)];
        assert!(verify_refinement_rows(&errors, &[17, 18, 40], &states).is_err());
        let mut changed = states;
        changed[1].coordinates[0] = f64::from_bits(3.0_f64.to_bits() + 1);
        assert!(verify_refinement_rows(&errors, &[17, 17, 40], &changed).is_err());
        changed = states;
        changed[1].velocity[1] = -0.0;
        assert!(verify_refinement_rows(&errors, &[17, 17, 40], &changed).is_err());
    }

    #[test]
    fn no_progress_or_reduced_step_count_is_not_refinement_evidence() {
        assert!(
            verify_refinement_rows(&[3.0, 3.0], &[17, 17], &[endpoint(3.0), endpoint(3.0)])
                .is_err()
        );
        assert!(
            verify_refinement_rows(&[3.0, 1.0], &[17, 17], &[endpoint(3.0), endpoint(1.0)])
                .is_err()
        );
        assert!(
            verify_refinement_rows(
                &[3.0, 2.0, 1.0],
                &[17, 16, 40],
                &[endpoint(3.0), endpoint(2.0), endpoint(1.0)]
            )
            .is_err()
        );
    }

    #[test]
    fn malformed_or_nonfinite_evidence_fails_closed() {
        assert!(verify_refinement_rows(&[], &[], &[]).is_err());
        assert!(verify_refinement_rows(&[1.0], &[1], &[endpoint(1.0)]).is_err());
        assert!(
            verify_refinement_rows(&[3.0, 1.0], &[17], &[endpoint(3.0), endpoint(1.0)]).is_err()
        );
        for invalid in [f64::NAN, f64::INFINITY, -1.0]
        {
            assert!(
                verify_refinement_rows(&[invalid, 1.0], &[17, 40], &[endpoint(3.0), endpoint(1.0)])
                    .is_err()
            );
        }
        assert!(
            verify_refinement_rows(&[3.0, 1.0], &[0, 40], &[endpoint(3.0), endpoint(1.0)]).is_err()
        );
        assert!(
            verify_refinement_rows(&[3.0, 1.0], &[17, 40], &[endpoint(f64::NAN), endpoint(1.0)])
                .is_err()
        );
    }
}

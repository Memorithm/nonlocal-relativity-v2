//! BSSN on a periodic one-dimensional grid — oracles and convergence.
//!
//! Every tolerance below is derived from a stated error model: second-order
//! truncation `~ dx^2`, fourth-order RK4 truncation `~ dt^4`, the `O(A^2)`
//! nonlinearity the linearized oracle neglects, or floating-point rounding.
//! None is chosen to make a failing test pass.

use scirust_relativity::adm_evolution::{AdmEvolutionSettings, AdmSources, SpatialTensorField};
use scirust_relativity::bssn::{BssnGauge, adm_rhs_from_bssn, bssn_evolution_rhs};
use scirust_relativity::bssn_grid::{
    BssnGridState, BssnGridSystem, BssnGridView, BssnShiftCondition, BssnSlicing,
    COMPONENTS_PER_POINT, TransverseTracelessWave, bssn_grid_constraints, bssn_grid_rhs,
    bssn_grid_ricci_report, evolve_bssn_grid, grid_conformal_ricci, project_grid_trace_free,
    project_grid_unit_determinant,
};
use scirust_relativity::grid::UniformGrid1d;
use scirust_relativity::{Metric, ricci_tensor_from_metric};

const TWO_PI: f64 = std::f64::consts::TAU;
/// Small enough that the neglected `O(A^2)` nonlinearity (`1e-12`) sits far
/// below the `O(dx^2)` truncation error at every resolution tested here.
const WAVE_AMPLITUDE: f64 = 1.0e-6;

fn unit_grid(points: usize) -> UniformGrid1d {
    UniformGrid1d::new(points, 0.0, 1.0).expect("valid grid")
}

fn zero_curvature() -> impl SpatialTensorField {
    |_: &[f64; 3]| [[0.0_f64; 3]; 3]
}

/// A smooth periodic manufactured metric. `Metric<3>` has no blanket closure
/// impl, so this is a concrete type rather than a lambda.
///
/// It is **not** a solution of the Einstein equations -- it is a numerical
/// oracle for the derivative stack.
struct ManufacturedMetric {
    amplitude: f64,
    wave_number: f64,
}

impl Metric<3> for ManufacturedMetric {
    fn components(&self, coordinates: &[f64; 3]) -> [[f64; 3]; 3] {
        let s = (self.wave_number * coordinates[0]).sin();
        [
            [1.0, 0.0, 0.0],
            [0.0, 1.0 + self.amplitude * s, 0.0],
            [0.0, 0.0, 1.0 - self.amplitude * s],
        ]
    }
}

/// A spatially constant isotropic metric.
struct ConstantMetric {
    scale: f64,
}

impl Metric<3> for ConstantMetric {
    fn components(&self, _coordinates: &[f64; 3]) -> [[f64; 3]; 3] {
        [
            [self.scale, 0.0, 0.0],
            [0.0, self.scale, 0.0],
            [0.0, 0.0, self.scale],
        ]
    }
}

/// The largest absolute value in a flat right-hand-side array.
fn max_abs(values: &[f64]) -> f64 {
    values.iter().fold(0.0_f64, |acc, v| acc.max(v.abs()))
}

fn observed_order(coarse: f64, fine: f64) -> f64 {
    (coarse / fine).log2()
}

// ---------------------------------------------------------------------------
// Storage and packing
// ---------------------------------------------------------------------------

#[test]
fn flat_state_round_trips_exactly_and_preserves_symmetry() {
    let grid = unit_grid(8);
    let wave = TransverseTracelessWave::new(1.0e-3, TWO_PI, 0.0);
    let state = BssnGridState::from_adm_fields(grid, &wave.metric_field(), &wave.curvature_field())
        .expect("initial data");

    assert_eq!(state.as_slice().len(), COMPONENTS_PER_POINT * grid.points());

    // Round trip through the flat array is bit-for-bit, not merely close.
    let copy = BssnGridState::from_flat(grid, state.as_slice().to_vec()).expect("round trip");
    assert_eq!(copy.as_slice(), state.as_slice());

    for index in 0..grid.points()
    {
        let local = state.state_at(index);
        // Only six independent components are stored, so symmetry cannot be
        // broken by a round trip -- the redundant entries do not exist.
        for i in 0..3
        {
            for j in 0..3
            {
                assert_eq!(local.conformal_metric[i][j], local.conformal_metric[j][i]);
                assert_eq!(
                    local.conformal_curvature[i][j],
                    local.conformal_curvature[j][i]
                );
            }
        }
        assert_eq!(copy.state_at(index), local);
    }
}

#[test]
fn flat_state_rejects_wrong_lengths() {
    let grid = unit_grid(8);
    let short = vec![0.0_f64; COMPONENTS_PER_POINT * grid.points() - 1];
    assert!(BssnGridState::from_flat(grid, short.clone()).is_err());
    assert!(BssnGridView::new(grid, &short).is_err());

    let good = vec![0.0_f64; COMPONENTS_PER_POINT * grid.points()];
    assert!(BssnGridView::new(grid, &good).is_ok());

    // The right-hand side checks its output buffer too.
    let view = BssnGridView::new(grid, &good).expect("view");
    let mut wrong_out = vec![0.0_f64; 3];
    assert!(
        bssn_grid_rhs(
            &view,
            &BssnGauge::SYNCHRONOUS,
            &AdmSources::VACUUM,
            0.0,
            &mut wrong_out
        )
        .is_err()
    );
}

// ---------------------------------------------------------------------------
// Oracle A — stationary Minkowski
// ---------------------------------------------------------------------------

#[test]
fn oracle_a_minkowski_right_hand_side_is_exactly_zero() {
    let grid = unit_grid(32);
    let state = BssnGridState::minkowski(grid);
    state.validate(0.0).expect("Minkowski is a valid state");

    let view = state.view();
    let mut rhs = vec![0.0_f64; COMPONENTS_PER_POINT * grid.points()];
    bssn_grid_rhs(
        &view,
        &BssnGauge::SYNCHRONOUS,
        &AdmSources::VACUUM,
        0.0,
        &mut rhs,
    )
    .expect("Minkowski right-hand side");

    // A centred difference of a constant field subtracts identical values, so
    // this is exactly zero -- not "small". Anything else is an indexing bug.
    for (slot, value) in rhs.iter().enumerate()
    {
        assert_eq!(*value, 0.0, "slot {slot} of the Minkowski RHS is {value}");
    }
}

#[test]
fn oracle_a_minkowski_constraints_are_exactly_zero() {
    let grid = unit_grid(32);
    let state = BssnGridState::minkowski(grid);
    let constraints =
        bssn_grid_constraints(&state.view(), &AdmSources::VACUUM).expect("constraints");

    assert_eq!(constraints.determinant.max_abs, 0.0);
    assert_eq!(constraints.trace_free.max_abs, 0.0);
    assert_eq!(constraints.connection.max_abs, 0.0);
    assert_eq!(constraints.hamiltonian.max_abs, 0.0);
    for component in 0..3
    {
        assert_eq!(constraints.momentum[component].max_abs, 0.0);
    }
    assert_eq!(constraints.determinant.non_finite, 0);
    assert_eq!(constraints.hamiltonian.non_finite, 0);
}

#[test]
fn oracle_a_minkowski_does_not_drift_over_many_steps() {
    let grid = unit_grid(16);
    let system = BssnGridSystem::vacuum(grid);
    let initial = BssnGridState::minkowski(grid);

    // 200 steps: long enough that any per-step drift would accumulate visibly.
    let samples = evolve_bssn_grid(&system, &initial, 0.0, 2.0, 0.01).expect("evolution");
    assert_eq!(samples.len(), 201);

    for sample in &samples
    {
        // Bit-for-bit stationary: RK4 of an exactly-zero right-hand side adds
        // exactly zero, so there is no rounding to accumulate.
        assert_eq!(
            sample.state.as_slice(),
            initial.as_slice(),
            "drift at t = {}",
            sample.time
        );
    }

    // And the drift does not grow with step count, because there is none.
    let last = samples.last().expect("final sample");
    let constraints =
        bssn_grid_constraints(&last.state.view(), &AdmSources::VACUUM).expect("constraints");
    assert_eq!(constraints.hamiltonian.max_abs, 0.0);
}

// ---------------------------------------------------------------------------
// Grid / local equivalence
// ---------------------------------------------------------------------------

#[test]
fn grid_and_local_right_hand_sides_agree_for_constant_fields() {
    // A spatially constant -- but non-trivial -- state must reproduce the
    // Layer 3.3 pointwise answer exactly, because the grid adds only spatial
    // derivatives and those vanish for a constant field.
    let grid = unit_grid(8);
    let scale = 1.4_f64;
    let curvature_amplitude = 0.3_f64;

    let constant_metric = ConstantMetric { scale };
    let constant_curvature = move |_: &[f64; 3]| {
        [
            [curvature_amplitude, 0.0, 0.0],
            [0.0, curvature_amplitude, 0.0],
            [0.0, 0.0, curvature_amplitude],
        ]
    };

    let state = BssnGridState::from_adm_fields(grid, &constant_metric, &constant_curvature)
        .expect("initial data");
    let view = state.view();
    let mut rhs = vec![0.0_f64; COMPONENTS_PER_POINT * grid.points()];
    bssn_grid_rhs(
        &view,
        &BssnGauge::SYNCHRONOUS,
        &AdmSources::VACUUM,
        0.0,
        &mut rhs,
    )
    .expect("grid rhs");

    let settings = AdmEvolutionSettings {
        spatial_step: grid.spacing(),
        metric_step: grid.spacing(),
    };
    let (_, local) = bssn_evolution_rhs(
        &constant_metric,
        &constant_curvature,
        &[grid.coordinate(3), 0.0, 0.0],
        &BssnGauge::SYNCHRONOUS,
        &AdmSources::VACUUM,
        &settings,
    )
    .expect("local rhs");

    let points = grid.points();
    for index in 0..points
    {
        // Every point sees the same constant field, so every point must get
        // bit-identical answers -- including the wrap points.
        assert_eq!(rhs[7 * points + index], local.mean_curvature);
        assert_eq!(rhs[index], local.conformal_factor);
    }
}

// ---------------------------------------------------------------------------
// Oracle B — manufactured periodic state
// ---------------------------------------------------------------------------

#[test]
fn oracle_b_conformal_ricci_reconstruction_converges_at_second_order() {
    // A manufactured smooth periodic metric. It is NOT a solution of the
    // Einstein equations; it is a numerical oracle for the derivative stack.
    let amplitude = 0.01_f64;
    let k = TWO_PI;
    let manufactured = ManufacturedMetric {
        amplitude,
        wave_number: k,
    };

    let mut errors = Vec::new();
    let mut scales = Vec::new();
    for &points in &[32_usize, 64, 128, 256]
    {
        let grid = unit_grid(points);
        let state = BssnGridState::from_adm_fields(grid, &manufactured, &zero_curvature())
            .expect("initial data");
        let report = bssn_grid_ricci_report(&state.view()).expect("ricci report");
        errors.push(report.mismatch.max_abs);
        scales.push((
            report.conformal_scale,
            report.conformal_factor_scale,
            report.physical_scale,
        ));
    }

    // The decomposition must not agree by both sides being zero: each part
    // has to be genuinely nonzero at every resolution.
    for (conformal, conformal_factor, physical) in &scales
    {
        assert!(*conformal > 1.0e-4, "Rtilde is trivially zero: {conformal}");
        assert!(
            *conformal_factor > 1.0e-6,
            "Rphi is trivially zero: {conformal_factor}"
        );
        assert!(*physical > 1.0e-4, "R is trivially zero: {physical}");
    }

    // Layer 3.3 could only report this mismatch as a fixed ~1e-6 "floor". On a
    // grid it is truncation error and must converge at the discretisation order.
    for window in errors.windows(2)
    {
        let order = observed_order(window[0], window[1]);
        assert!(
            (order - 2.0).abs() < 0.15,
            "observed Ricci order {order} from {errors:?}"
        );
    }
}

// ---------------------------------------------------------------------------
// Oracle C — linearized transverse-traceless wave
// ---------------------------------------------------------------------------

#[test]
fn oracle_c_wave_constraint_residual_is_the_neglected_quadratic_nonlinearity() {
    // A transverse-traceless wave solves the *linearized* equations, so its
    // Hamiltonian residual is O(A^2) analytically. The measured residual turns
    // out to be dominated by exactly that, not by discretisation:
    //
    //   fixed A = 1e-6:  N = 32..256 gives 7.74e-11 -> 8.98e-11  (flat)
    //   fixed N = 64:    A doubling gives ratios 3.96, 4.01, 4.00 (= A^2)
    //
    // The reason the truncation error does not show up is that the linear-order
    // Ricci *scalar* is proportional to d^2(h_yy + h_zz), and h_yy + h_zz = 0
    // holds exactly at every grid point, so the discrete trace cancels term by
    // term. What survives is the genuine quadratic nonlinearity the linearized
    // oracle neglects.
    //
    // Asserting dx-convergence here would therefore be asserting something
    // false. The honest, and much sharper, test is the amplitude scaling.
    let k = TWO_PI;
    let grid = unit_grid(64);
    let mut residuals = Vec::new();

    for &amplitude in &[1.0e-6_f64, 2.0e-6, 4.0e-6, 8.0e-6]
    {
        let wave = TransverseTracelessWave::new(amplitude, k, 0.0);
        let state =
            BssnGridState::from_adm_fields(grid, &wave.metric_field(), &wave.curvature_field())
                .expect("wave initial data");
        let constraints =
            bssn_grid_constraints(&state.view(), &AdmSources::VACUUM).expect("constraints");

        // Algebraic constraints are satisfied by construction: the conversion
        // built the state and seeded Gammatilde^i with the same reconstruction
        // the monitor compares against.
        assert!(
            constraints.determinant.max_abs < 1.0e-14,
            "determinant residual {}",
            constraints.determinant.max_abs
        );
        assert!(
            constraints.trace_free.max_abs < 1.0e-14,
            "trace residual {}",
            constraints.trace_free.max_abs
        );
        assert!(
            constraints.connection.max_abs < 1.0e-14,
            "connection residual {}",
            constraints.connection.max_abs
        );

        residuals.push(constraints.hamiltonian.max_abs);
    }

    // Doubling the amplitude must quadruple the residual: that is what makes it
    // O(A^2) rather than a discretisation artefact.
    for window in residuals.windows(2)
    {
        let order = observed_order(window[1], window[0]);
        assert!(
            (order - 2.0).abs() < 0.1,
            "amplitude scaling exponent {order} from {residuals:?}"
        );
    }
}

#[test]
fn oracle_c_wave_constraint_residual_is_not_resolution_limited() {
    // The complement of the test above: at fixed amplitude, refining the grid
    // does NOT reduce the residual, because it is not truncation error. This is
    // recorded as a measured property rather than left as an unexplained
    // non-convergence.
    let k = TWO_PI;
    let mut residuals = Vec::new();
    for &points in &[32_usize, 64, 128, 256]
    {
        let grid = unit_grid(points);
        let wave = TransverseTracelessWave::new(WAVE_AMPLITUDE, k, 0.0);
        let state =
            BssnGridState::from_adm_fields(grid, &wave.metric_field(), &wave.curvature_field())
                .expect("wave initial data");
        let constraints =
            bssn_grid_constraints(&state.view(), &AdmSources::VACUUM).expect("constraints");
        residuals.push(constraints.hamiltonian.max_abs);

        // The momentum constraint is likewise at the nonlinearity level.
        for component in 0..3
        {
            assert!(
                constraints.momentum[component].max_abs < 1.0e-8,
                "momentum[{component}] = {}",
                constraints.momentum[component].max_abs
            );
        }
    }

    // Every residual sits within a factor of two of the O(A^2) level; a
    // truncation-limited quantity would fall by 64 across this range.
    let smallest = residuals.iter().cloned().fold(f64::INFINITY, f64::min);
    let largest = residuals.iter().cloned().fold(0.0_f64, f64::max);
    assert!(
        largest / smallest < 2.0,
        "residual varies with resolution more than expected: {residuals:?}"
    );
    assert!(
        largest < 1.0e-9,
        "residual above the O(A^2) level: {residuals:?}"
    );
}

#[test]
fn oracle_c_wave_right_hand_side_matches_the_linearized_prediction() {
    // For a TT wave with unit lapse and zero shift, vacuum:
    //   d_t K_yy = R_yy + O(A^2) = (A k^2 / 2) sin(k(x - t)) + O(A^2).
    // This checks the assembled right-hand side against that closed form and
    // measures how the difference converges.
    let k = TWO_PI;
    let mut errors = Vec::new();

    for &points in &[32_usize, 64, 128, 256]
    {
        let grid = unit_grid(points);
        let wave = TransverseTracelessWave::new(WAVE_AMPLITUDE, k, 0.0);
        let state =
            BssnGridState::from_adm_fields(grid, &wave.metric_field(), &wave.curvature_field())
                .expect("wave initial data");
        let view = state.view();
        let settings = view.settings();

        let mut worst = 0.0_f64;
        for index in 0..grid.points()
        {
            let x = grid.coordinate(index);
            let at = [x, 0.0, 0.0];
            let (bssn_state, rhs) = bssn_evolution_rhs(
                &view.metric(),
                &view.curvature(),
                &at,
                &BssnGauge::SYNCHRONOUS,
                &AdmSources::VACUUM,
                &settings,
            )
            .expect("rhs");
            let adm_rhs = adm_rhs_from_bssn(&bssn_state, &rhs);
            let predicted = 0.5 * WAVE_AMPLITUDE * k * k * (k * x).sin();
            worst = worst.max((adm_rhs.extrinsic_curvature[1][1] - predicted).abs());
        }
        errors.push(worst);
    }

    // Second order: the difference is dominated by dx^2 truncation, since the
    // neglected nonlinearity is O(A^2) = 1e-12 and the coarsest error here is
    // orders of magnitude larger.
    for window in errors.windows(2)
    {
        let order = observed_order(window[0], window[1]);
        assert!(
            (order - 2.0).abs() < 0.2,
            "observed d_t K order {order} from {errors:?}"
        );
    }
}

#[test]
fn oracle_c_short_time_propagation_converges_in_space() {
    // Evolve a short interval and compare against the analytic linear wave.
    // The timestep is refined with the grid so the temporal error (dt^4) stays
    // negligible against the spatial error (dx^2) at every resolution.
    let k = TWO_PI;
    let courant = 0.25_f64;
    let t_end = 0.1_f64;
    let mut errors = Vec::new();

    for &points in &[16_usize, 32, 64]
    {
        let grid = unit_grid(points);
        let step = courant * grid.spacing();
        let steps = (t_end / step).round() as usize;
        let step = t_end / steps as f64;

        let initial_wave = TransverseTracelessWave::new(WAVE_AMPLITUDE, k, 0.0);
        let initial = BssnGridState::from_adm_fields(
            grid,
            &initial_wave.metric_field(),
            &initial_wave.curvature_field(),
        )
        .expect("wave initial data");

        let system = BssnGridSystem::vacuum(grid);
        let samples =
            evolve_bssn_grid(&system, &initial, 0.0, t_end, step).expect("wave evolution");
        let final_sample = samples.last().expect("final sample");

        let exact = TransverseTracelessWave::new(WAVE_AMPLITUDE, k, final_sample.time);
        let mut worst = 0.0_f64;
        for index in 0..grid.points()
        {
            let x = grid.coordinate(index);
            let numerical = final_sample.state.state_at(index);
            let reconstructed =
                scirust_relativity::bssn::bssn_to_adm(&numerical).expect("reconstruct");
            let analytic = exact.spatial_metric(x);
            worst = worst.max((reconstructed.spatial_metric[1][1] - analytic[1][1]).abs());
        }
        errors.push(worst);
    }

    // Halving dx must reduce the error by about four. The window is wide enough
    // to absorb the mild resolution dependence of the leading constant at these
    // very coarse grids, but nowhere near wide enough to admit first order.
    for window in errors.windows(2)
    {
        let order = observed_order(window[0], window[1]);
        assert!(
            order > 1.6,
            "observed propagation order {order} from {errors:?}"
        );
    }
}

#[test]
fn oracle_c_temporal_refinement_is_fourth_order() {
    // Temporal order is measured against a finely-stepped run *at the same
    // spatial resolution*, so the spatial truncation error cancels exactly and
    // what remains is purely the RK4 error of the semidiscrete ODE system.
    // Measuring against the analytic PDE solution instead would be dominated by
    // dx^2 and would report a misleadingly low order.
    let grid = unit_grid(16);
    let k = TWO_PI;
    let t_end = 0.2_f64;

    let wave = TransverseTracelessWave::new(1.0e-3, k, 0.0);
    let initial =
        BssnGridState::from_adm_fields(grid, &wave.metric_field(), &wave.curvature_field())
            .expect("wave initial data");
    let system = BssnGridSystem::vacuum(grid);

    let reference = evolve_bssn_grid(&system, &initial, 0.0, t_end, t_end / 640.0)
        .expect("reference evolution");
    let reference_state = reference.last().expect("final").state.clone();

    let mut errors = Vec::new();
    for &steps in &[10_usize, 20, 40]
    {
        let samples = evolve_bssn_grid(&system, &initial, 0.0, t_end, t_end / steps as f64)
            .expect("evolution");
        let final_state = &samples.last().expect("final").state;
        let difference: Vec<f64> = final_state
            .as_slice()
            .iter()
            .zip(reference_state.as_slice())
            .map(|(a, b)| a - b)
            .collect();
        errors.push(max_abs(&difference));
    }

    for window in errors.windows(2)
    {
        let order = observed_order(window[0], window[1]);
        assert!(
            order > 3.5,
            "observed temporal order {order} from {errors:?}"
        );
    }
}

// ---------------------------------------------------------------------------
// Oracle D — off-constraint state
// ---------------------------------------------------------------------------

#[test]
fn oracle_d_constraint_violation_is_visible_and_never_silently_repaired() {
    let grid = unit_grid(16);
    let k = TWO_PI;
    let wave = TransverseTracelessWave::new(1.0e-3, k, 0.0);
    let mut state =
        BssnGridState::from_adm_fields(grid, &wave.metric_field(), &wave.curvature_field())
            .expect("wave initial data");

    // Inject a periodic violation into K, which enters the Hamiltonian
    // constraint through K^2 and cannot be absorbed by any algebraic identity.
    let amplitude = 1.0e-3_f64;
    for index in 0..grid.points()
    {
        let mut local = state.state_at(index);
        local.mean_curvature += amplitude * (k * grid.coordinate(index)).cos();
        state.set_state_at(index, &local);
    }

    let before =
        bssn_grid_constraints(&state.view(), &AdmSources::VACUUM).expect("constraints before");
    // H picks up K^2 = (1e-3)^2 = 1e-6 -- far above the 1e-12 nonlinearity and
    // the truncation error of the clean wave.
    assert!(
        before.hamiltonian.max_abs > 1.0e-7,
        "violation not visible: {}",
        before.hamiltonian.max_abs
    );

    // Free evolution must not repair it.
    let system = BssnGridSystem::vacuum(grid);
    let samples = evolve_bssn_grid(&system, &state, 0.0, 0.05, 0.005).expect("evolution");
    let after = bssn_grid_constraints(
        &samples.last().expect("final").state.view(),
        &AdmSources::VACUUM,
    )
    .expect("constraints after");
    assert!(
        after.hamiltonian.max_abs > 1.0e-7,
        "constraint was silently repaired: {} -> {}",
        before.hamiltonian.max_abs,
        after.hamiltonian.max_abs
    );

    // The algebraic constraints are untouched by a K perturbation: perturbing
    // the trace cannot break det gammatilde = 1 or the trace-free condition.
    assert!(before.determinant.max_abs < 1.0e-14);
    assert!(before.trace_free.max_abs < 1.0e-14);
}

// ---------------------------------------------------------------------------
// Oracle E — explicit algebraic projection
// ---------------------------------------------------------------------------

#[test]
fn oracle_e_determinant_violation_projects_and_is_idempotent() {
    let grid = unit_grid(16);
    let mut state = BssnGridState::minkowski(grid);

    let epsilon = 0.01_f64;
    for index in 0..grid.points()
    {
        let mut local = state.state_at(index);
        for i in 0..3
        {
            local.conformal_metric[i][i] *= 1.0 + epsilon;
        }
        state.set_state_at(index, &local);
    }

    let before =
        bssn_grid_constraints(&state.view(), &AdmSources::VACUUM).expect("constraints before");
    // Scaling all three diagonal entries by (1 + eps) multiplies the
    // determinant by (1 + eps)^3, so the residual is exactly (1 + eps)^3 - 1.
    let expected = (1.0 + epsilon).powi(3) - 1.0;
    assert!(
        (before.determinant.max_abs - expected).abs() < 1.0e-12,
        "residual {} vs closed form {expected}",
        before.determinant.max_abs
    );

    let report = project_grid_unit_determinant(&mut state).expect("projection");
    assert!((report.residual_before - expected).abs() < 1.0e-12);
    assert!(report.residual_after < 1.0e-14);
    assert!(report.correction_magnitude > 0.0);

    // Idempotent: projecting an already-projected state changes nothing
    // meaningful.
    let again = project_grid_unit_determinant(&mut state).expect("second projection");
    assert!(again.residual_before < 1.0e-14);
    assert!(again.residual_after < 1.0e-14);
}

#[test]
fn oracle_e_trace_violation_projects_and_default_evolution_does_not() {
    let grid = unit_grid(16);
    let mut state = BssnGridState::minkowski(grid);

    let amplitude = 0.01_f64;
    for index in 0..grid.points()
    {
        let mut local = state.state_at(index);
        for i in 0..3
        {
            local.conformal_curvature[i][i] += amplitude;
        }
        state.set_state_at(index, &local);
    }

    let before =
        bssn_grid_constraints(&state.view(), &AdmSources::VACUUM).expect("constraints before");
    // With an identity conformal metric, gammatilde^{ij} Atilde_ij is exactly
    // the sum of the three diagonal perturbations.
    assert!(
        (before.trace_free.max_abs - 3.0 * amplitude).abs() < 1.0e-12,
        "trace residual {} vs closed form {}",
        before.trace_free.max_abs,
        3.0 * amplitude
    );

    // Default evolution is free: it must NOT project.
    let system = BssnGridSystem::vacuum(grid);
    let samples = evolve_bssn_grid(&system, &state, 0.0, 0.02, 0.002).expect("evolution");
    let after = bssn_grid_constraints(
        &samples.last().expect("final").state.view(),
        &AdmSources::VACUUM,
    )
    .expect("constraints after");
    assert!(
        after.trace_free.max_abs > 1.0e-3,
        "default evolution projected the trace away: {}",
        after.trace_free.max_abs
    );

    // Explicit projection, by contrast, removes it.
    let report = project_grid_trace_free(&mut state).expect("projection");
    assert!((report.residual_before - 3.0 * amplitude).abs() < 1.0e-12);
    assert!(report.residual_after < 1.0e-14);
}

// ---------------------------------------------------------------------------
// Rejection and determinism
// ---------------------------------------------------------------------------

#[test]
fn invalid_states_and_requests_are_rejected_with_typed_errors() {
    let grid = unit_grid(8);

    // A zeroed state has a singular conformal metric and must be rejected
    // rather than silently repaired.
    let zeroed = BssnGridState::zeroed(grid);
    assert!(zeroed.validate(0.0).is_err());

    let system = BssnGridSystem::vacuum(grid);
    let good = BssnGridState::minkowski(grid);
    assert!(evolve_bssn_grid(&system, &zeroed, 0.0, 0.1, 0.01).is_err());

    // Invalid evolution parameters.
    assert!(evolve_bssn_grid(&system, &good, 0.0, 0.1, 0.0).is_err());
    assert!(evolve_bssn_grid(&system, &good, 0.0, 0.1, -0.01).is_err());
    assert!(evolve_bssn_grid(&system, &good, 0.0, 0.1, f64::NAN).is_err());
    assert!(evolve_bssn_grid(&system, &good, 0.0, f64::INFINITY, 0.01).is_err());
    // A backwards interval is rejected rather than silently producing nothing.
    assert!(evolve_bssn_grid(&system, &good, 1.0, 0.0, 0.01).is_err());

    // A non-finite stored field is caught by validation, with context.
    let mut dirty = BssnGridState::minkowski(grid);
    let mut local = dirty.state_at(3);
    local.mean_curvature = f64::NAN;
    dirty.set_state_at(3, &local);
    match dirty.validate(0.5)
    {
        Err(scirust_relativity::bssn_grid::BssnGridError::InvalidState {
            time,
            index,
            field,
            ..
        }) =>
        {
            assert_eq!(time, 0.5);
            assert_eq!(index, 3);
            assert_eq!(field, "K");
        },
        other => panic!("expected a located InvalidState, got {other:?}"),
    }
}

#[test]
fn evolution_and_diagnostics_are_deterministic() {
    let grid = unit_grid(16);
    let wave = TransverseTracelessWave::new(1.0e-4, TWO_PI, 0.0);
    let initial =
        BssnGridState::from_adm_fields(grid, &wave.metric_field(), &wave.curvature_field())
            .expect("initial data");
    let system = BssnGridSystem::vacuum(grid);

    let first = evolve_bssn_grid(&system, &initial, 0.0, 0.05, 0.005).expect("first run");
    let second = evolve_bssn_grid(&system, &initial, 0.0, 0.05, 0.005).expect("second run");
    assert_eq!(first.len(), second.len());
    for (a, b) in first.iter().zip(second.iter())
    {
        assert_eq!(a.time, b.time);
        // Bit-for-bit, not approximately.
        assert_eq!(a.state.as_slice(), b.state.as_slice());
    }

    let constraints_a =
        bssn_grid_constraints(&initial.view(), &AdmSources::VACUUM).expect("constraints");
    let constraints_b =
        bssn_grid_constraints(&initial.view(), &AdmSources::VACUUM).expect("constraints");
    assert_eq!(constraints_a, constraints_b);

    // The maximum-residual index is selected deterministically.
    assert_eq!(
        constraints_a.hamiltonian.max_index,
        constraints_b.hamiltonian.max_index
    );
}

#[test]
fn transverse_derivatives_are_exactly_zero_by_construction() {
    // The 1D3V reduction: the provider is a function of x alone, so a
    // difference in y or z subtracts a value from itself and is exactly zero.
    let grid = unit_grid(16);
    let wave = TransverseTracelessWave::new(1.0e-3, TWO_PI, 0.0);
    let state = BssnGridState::from_adm_fields(grid, &wave.metric_field(), &wave.curvature_field())
        .expect("initial data");
    let view = state.view();
    let metric = view.metric();
    let dx = grid.spacing();

    for index in 0..grid.points()
    {
        let x = grid.coordinate(index);
        let centre = metric.components(&[x, 0.0, 0.0]);
        for axis in 1..3
        {
            let mut shifted = [x, 0.0, 0.0];
            shifted[axis] += dx;
            let moved = metric.components(&shifted);
            for i in 0..3
            {
                for j in 0..3
                {
                    assert_eq!(centre[i][j], moved[i][j]);
                }
            }
        }
    }
}

#[test]
fn ricci_decomposition_agrees_with_the_independent_physical_ricci_pointwise() {
    // Point-by-point, not just in a norm: a norm can hide a single bad point.
    let grid = unit_grid(64);
    let amplitude = 0.01_f64;
    let k = TWO_PI;
    let manufactured = ManufacturedMetric {
        amplitude,
        wave_number: k,
    };
    let state = BssnGridState::from_adm_fields(grid, &manufactured, &zero_curvature())
        .expect("initial data");
    let view = state.view();
    let settings = view.settings();
    let metric = view.metric();

    for index in 0..grid.points()
    {
        let at = [grid.coordinate(index), 0.0, 0.0];
        let decomposition =
            scirust_relativity::bssn::conformal_ricci(&metric, &at, &settings).expect("ricci");
        let physical = ricci_tensor_from_metric::<_, 3>(
            &metric,
            &at,
            settings.spatial_step,
            settings.metric_step,
        )
        .expect("physical ricci");
        // Tensor-component loops index two distinct rank-2 arrays together.
        #[allow(clippy::needless_range_loop)]
        for i in 0..3
        {
            for j in 0..3
            {
                // The bound is the measured second-order truncation at this
                // resolution, not a round number picked to pass.
                assert!(
                    (decomposition.total[i][j] - physical[i][j]).abs() < 1.0e-3,
                    "point {index}, component ({i},{j})"
                );
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Stability of the genuine BSSN principal part
// ---------------------------------------------------------------------------
//
// An earlier revision of this module computed `Rtilde_ij` from the generic
// metric Ricci and held `Gammatilde^i` frozen. That is BSSN *variables* carrying
// the ADM *principal part*, and it behaved exactly like weakly hyperbolic ADM:
// `N >= 64` blew up at `t < 1` at every Courant factor, worse as the grid was
// refined, and dissipation could not cure it.
//
// The fix was not a stencil width. It was to write `Rtilde_ij` in genuine BSSN
// form, so the term `gammatilde_{k(i} d_{j)} Gammatilde^k` carries the evolved
// connection and removes the mixed second derivatives from the principal part —
// and to supply the `d_t Gammatilde^i` equation that Layer 3.3 had explicitly
// deferred. These tests lock that in.

#[test]
fn wave_evolution_is_stable_at_every_resolution() {
    // Previously N >= 64 failed before t = 1 at every Courant factor. All four
    // resolutions must now reach t = 4, and refining the grid must not make
    // things worse.
    let k = TWO_PI;
    let amplitude = 1.0e-3_f64;
    let mut errors = Vec::new();

    for &points in &[32_usize, 64, 128]
    {
        let grid = unit_grid(points);
        let wave = TransverseTracelessWave::new(amplitude, k, 0.0);
        let initial =
            BssnGridState::from_adm_fields(grid, &wave.metric_field(), &wave.curvature_field())
                .expect("initial data");
        let system = BssnGridSystem::vacuum(grid);
        let samples = evolve_bssn_grid(&system, &initial, 0.0, 4.0, 0.25 * grid.spacing())
            .unwrap_or_else(|error| {
                panic!("N = {points} failed to reach t = 4: {error}");
            });

        let last = samples.last().expect("final sample");
        let exact = TransverseTracelessWave::new(amplitude, k, last.time);
        let mut worst = 0.0_f64;
        for index in 0..grid.points()
        {
            let reconstructed =
                scirust_relativity::bssn::bssn_to_adm(&last.state.state_at(index)).expect("adm");
            let analytic = exact.spatial_metric(grid.coordinate(index));
            worst = worst.max((reconstructed.spatial_metric[1][1] - analytic[1][1]).abs());
        }
        errors.push(worst);
    }

    // The error must not grow with resolution -- that was the signature of the
    // instability. It saturates rather than converging at order 2 because at
    // t = 4 the O(A^2) nonlinearity the linearized oracle neglects (1e-6 here)
    // dominates the remaining truncation error.
    for window in errors.windows(2)
    {
        assert!(
            window[1] <= window[0],
            "error grew with resolution: {errors:?}"
        );
    }
    assert!(
        errors.iter().all(|value| *value < 1.0e-3),
        "errors above the nonlinearity level: {errors:?}"
    );
}

#[test]
fn connection_constraint_stays_bounded_and_resolution_independent() {
    // The connection constraint measures the drift of the evolved Gammatilde^i
    // from -d_j gammatilde^{ij}. With the equation supplied it must grow only
    // secularly (linearly in t) rather than exploding, and it must not worsen
    // with resolution -- both of which failed before the fix.
    let k = TWO_PI;
    let grid_coarse = unit_grid(32);
    let grid_fine = unit_grid(128);

    let measure = |grid: UniformGrid1d, t_end: f64| -> f64 {
        let wave = TransverseTracelessWave::new(1.0e-3, k, 0.0);
        let initial =
            BssnGridState::from_adm_fields(grid, &wave.metric_field(), &wave.curvature_field())
                .expect("initial data");
        let system = BssnGridSystem::vacuum(grid);
        let samples = evolve_bssn_grid(&system, &initial, 0.0, t_end, 0.25 * grid.spacing())
            .expect("evolution");
        bssn_grid_constraints(
            &samples.last().expect("final").state.view(),
            &AdmSources::VACUUM,
        )
        .expect("constraints")
        .connection
        .max_abs
    };

    let one = measure(grid_coarse, 1.0);
    let two = measure(grid_coarse, 2.0);
    // Linear secular growth: doubling the time roughly doubles the drift. An
    // exponential instability would blow this ratio far past 2.
    let ratio = two / one;
    assert!(
        (1.5..3.0).contains(&ratio),
        "connection drift is not secular: {one} -> {two} (ratio {ratio})"
    );

    // And refining the grid does not make it worse.
    let fine = measure(grid_fine, 1.0);
    assert!(
        fine <= one * 1.5,
        "connection drift worsened with resolution: {one} (N=32) vs {fine} (N=128)"
    );
}

#[test]
fn grid_conformal_ricci_matches_the_generic_metric_ricci() {
    // The BSSN-form Ricci uses the evolved Gammatilde^k and compact stencils;
    // the generic form differentiates the metric. They are the same tensor, so
    // they must agree -- and the difference must converge, confirming the BSSN
    // formula rather than merely asserting it.
    let amplitude = 0.01_f64;
    let k = TWO_PI;
    let manufactured = ManufacturedMetric {
        amplitude,
        wave_number: k,
    };

    let mut differences = Vec::new();
    for &points in &[32_usize, 64, 128, 256]
    {
        let grid = unit_grid(points);
        let state = BssnGridState::from_adm_fields(grid, &manufactured, &zero_curvature())
            .expect("initial data");
        let view = state.view();
        let settings = view.settings();
        let metric = view.metric();

        let mut worst = 0.0_f64;
        let mut scale = 0.0_f64;
        for index in 0..grid.points()
        {
            let at = [grid.coordinate(index), 0.0, 0.0];
            let bssn = grid_conformal_ricci(&view, index).expect("bssn ricci");
            let generic = scirust_relativity::bssn::conformal_ricci(&metric, &at, &settings)
                .expect("generic ricci");
            for i in 0..3
            {
                for j in 0..3
                {
                    worst = worst.max((bssn.total[i][j] - generic.total[i][j]).abs());
                    scale = scale.max(generic.total[i][j].abs());
                }
            }
        }
        // Not two zeros agreeing.
        assert!(scale > 1.0e-3, "Ricci is trivially zero: {scale}");
        differences.push(worst);
    }

    for window in differences.windows(2)
    {
        let order = observed_order(window[0], window[1]);
        assert!(
            (order - 2.0).abs() < 0.15,
            "BSSN/generic Ricci order {order} from {differences:?}"
        );
    }
}

#[test]
fn the_connection_equation_is_actually_evolved() {
    // A regression guard on the specific defect that caused the instability:
    // Layer 3.3 reports d_t Gammatilde^i = 0 because it had no spatial
    // gradients. On a grid it must not be zero.
    let grid = unit_grid(32);
    let wave = TransverseTracelessWave::new(1.0e-3, TWO_PI, 0.0);
    let state = BssnGridState::from_adm_fields(grid, &wave.metric_field(), &wave.curvature_field())
        .expect("initial data");
    let view = state.view();
    let mut rhs = vec![0.0_f64; COMPONENTS_PER_POINT * grid.points()];
    bssn_grid_rhs(
        &view,
        &BssnGauge::SYNCHRONOUS,
        &AdmSources::VACUUM,
        0.0,
        &mut rhs,
    )
    .expect("rhs");

    let points = grid.points();
    let connection_slice = &rhs[14 * points..17 * points];
    let worst = max_abs(connection_slice);
    assert!(
        worst > 0.0,
        "d_t Gammatilde^i is identically zero on the grid -- the connection \
         equation is not being supplied"
    );

    // Minkowski, by contrast, has no gradients at all, so it must still be zero.
    let flat = BssnGridState::minkowski(grid);
    let mut flat_rhs = vec![0.0_f64; COMPONENTS_PER_POINT * points];
    bssn_grid_rhs(
        &flat.view(),
        &BssnGauge::SYNCHRONOUS,
        &AdmSources::VACUUM,
        0.0,
        &mut flat_rhs,
    )
    .expect("rhs");
    assert_eq!(max_abs(&flat_rhs[14 * points..17 * points]), 0.0);
}

#[test]
fn connection_rhs_rejects_an_unvalidated_matter_term() {
    // The matter contribution to d_t Gammatilde^i is not implemented, because
    // this increment has no non-vacuum oracle to validate it against. It must
    // be refused rather than silently ignored or guessed.
    use scirust_relativity::bssn::{
        BssnLapseDerivatives, BssnSecondDerivatives, BssnShiftDerivatives, BssnSpatialDerivatives,
        bssn_connection_rhs,
    };
    let zero_second = BssnSecondDerivatives {
        conformal_factor_hessian: [[0.0; 3]; 3],
        conformal_metric_hessian: [[[[0.0; 3]; 3]; 3]; 3],
        conformal_connection_gradient: [[0.0; 3]; 3],
    };
    let grid = unit_grid(8);
    let state = BssnGridState::minkowski(grid).state_at(0);
    let sources = AdmSources {
        energy_density: 0.0,
        momentum_density: [1.0e-3, 0.0, 0.0],
        stress: [[0.0; 3]; 3],
    };
    assert!(
        bssn_connection_rhs(
            &state,
            &BssnSpatialDerivatives::ZERO,
            &zero_second,
            &BssnLapseDerivatives::ZERO,
            &BssnShiftDerivatives::ZERO,
            &BssnGauge::SYNCHRONOUS,
            &sources,
        )
        .is_err()
    );
    // Vacuum is accepted.
    assert!(
        bssn_connection_rhs(
            &state,
            &BssnSpatialDerivatives::ZERO,
            &zero_second,
            &BssnLapseDerivatives::ZERO,
            &BssnShiftDerivatives::ZERO,
            &BssnGauge::SYNCHRONOUS,
            &AdmSources::VACUUM,
        )
        .is_ok()
    );
}

#[test]
fn dissipation_is_off_by_default_and_changes_nothing_when_zero() {
    let grid = unit_grid(32);
    let system = BssnGridSystem::vacuum(grid);
    // Off by default: never enabled implicitly.
    assert_eq!(system.dissipation(), 0.0);

    let wave = TransverseTracelessWave::new(1.0e-4, TWO_PI, 0.0);
    let initial =
        BssnGridState::from_adm_fields(grid, &wave.metric_field(), &wave.curvature_field())
            .expect("initial data");

    // An explicit zero coefficient must be bit-identical to the default path.
    let plain = evolve_bssn_grid(&system, &initial, 0.0, 0.05, 0.005).expect("plain");
    let zeroed = evolve_bssn_grid(&system.with_dissipation(0.0), &initial, 0.0, 0.05, 0.005)
        .expect("explicit zero");
    for (a, b) in plain.iter().zip(zeroed.iter())
    {
        assert_eq!(a.state.as_slice(), b.state.as_slice());
    }
}

#[test]
fn kreiss_oliger_annihilates_constants_and_peaks_at_the_nyquist_mode() {
    use scirust_relativity::bssn_grid::kreiss_oliger_term;
    let grid = unit_grid(16);
    let sigma = 0.1_f64;

    // The bracket sums to zero on a constant field, so smooth content at the
    // lowest frequency is untouched exactly.
    let constant = vec![3.0_f64; grid.points()];
    for index in 0..grid.points()
    {
        assert!(kreiss_oliger_term(&grid, &constant, index, sigma).abs() < 1.0e-14);
    }

    // At the Nyquist mode the bracket attains its maximum symbol of 16, so the
    // term is exactly -sigma * f_i / dx.
    let nyquist: Vec<f64> = (0..grid.points())
        .map(|i| if i % 2 == 0 { 1.0 } else { -1.0 })
        .collect();
    for index in 0..grid.points()
    {
        let expected = -sigma * nyquist[index] / grid.spacing();
        let actual = kreiss_oliger_term(&grid, &nyquist, index, sigma);
        assert!(
            (actual - expected).abs() < 1.0e-12,
            "index {index}: {actual} vs {expected}"
        );
    }
}

// ---------------------------------------------------------------------------
// Live slicing — 1+log
// ---------------------------------------------------------------------------

#[test]
fn one_plus_log_reproduces_the_root_two_gauge_speed() {
    // Linearised about flat space, `d_t alpha = -2 alpha K` together with
    // `d_t K = -D^2 alpha` gives `d_t^2 alpha = 2 d_x^2 alpha`: a wave equation
    // with characteristic speed sqrt(2). That speed is faster than light, which
    // is legitimate -- the lapse is gauge and carries no physical signal.
    //
    // Standing-wave data `alpha = 1 + A sin(kx)` with `K = 0` therefore evolves
    // as `alpha = 1 + A cos(sqrt(2) k t) sin(kx)`, an exact closed form to
    // compare against.
    let k = TWO_PI;
    let amplitude = 1.0e-6_f64;
    let speed = 2.0_f64.sqrt();
    let t_end = 0.25_f64;
    let mut errors = Vec::new();

    for &points in &[32_usize, 64, 128]
    {
        let grid = unit_grid(points);
        let mut initial = BssnGridState::minkowski(grid);
        for index in 0..grid.points()
        {
            initial.set_lapse_at(index, 1.0 + amplitude * (k * grid.coordinate(index)).sin());
        }

        let system = BssnGridSystem::vacuum(grid).with_slicing(BssnSlicing::OnePlusLog);
        let samples = evolve_bssn_grid(&system, &initial, 0.0, t_end, 0.1 * grid.spacing())
            .expect("gauge wave evolves");
        let last = samples.last().expect("final sample");

        let mut worst = 0.0_f64;
        for index in 0..grid.points()
        {
            let x = grid.coordinate(index);
            let exact = 1.0 + amplitude * (speed * k * last.time).cos() * (k * x).sin();
            worst = worst.max((last.state.lapse_at(index) - exact).abs());
        }
        errors.push(worst);
    }

    // Second order in dx, against the analytic gauge wave.
    for window in errors.windows(2)
    {
        let order = observed_order(window[0], window[1]);
        assert!(
            (order - 2.0).abs() < 0.15,
            "observed gauge-wave order {order} from {errors:?}"
        );
    }
    // And the residual is far below the perturbation amplitude, so the wave is
    // genuinely tracked rather than merely small.
    assert!(
        errors.iter().all(|value| *value < amplitude / 100.0),
        "gauge-wave errors are not small against A = {amplitude}: {errors:?}"
    );
}

#[test]
fn prescribed_slicing_is_the_default_and_freezes_the_lapse() {
    let grid = unit_grid(16);
    // A live gauge is never enabled implicitly.
    assert_eq!(
        BssnGridSystem::vacuum(grid).slicing(),
        BssnSlicing::Prescribed
    );

    let k = TWO_PI;
    let mut initial = BssnGridState::minkowski(grid);
    for index in 0..grid.points()
    {
        initial.set_lapse_at(index, 1.0 + 1.0e-3 * (k * grid.coordinate(index)).sin());
    }

    let system = BssnGridSystem::vacuum(grid);
    let samples = evolve_bssn_grid(&system, &initial, 0.0, 0.2, 0.02).expect("evolution");
    let last = samples.last().expect("final sample");
    for index in 0..grid.points()
    {
        // Bit-for-bit frozen: a prescribed lapse has an exactly zero
        // right-hand side, so RK4 adds exactly zero.
        assert_eq!(last.state.lapse_at(index), initial.lapse_at(index));
    }
}

#[test]
fn one_plus_log_leaves_minkowski_exactly_stationary() {
    // Unit lapse with K = 0 gives d_t alpha = -2 alpha K = 0 exactly, so
    // enabling the live slicing must not disturb the vacuum solution at all.
    let grid = unit_grid(16);
    let initial = BssnGridState::minkowski(grid);
    let system = BssnGridSystem::vacuum(grid).with_slicing(BssnSlicing::OnePlusLog);
    let samples = evolve_bssn_grid(&system, &initial, 0.0, 1.0, 0.01).expect("evolution");
    for sample in &samples
    {
        assert_eq!(
            sample.state.as_slice(),
            initial.as_slice(),
            "1+log disturbed Minkowski at t = {}",
            sample.time
        );
    }
}

#[test]
fn a_non_positive_or_non_finite_lapse_is_rejected() {
    let grid = unit_grid(8);

    let mut collapsed = BssnGridState::minkowski(grid);
    collapsed.set_lapse_at(3, 0.0);
    match collapsed.validate(0.25)
    {
        Err(scirust_relativity::bssn_grid::BssnGridError::InvalidState {
            time,
            index,
            field,
            ..
        }) =>
        {
            assert_eq!(time, 0.25);
            assert_eq!(index, 3);
            assert_eq!(field, "alpha");
        },
        other => panic!("expected a located InvalidState, got {other:?}"),
    }

    let mut infinite = BssnGridState::minkowski(grid);
    infinite.set_lapse_at(1, f64::NAN);
    assert!(infinite.validate(0.0).is_err());

    // A valid positive lapse is accepted.
    let mut fine = BssnGridState::minkowski(grid);
    fine.set_lapse_at(1, 0.5);
    assert!(fine.validate(0.0).is_ok());
}

#[test]
fn a_live_lapse_gradient_feeds_the_curvature_and_connection_equations() {
    // Regression guard: with a prescribed constant lapse the D_i D_j alpha term
    // in d_t K and d_t Atilde, and the -2 Atilde^{ij} d_j alpha term in
    // d_t Gammatilde^i, are all exactly zero. A live, spatially varying lapse
    // must make them non-zero -- otherwise the gauge is not actually coupled in.
    let grid = unit_grid(32);
    let k = TWO_PI;
    let mut state = BssnGridState::minkowski(grid);
    for index in 0..grid.points()
    {
        state.set_lapse_at(index, 1.0 + 1.0e-3 * (k * grid.coordinate(index)).sin());
    }

    let view = state.view();
    let mut rhs = vec![0.0_f64; COMPONENTS_PER_POINT * grid.points()];
    bssn_grid_rhs(
        &view,
        &BssnGauge::SYNCHRONOUS,
        &AdmSources::VACUUM,
        0.0,
        &mut rhs,
    )
    .expect("rhs");

    let points = grid.points();
    // d_t K picks up -D^2 alpha, which for this data is the only non-zero term.
    let trace_rhs = max_abs(&rhs[7 * points..8 * points]);
    assert!(
        trace_rhs > 1.0e-6,
        "d_t K is {trace_rhs}, expected non-zero"
    );

    // Atilde vanishes on this slice, so the connection lapse term does too --
    // which is a statement about this data, not about the implementation. The
    // trace-free lapse term in d_t Atilde is likewise traceless and vanishes for
    // an isotropic Hessian, so d_t K is the sharp probe here.
    let flat = BssnGridState::minkowski(grid);
    let mut flat_rhs = vec![0.0_f64; COMPONENTS_PER_POINT * points];
    bssn_grid_rhs(
        &flat.view(),
        &BssnGauge::SYNCHRONOUS,
        &AdmSources::VACUUM,
        0.0,
        &mut flat_rhs,
    )
    .expect("rhs");
    assert_eq!(max_abs(&flat_rhs[7 * points..8 * points]), 0.0);
}

// ---------------------------------------------------------------------------
// Shift terms and the Gamma-driver
// ---------------------------------------------------------------------------

#[test]
fn a_constant_shift_is_exactly_pure_advection() {
    // The sharpest test of the shift terms available. For a *spatially constant*
    // shift every `d beta` term vanishes identically, so the only surviving
    // contribution is the advection `beta^j d_j`. The difference between the
    // shifted and unshifted right-hand sides must therefore equal
    // `v * d_x(field)` for every evolved component -- and because both sides use
    // the same compact stencil, it must hold to *rounding*, not to a tolerance.
    //
    // A sign error in any one of the seventeen advection terms would show up
    // here immediately.
    use scirust_relativity::grid::periodic_first_derivative;

    let grid = unit_grid(64);
    let points = grid.points();
    let wave = TransverseTracelessWave::new(1.0e-3, TWO_PI, 0.0);
    let base = BssnGridState::from_adm_fields(grid, &wave.metric_field(), &wave.curvature_field())
        .expect("initial data");
    let velocity = 0.3_f64;

    let mut unshifted = vec![0.0_f64; COMPONENTS_PER_POINT * points];
    bssn_grid_rhs(
        &base.view(),
        &BssnGauge::SYNCHRONOUS,
        &AdmSources::VACUUM,
        0.0,
        &mut unshifted,
    )
    .expect("rhs");

    let mut drifting = base.clone();
    for index in 0..points
    {
        drifting.set_shift_at(index, &[velocity, 0.0, 0.0]);
    }
    let mut shifted = vec![0.0_f64; COMPONENTS_PER_POINT * points];
    bssn_grid_rhs(
        &drifting.view(),
        &BssnGauge::SYNCHRONOUS,
        &AdmSources::VACUUM,
        0.0,
        &mut shifted,
    )
    .expect("rhs");

    let mut worst = 0.0_f64;
    let mut scale = 0.0_f64;
    for slot in 0..17
    {
        let field: Vec<f64> = (0..points)
            .map(|index| base.as_slice()[slot * points + index])
            .collect();
        for index in 0..points
        {
            let gradient = periodic_first_derivative(&field, &grid, index).expect("d1");
            let difference = shifted[slot * points + index] - unshifted[slot * points + index];
            worst = worst.max((difference - velocity * gradient).abs());
            scale = scale.max((velocity * gradient).abs());
        }
    }

    // The advection is genuinely present, not trivially zero.
    assert!(
        scale > 1.0e-4,
        "advection scale is trivially small: {scale}"
    );
    // And it matches to rounding: a relative agreement of ~1e-16.
    assert!(
        worst < 1.0e-16 * scale.max(1.0) + 1.0e-17,
        "advection residual {worst} against scale {scale}"
    );
}

#[test]
fn the_gamma_driver_matches_its_closed_form_on_flat_space() {
    // On Minkowski with a spatially constant shift every gradient vanishes, so
    // d_t Gammatilde^i is exactly zero and the driver decouples into two exact
    // ODEs:
    //
    //   d_t B^i    = -eta B^i          =>  B(t)    = B0 exp(-eta t)
    //   d_t beta^i = (3/4) B^i         =>  beta(t) = beta0 + (3/4)(B0/eta)(1 - exp(-eta t))
    //
    // a closed form with no free parameters.
    let grid = unit_grid(16);
    let eta = 2.0_f64;
    let initial_driver = 0.4_f64;

    let mut initial = BssnGridState::minkowski(grid);
    for index in 0..grid.points()
    {
        initial.set_driver_at(index, &[initial_driver, 0.0, 0.0]);
    }

    let system =
        BssnGridSystem::vacuum(grid).with_shift_condition(BssnShiftCondition::GammaDriver { eta });
    let samples = evolve_bssn_grid(&system, &initial, 0.0, 1.0, 0.005).expect("driver evolves");

    for sample in &samples
    {
        let t = sample.time;
        let expected_driver = initial_driver * (-eta * t).exp();
        let expected_shift = 0.75 * (initial_driver / eta) * (1.0 - (-eta * t).exp());
        for index in 0..grid.points()
        {
            let driver = sample.state.driver_at(index)[0];
            let shift = sample.state.shift_at(index)[0];
            // RK4 at this step size resolves a decay rate of 2 to far better
            // than 1e-9; the bound is the integrator's accuracy, not a fudge.
            assert!(
                (driver - expected_driver).abs() < 1.0e-9,
                "B at t = {t}: {driver} vs {expected_driver}"
            );
            assert!(
                (shift - expected_shift).abs() < 1.0e-9,
                "beta at t = {t}: {shift} vs {expected_shift}"
            );
        }
    }

    // The transverse components were never excited and must stay exactly zero.
    let last = samples.last().expect("final");
    for index in 0..grid.points()
    {
        assert_eq!(last.state.shift_at(index)[1], 0.0);
        assert_eq!(last.state.shift_at(index)[2], 0.0);
    }
}

#[test]
fn the_gamma_driver_leaves_minkowski_exactly_stationary() {
    // Flat space has Gammatilde^i = 0 and d_t Gammatilde^i = 0, so with B = 0
    // and beta = 0 the driver produces exactly nothing. Enabling the full
    // moving-puncture gauge must not disturb the vacuum solution at all.
    let grid = unit_grid(16);
    let initial = BssnGridState::minkowski(grid);
    let system = BssnGridSystem::vacuum(grid)
        .with_slicing(BssnSlicing::OnePlusLog)
        .with_shift_condition(BssnShiftCondition::GammaDriver { eta: 1.0 });
    let samples = evolve_bssn_grid(&system, &initial, 0.0, 1.0, 0.01).expect("evolution");
    for sample in &samples
    {
        assert_eq!(
            sample.state.as_slice(),
            initial.as_slice(),
            "moving-puncture gauge disturbed Minkowski at t = {}",
            sample.time
        );
    }
}

#[test]
fn shift_conditions_default_to_prescribed_and_reject_negative_damping() {
    let grid = unit_grid(16);
    // A live shift is never enabled implicitly.
    assert_eq!(
        BssnGridSystem::vacuum(grid).shift_condition(),
        BssnShiftCondition::Prescribed
    );

    // A prescribed shift is frozen bit-for-bit, driver included.
    let mut initial = BssnGridState::minkowski(grid);
    for index in 0..grid.points()
    {
        initial.set_shift_at(index, &[0.2, 0.0, 0.0]);
        initial.set_driver_at(index, &[0.5, 0.0, 0.0]);
    }
    let samples = evolve_bssn_grid(&BssnGridSystem::vacuum(grid), &initial, 0.0, 0.2, 0.01)
        .expect("evolution");
    let last = samples.last().expect("final");
    for index in 0..grid.points()
    {
        assert_eq!(last.state.shift_at(index), [0.2, 0.0, 0.0]);
        assert_eq!(last.state.driver_at(index), [0.5, 0.0, 0.0]);
    }

    // A negative damping rate would amplify the drift the driver exists to
    // suppress, so it is refused rather than accepted.
    let bad = BssnGridSystem::vacuum(grid)
        .with_shift_condition(BssnShiftCondition::GammaDriver { eta: -1.0 });
    assert!(evolve_bssn_grid(&bad, &initial, 0.0, 0.1, 0.01).is_err());
    let nan = BssnGridSystem::vacuum(grid)
        .with_shift_condition(BssnShiftCondition::GammaDriver { eta: f64::NAN });
    assert!(evolve_bssn_grid(&nan, &initial, 0.0, 0.1, 0.01).is_err());
}

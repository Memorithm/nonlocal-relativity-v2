use scirust_nonlocal_relativity::{
    AdaptiveNonlocalConfig, CaputoCoordinateMemory, CompleteUniformHistory, HeunPeceStepper,
    IdentityHistoryTransport, NonlocalConfig, NonlocalRelativityError, NonlocalSimulationPolicy,
    WorldlineState, simulate_nonlocal_worldline_adaptive, simulate_nonlocal_worldline_with_policy,
};
use scirust_relativity::{Minkowski, Schwarzschild};
use std::f64::consts::FRAC_PI_2;

fn assert_close(actual: f64, expected: f64, tolerance: f64) {
    let scale = expected.abs().max(1.0);
    let relative_error = (actual - expected).abs() / scale;

    assert!(
        relative_error <= tolerance,
        "actual={actual:.17e}, expected={expected:.17e}, \
         relative_error={relative_error:.17e}, tolerance={tolerance:.17e}"
    );
}

fn circular_schwarzschild_state(mass: f64, radius: f64) -> WorldlineState<4> {
    let denominator = (1.0 - 3.0 * mass / radius).sqrt();
    let u_t = 1.0 / denominator;
    let u_phi = (mass / (radius * radius * radius)).sqrt() / denominator;

    WorldlineState::new([0.0, radius, FRAC_PI_2, 0.0], [u_t, 0.0, 0.0, u_phi])
}

#[test]
fn constructor_rejects_invalid_parameters() {
    assert!(matches!(
        AdaptiveNonlocalConfig::new(1.5, 0.02, 0.05, 0.001, 0.5, 1.0e-9, 1.0e-8, 2.0, 500, 20),
        Err(NonlocalRelativityError::InvalidFractionalOrder(_))
    ));
    assert!(matches!(
        AdaptiveNonlocalConfig::new(0.5, -0.1, 0.05, 0.001, 0.5, 1.0e-9, 1.0e-8, 2.0, 500, 20),
        Err(NonlocalRelativityError::InvalidCoupling(_))
    ));
    assert!(matches!(
        AdaptiveNonlocalConfig::new(0.5, 0.02, 0.05, 0.0, 0.5, 1.0e-9, 1.0e-8, 2.0, 500, 20),
        Err(NonlocalRelativityError::InvalidAdaptiveConfiguration {
            field: "min_step",
            ..
        })
    ));
    assert!(matches!(
        AdaptiveNonlocalConfig::new(0.5, 0.02, 0.05, 0.1, 0.05, 1.0e-9, 1.0e-8, 2.0, 500, 20),
        Err(NonlocalRelativityError::InvalidAdaptiveConfiguration {
            field: "max_step",
            ..
        })
    ));
    assert!(matches!(
        AdaptiveNonlocalConfig::new(0.5, 0.02, 1.0, 0.001, 0.5, 1.0e-9, 1.0e-8, 2.0, 500, 20),
        Err(NonlocalRelativityError::InvalidAdaptiveConfiguration {
            field: "initial_step",
            ..
        })
    ));
    assert!(matches!(
        AdaptiveNonlocalConfig::new(0.5, 0.02, 0.05, 0.001, 0.5, 0.0, 1.0e-8, 2.0, 500, 20),
        Err(NonlocalRelativityError::InvalidAdaptiveConfiguration {
            field: "error_tolerance",
            ..
        })
    ));
    assert!(matches!(
        AdaptiveNonlocalConfig::new(0.5, 0.02, 0.05, 0.001, 0.5, 1.0e-9, 0.0, 2.0, 500, 20),
        Err(NonlocalRelativityError::InvalidMetricNormFloor(_))
    ));
    assert!(matches!(
        AdaptiveNonlocalConfig::new(0.5, 0.02, 0.05, 0.001, 0.5, 1.0e-9, 1.0e-8, 0.0, 500, 20),
        Err(NonlocalRelativityError::InvalidAdaptiveConfiguration {
            field: "target_affine_parameter",
            ..
        })
    ));
    assert!(matches!(
        AdaptiveNonlocalConfig::new(0.5, 0.02, 0.05, 0.001, 0.5, 1.0e-9, 1.0e-8, 2.0, 0, 20),
        Err(NonlocalRelativityError::InvalidStepCount(0))
    ));
    assert!(matches!(
        AdaptiveNonlocalConfig::new(0.5, 0.02, 0.05, 0.001, 0.5, 1.0e-9, 1.0e-8, 2.0, 500, 0),
        Err(NonlocalRelativityError::InvalidAdaptiveConfiguration {
            field: "max_rejections_per_step",
            ..
        })
    ));
    assert!(
        AdaptiveNonlocalConfig::new(0.5, 0.02, 0.05, 0.001, 0.5, 1.0e-9, 1.0e-8, 2.0, 500, 20)
            .is_ok()
    );
}

#[test]
fn adaptive_simulation_reaches_the_target_affine_parameter() {
    let mass = 1.0;
    let background = Schwarzschild::try_new(mass).unwrap();
    let mut initial = circular_schwarzschild_state(mass, 10.0);
    initial.velocity[1] = -0.01;
    let target = 1.6;
    let config = AdaptiveNonlocalConfig::new(
        0.55, 0.02, 0.02, 0.001, 0.2, 1.0e-9, 1.0e-8, target, 5_000, 30,
    )
    .unwrap();

    let trajectory = simulate_nonlocal_worldline_adaptive(&background, initial, config).unwrap();

    let final_parameter = trajectory.final_diagnostics().unwrap().affine_parameter;
    assert_close(final_parameter, target, 1.0e-9);
}

#[test]
fn adaptive_trajectory_uses_non_uniform_steps() {
    let mass = 1.0;
    let background = Schwarzschild::try_new(mass).unwrap();
    let mut initial = circular_schwarzschild_state(mass, 10.0);
    initial.velocity[1] = -0.01;
    // A deliberately generous initial step relative to the tolerance forces
    // early rejections and shrinkage.
    let config = AdaptiveNonlocalConfig::new(
        0.55, 0.02, 0.15, 0.00005, 0.2, 1.0e-8, 1.0e-8, 1.0, 5_000, 60,
    )
    .unwrap();

    let trajectory = simulate_nonlocal_worldline_adaptive(&background, initial, config).unwrap();
    let parameters: Vec<f64> = trajectory
        .diagnostics()
        .iter()
        .map(|d| d.affine_parameter)
        .collect();

    assert!(parameters.len() > 2, "expected more than two samples");

    let mut increments: Vec<f64> = parameters.windows(2).map(|w| w[1] - w[0]).collect();
    increments.dedup_by(|a, b| (*a - *b).abs() < 1.0e-15);

    assert!(
        increments.len() > 1,
        "expected non-uniform step sizes, got a single repeated increment"
    );
}

#[test]
fn adaptive_is_deterministic_bit_for_bit() {
    let mass = 1.0;
    let background = Schwarzschild::try_new(mass).unwrap();
    let mut initial = circular_schwarzschild_state(mass, 10.0);
    initial.velocity[1] = -0.01;
    let config =
        AdaptiveNonlocalConfig::new(0.55, 0.02, 0.05, 0.001, 0.2, 1.0e-9, 1.0e-8, 0.8, 5_000, 30)
            .unwrap();

    let first = simulate_nonlocal_worldline_adaptive(&background, initial, config).unwrap();
    let second = simulate_nonlocal_worldline_adaptive(&background, initial, config).unwrap();

    assert_eq!(first.len(), second.len());
    for (left, right) in first.states().iter().zip(second.states())
    {
        for component in 0..4
        {
            assert_eq!(
                left.coordinates[component].to_bits(),
                right.coordinates[component].to_bits()
            );
            assert_eq!(
                left.velocity[component].to_bits(),
                right.velocity[component].to_bits()
            );
        }
    }
}

#[test]
fn constant_velocity_in_flat_spacetime_has_exactly_zero_memory() {
    let initial = WorldlineState::new([1.0, -2.0, 3.0, -4.0], [2.0, 0.25, -0.5, 0.75]);
    let config =
        AdaptiveNonlocalConfig::new(0.5, 0.03, 0.05, 0.001, 0.2, 1.0e-9, 1.0e-8, 1.0, 5_000, 30)
            .unwrap();

    let trajectory = simulate_nonlocal_worldline_adaptive(&Minkowski, initial, config).unwrap();

    for state in trajectory.states()
    {
        for component in 0..4
        {
            assert_eq!(
                state.velocity[component].to_bits(),
                initial.velocity[component].to_bits()
            );
        }
    }

    for diagnostics in trajectory.diagnostics()
    {
        assert_eq!(diagnostics.memory_l2_norm.to_bits(), 0.0_f64.to_bits());
    }
}

#[test]
fn step_budget_exhaustion_is_reported_not_silently_truncated() {
    let mass = 1.0;
    let background = Schwarzschild::try_new(mass).unwrap();
    let mut initial = circular_schwarzschild_state(mass, 10.0);
    initial.velocity[1] = -0.01;
    // Deliberately too few accepted steps to reach the target.
    let config =
        AdaptiveNonlocalConfig::new(0.55, 0.02, 0.02, 0.001, 0.02, 1.0e-9, 1.0e-8, 5.0, 3, 30)
            .unwrap();

    let result = simulate_nonlocal_worldline_adaptive(&background, initial, config);

    assert!(matches!(
        result,
        Err(NonlocalRelativityError::AdaptiveStepBudgetExhausted {
            accepted_steps: 3,
            ..
        })
    ));
}

#[test]
fn rejects_non_finite_initial_state() {
    let background = Minkowski;
    let initial = WorldlineState::new([f64::NAN, 0.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0]);
    let config =
        AdaptiveNonlocalConfig::new(0.5, 0.02, 0.05, 0.001, 0.2, 1.0e-9, 1.0e-8, 1.0, 500, 20)
            .unwrap();

    let result = simulate_nonlocal_worldline_adaptive(&background, initial, config);

    assert!(matches!(
        result,
        Err(NonlocalRelativityError::NonFiniteInitialCoordinate { .. })
    ));
}

#[test]
fn adaptive_matches_fine_fixed_step_heun_pece_closely() {
    let mass = 1.0;
    let background = Schwarzschild::try_new(mass).unwrap();
    let mut initial = circular_schwarzschild_state(mass, 10.0);
    initial.velocity[1] = -0.01;
    let alpha = 0.55;
    let coupling = 0.02;
    let target = 0.8;

    let adaptive_config = AdaptiveNonlocalConfig::new(
        alpha, coupling, 0.01, 0.00002, 0.05, 1.0e-8, 1.0e-8, target, 20_000, 60,
    )
    .unwrap();
    let adaptive_trajectory =
        simulate_nonlocal_worldline_adaptive(&background, initial, adaptive_config).unwrap();

    let fine_step = 0.005;
    let fine_steps = (target / fine_step).round() as usize;
    let fixed_config = NonlocalConfig::new(alpha, coupling, fine_step, fine_steps, 1.0e-8).unwrap();
    let fixed_trajectory = simulate_nonlocal_worldline_with_policy(
        &background,
        initial,
        fixed_config,
        NonlocalSimulationPolicy::new(
            CompleteUniformHistory::<4>::with_capacity(fine_steps + 1),
            CaputoCoordinateMemory,
            IdentityHistoryTransport,
            HeunPeceStepper,
        ),
    )
    .unwrap();

    let adaptive_final = adaptive_trajectory.final_state().unwrap();
    let fixed_final = fixed_trajectory.final_state().unwrap();

    for component in 0..4
    {
        assert_close(
            adaptive_final.coordinates[component],
            fixed_final.coordinates[component],
            1.0e-4,
        );
        assert_close(
            adaptive_final.velocity[component],
            fixed_final.velocity[component],
            1.0e-4,
        );
    }
}

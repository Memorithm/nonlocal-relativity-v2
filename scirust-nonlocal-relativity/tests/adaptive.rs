use scirust_nonlocal_relativity::{
    AdaptiveNonlocalConfig, AdaptiveSimulationPolicy, BoundedShortMemoryHistory,
    CaputoCoordinateMemory, CompleteUniformHistory, DiscreteConnectionTransport, HeunPeceStepper,
    IdentityHistoryModulator, IdentityHistoryTransport, NonlocalConfig, NonlocalRelativityError,
    NonlocalSimulationPolicy, ReissnerNordstromFieldModulator, SchwarzschildKretschmannModulator,
    WorldlineState, simulate_nonlocal_worldline_adaptive,
    simulate_nonlocal_worldline_adaptive_with_policy, simulate_nonlocal_worldline_with_policy,
};
use scirust_relativity::{Minkowski, ReissnerNordstrom, Schwarzschild};
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

/// Golden bit-for-bit regression captured from the pre-refactor
/// implementation of `simulate_nonlocal_worldline_adaptive` (a hand-rolled
/// history loop with no `HistoryBackend`/`HistoryTransport`/
/// `HistoryModulator` involvement), before it was reimplemented as a thin
/// wrapper over `simulate_nonlocal_worldline_adaptive_with_policy` with
/// `CompleteUniformHistory`, `IdentityHistoryTransport`, and
/// `IdentityHistoryModulator`. Composing these identity components must
/// reproduce the original numbers exactly, not merely approximately.
#[test]
fn adaptive_reproduces_pre_refactor_golden_values_bit_for_bit() {
    let mass = 1.0;
    let background = Schwarzschild::try_new(mass).unwrap();
    let mut initial = circular_schwarzschild_state(mass, 10.0);
    initial.velocity[1] = -0.01;
    let config =
        AdaptiveNonlocalConfig::new(0.55, 0.02, 0.05, 0.001, 0.2, 1.0e-9, 1.0e-8, 0.8, 5_000, 30)
            .unwrap();

    let trajectory = simulate_nonlocal_worldline_adaptive(&background, initial, config).unwrap();

    assert_eq!(trajectory.len(), 353);

    let final_state = trajectory.final_state().unwrap();
    let expected_coordinates_bits: [u64; 4] = [
        0x3fee99d338548003,
        0x4023fbe788e8c53d,
        0x3ff921fb54442d18,
        0x3f9efcc35f493b88,
    ];
    let expected_velocity_bits: [u64; 4] = [
        0x3ff3209f6d5766ff,
        0xbf84797d45a810cd,
        0x3bf46a6831b53a3a,
        0x3fa361e1756df532,
    ];
    for component in 0..4
    {
        assert_eq!(
            final_state.coordinates[component].to_bits(),
            expected_coordinates_bits[component],
            "coordinate {component}"
        );
        assert_eq!(
            final_state.velocity[component].to_bits(),
            expected_velocity_bits[component],
            "velocity {component}"
        );
    }

    let final_diagnostics = trajectory.final_diagnostics().unwrap();
    assert_eq!(
        final_diagnostics.affine_parameter.to_bits(),
        0x3fe999999999999a
    );
    assert_eq!(
        final_diagnostics.memory_l2_norm.to_bits(),
        0x3f345109f7290677
    );
}

#[test]
fn adaptive_with_identity_policy_matches_plain_entry_point_bit_for_bit() {
    let mass = 1.0;
    let background = Schwarzschild::try_new(mass).unwrap();
    let mut initial = circular_schwarzschild_state(mass, 10.0);
    initial.velocity[1] = -0.01;
    let config =
        AdaptiveNonlocalConfig::new(0.55, 0.02, 0.05, 0.001, 0.2, 1.0e-9, 1.0e-8, 0.8, 5_000, 30)
            .unwrap();

    let plain = simulate_nonlocal_worldline_adaptive(&background, initial, config).unwrap();
    let explicit = simulate_nonlocal_worldline_adaptive_with_policy(
        &background,
        initial,
        config,
        AdaptiveSimulationPolicy::new(
            CompleteUniformHistory::<4>::new(),
            IdentityHistoryTransport,
            IdentityHistoryModulator,
        ),
    )
    .unwrap();

    assert_eq!(plain.len(), explicit.len());
    for (left, right) in plain.states().iter().zip(explicit.states())
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
fn adaptive_composes_with_discrete_connection_transport() {
    let mass = 1.0;
    let background = Schwarzschild::try_new(mass).unwrap();
    let mut initial = circular_schwarzschild_state(mass, 10.0);
    initial.velocity[1] = -0.01;
    let config =
        AdaptiveNonlocalConfig::new(0.55, 0.02, 0.05, 0.001, 0.2, 1.0e-9, 1.0e-8, 0.6, 5_000, 30)
            .unwrap();

    let transported = simulate_nonlocal_worldline_adaptive_with_policy(
        &background,
        initial,
        config,
        AdaptiveSimulationPolicy::new(
            CompleteUniformHistory::<4>::new(),
            DiscreteConnectionTransport,
            IdentityHistoryModulator,
        ),
    )
    .unwrap();

    for state in transported.states()
    {
        assert!(state.coordinates.iter().all(|value| value.is_finite()));
        assert!(state.velocity.iter().all(|value| value.is_finite()));
    }

    let identity = simulate_nonlocal_worldline_adaptive_with_policy(
        &background,
        initial,
        config,
        AdaptiveSimulationPolicy::new(
            CompleteUniformHistory::<4>::new(),
            IdentityHistoryTransport,
            IdentityHistoryModulator,
        ),
    )
    .unwrap();

    // DiscreteConnectionTransport must actually change the numerical result
    // here (Schwarzschild is curved), not silently fall back to identity
    // behavior.
    let transported_final = transported.final_state().unwrap();
    let identity_final = identity.final_state().unwrap();
    let differs = (0..4).any(|component| {
        transported_final.coordinates[component] != identity_final.coordinates[component]
    });
    assert!(differs, "transported and identity results were identical");
}

#[test]
fn adaptive_composes_with_schwarzschild_kretschmann_modulator() {
    let mass = 1.0;
    let background = Schwarzschild::try_new(mass).unwrap();
    let mut initial = circular_schwarzschild_state(mass, 10.0);
    initial.velocity[1] = -0.01;
    let config =
        AdaptiveNonlocalConfig::new(0.55, 0.02, 0.05, 0.001, 0.2, 1.0e-9, 1.0e-8, 0.6, 5_000, 30)
            .unwrap();

    let baseline_modulator = SchwarzschildKretschmannModulator::try_new(mass, 1.0, 0.0).unwrap();
    let coupled_modulator = SchwarzschildKretschmannModulator::try_new(mass, 1.0, 0.5).unwrap();

    let baseline = simulate_nonlocal_worldline_adaptive_with_policy(
        &background,
        initial,
        config,
        AdaptiveSimulationPolicy::new(
            CompleteUniformHistory::<4>::new(),
            IdentityHistoryTransport,
            baseline_modulator,
        ),
    )
    .unwrap();
    let coupled = simulate_nonlocal_worldline_adaptive_with_policy(
        &background,
        initial,
        config,
        AdaptiveSimulationPolicy::new(
            CompleteUniformHistory::<4>::new(),
            IdentityHistoryTransport,
            coupled_modulator,
        ),
    )
    .unwrap();

    for state in coupled.states()
    {
        assert!(state.coordinates.iter().all(|value| value.is_finite()));
    }

    let baseline_final = baseline.final_state().unwrap();
    let coupled_final = coupled.final_state().unwrap();
    let differs = (0..4).any(|component| {
        baseline_final.coordinates[component] != coupled_final.coordinates[component]
    });
    assert!(differs, "modulated and baseline results were identical");
}

#[test]
fn adaptive_beta_zero_modulator_matches_identity_modulator_bit_for_bit() {
    let mass = 1.0;
    let background = Schwarzschild::try_new(mass).unwrap();
    let mut initial = circular_schwarzschild_state(mass, 10.0);
    initial.velocity[1] = -0.01;
    let config =
        AdaptiveNonlocalConfig::new(0.55, 0.02, 0.05, 0.001, 0.2, 1.0e-9, 1.0e-8, 0.6, 5_000, 30)
            .unwrap();
    let beta_zero_modulator = SchwarzschildKretschmannModulator::try_new(mass, 1.0, 0.0).unwrap();

    let with_beta_zero = simulate_nonlocal_worldline_adaptive_with_policy(
        &background,
        initial,
        config,
        AdaptiveSimulationPolicy::new(
            CompleteUniformHistory::<4>::new(),
            IdentityHistoryTransport,
            beta_zero_modulator,
        ),
    )
    .unwrap();
    let with_identity = simulate_nonlocal_worldline_adaptive_with_policy(
        &background,
        initial,
        config,
        AdaptiveSimulationPolicy::new(
            CompleteUniformHistory::<4>::new(),
            IdentityHistoryTransport,
            IdentityHistoryModulator,
        ),
    )
    .unwrap();

    assert_eq!(with_beta_zero.len(), with_identity.len());
    for (left, right) in with_beta_zero.states().iter().zip(with_identity.states())
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
fn adaptive_composes_with_reissner_nordstrom_field_modulator() {
    let mass = 1.0;
    let charge = 0.4;
    let background = ReissnerNordstrom::try_new(mass, charge).unwrap();
    let mut initial = circular_schwarzschild_state(mass, 10.0);
    initial.velocity[1] = -0.01;
    let config =
        AdaptiveNonlocalConfig::new(0.55, 0.02, 0.05, 0.001, 0.2, 1.0e-9, 1.0e-8, 0.6, 5_000, 30)
            .unwrap();

    let modulator = ReissnerNordstromFieldModulator::try_new(background, 1.0, 0.5).unwrap();

    let modulated = simulate_nonlocal_worldline_adaptive_with_policy(
        &background,
        initial,
        config,
        AdaptiveSimulationPolicy::new(
            CompleteUniformHistory::<4>::new(),
            IdentityHistoryTransport,
            modulator,
        ),
    )
    .unwrap();

    for state in modulated.states()
    {
        assert!(state.coordinates.iter().all(|value| value.is_finite()));
        assert!(state.velocity.iter().all(|value| value.is_finite()));
    }
}

#[test]
fn adaptive_composes_with_both_transport_and_modulation_together() {
    let mass = 1.0;
    let background = Schwarzschild::try_new(mass).unwrap();
    let mut initial = circular_schwarzschild_state(mass, 10.0);
    initial.velocity[1] = -0.01;
    let config =
        AdaptiveNonlocalConfig::new(0.55, 0.02, 0.05, 0.001, 0.2, 1.0e-9, 1.0e-8, 0.6, 5_000, 30)
            .unwrap();
    let modulator = SchwarzschildKretschmannModulator::try_new(mass, 1.0, 0.3).unwrap();

    let trajectory = simulate_nonlocal_worldline_adaptive_with_policy(
        &background,
        initial,
        config,
        AdaptiveSimulationPolicy::new(
            CompleteUniformHistory::<4>::new(),
            DiscreteConnectionTransport,
            modulator,
        ),
    )
    .unwrap();

    for state in trajectory.states()
    {
        assert!(state.coordinates.iter().all(|value| value.is_finite()));
        assert!(state.velocity.iter().all(|value| value.is_finite()));
    }
    for diagnostics in trajectory.diagnostics()
    {
        assert!(diagnostics.memory_l2_norm.is_finite());
        assert!(diagnostics.memory_force_l2_norm.is_finite());
    }
}

#[test]
fn adaptive_composes_with_bounded_short_memory() {
    let mass = 1.0;
    let background = Schwarzschild::try_new(mass).unwrap();
    let mut initial = circular_schwarzschild_state(mass, 10.0);
    initial.velocity[1] = -0.01;
    let config =
        AdaptiveNonlocalConfig::new(0.55, 0.02, 0.05, 0.001, 0.2, 1.0e-9, 1.0e-8, 0.4, 5_000, 30)
            .unwrap();

    let trajectory = simulate_nonlocal_worldline_adaptive_with_policy(
        &background,
        initial,
        config,
        AdaptiveSimulationPolicy::new(
            BoundedShortMemoryHistory::<4>::new(6).unwrap(),
            IdentityHistoryTransport,
            IdentityHistoryModulator,
        ),
    )
    .unwrap();

    for state in trajectory.states()
    {
        assert!(state.coordinates.iter().all(|value| value.is_finite()));
    }
    assert!(
        trajectory
            .history_diagnostics()
            .iter()
            .all(|diagnostics| diagnostics.retained_samples <= 6),
        "bounded short memory retained more than its window"
    );
}

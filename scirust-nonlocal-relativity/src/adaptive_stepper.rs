//! Adaptive-step worldline integration that genuinely composes with
//! [`crate::MemoryLaw`] and [`crate::WorldlineStepper`]
//! ([`crate::SemiImplicitEulerStepper`] specifically), unlike
//! [`crate::simulate_nonlocal_worldline_adaptive_with_policy`]'s embedded
//! Heun-Euler controller, which cannot reuse either trait: both thread a
//! single fixed [`crate::NonlocalConfig`] step through their signatures
//! ([`crate::StepperContext`] for the latter), and
//! [`crate::CaputoCoordinateMemory`] specifically applies that one step to
//! the *entire* retained history via `caputo_l1_uniform`, which is only
//! correct when every accepted segment used the same step.
//!
//! This module closes that gap for [`crate::SemiImplicitEulerStepper`],
//! composed with the two non-uniform-aware memory laws in
//! [`crate::nonuniform_memory`] ([`crate::NonuniformCaputoCoordinateMemory`],
//! [`crate::NonuniformModulatedCaputoCoordinateMemory`]) that read each
//! retained sample's own recorded parameter instead of assuming uniform
//! spacing, via classical step-doubling rather than an embedded pair:
//!
//! 1. advance one full trial step of size `h` with
//!    [`crate::SemiImplicitEulerStepper`] from the accepted state (no memory
//!    evaluation needed for this branch: the accepted acceleration is
//!    already known from the previous accepted step);
//! 2. advance two steps of size `h/2`, evaluating the memory law at the
//!    midpoint against a throwaway provisional history clone (exactly like
//!    [`crate::HeunPeceStepper`]'s predictor push);
//! 3. the difference between the one-step and two-half-step results
//!    estimates the local error; the two-half-steps result (more accurate)
//!    is kept when it is accepted.
//!
//! [`crate::SemiImplicitEulerStepper`]'s `advance` is safe to reuse
//! unmodified here because its body only reads `context.state`,
//! `context.accepted_acceleration`, and `context.config.step` (the current
//! trial step size): it never reconstructs an absolute affine parameter from
//! `context.step_index`, unlike [`crate::HeunPeceStepper`].
//!
//! **[`crate::HeunPeceStepper`] is deliberately excluded and cannot be added
//! here without changing its existing body.** Its predictor pushes a
//! provisional history entry whose parameter it computes as
//! `(context.step_index + 1) as f64 * context.config.step`: an absolute
//! affine parameter reconstructed by multiplying an integer step count by
//! *one* step size. That is exact for the fixed-step architecture, where
//! every accepted step shares that size, but wrong the moment accepted step
//! sizes vary, which adaptive stepping does by construction — there is no
//! single `step` value for which `step_index * step` equals the true
//! accumulated parameter along a non-uniform trajectory. Correcting this
//! would require [`crate::StepperContext`] to carry the true accumulated
//! parameter directly instead of deriving it from `step_index`, which means
//! changing that existing, already-tested struct and
//! [`crate::HeunPeceStepper`]'s `advance` body itself — out of scope for an
//! additive change. Classical step-doubling with a first-order method (used
//! here) needs no predictor push in the first place, so it sidesteps the
//! problem rather than solving it for Heun-PECE: it is a different, standard
//! way to error-control a method with no natural embedded higher-order
//! partner, not a higher-accuracy replacement for the embedded pair.
//!
//! The step-size controller's growth/shrink exponent is `1 / (p_low + 1) =
//! 0.5`, exactly as in
//! [`crate::simulate_nonlocal_worldline_adaptive_with_policy`]'s embedded
//! pair, because semi-implicit Euler is *also* a first-order method
//! (`p_low = 1`, local truncation error `O(h^2)`): for a first-order method,
//! classical step-doubling's Richardson error estimate needs no rescaling
//! (the one-step/two-half-step difference divided by `2^1 - 1 = 1` is itself
//! already the estimate), so the raw difference computed below is used
//! directly, exactly as the embedded pair uses its raw predictor/corrector
//! difference directly.
//!
//! One further, independent difference from
//! [`crate::simulate_nonlocal_worldline_adaptive_with_policy`]:
//! [`crate::AdaptiveNonlocalConfig::max_rejections_per_step`] is validated
//! there but never consulted by that loop, which instead stops shrinking
//! only once the trial step falls below
//! [`crate::AdaptiveNonlocalConfig::min_step`]. This loop counts rejections
//! explicitly and returns
//! [`crate::NonlocalRelativityError::AdaptiveRejectionBudgetExhausted`] as
//! soon as either bound is hit, matching that field's documented meaning
//! ("the maximum number of consecutive rejections permitted") precisely.

use crate::{
    AdaptiveNonlocalConfig, CompleteUniformHistory, Connection, HistoryBackend, HistoryEntry,
    HistoryTransport, IdentityHistoryTransport, MemoryLaw, Metric, NonlocalConfig,
    NonlocalRelativityError, NonlocalResult, NonlocalTrajectory, NonuniformCaputoCoordinateMemory,
    SemiImplicitEulerStepper, StepEvaluationInput, StepperContext, WorldlineState,
    WorldlineStepper, evaluate_step_with_policy, scaled_local_error_norm, validate_initial_state,
    validated_metric, validated_metric_norm,
};

const STEP_SAFETY_FACTOR: f64 = 0.9;
const STEP_GROWTH_CAP: f64 = 4.0;
const STEP_SHRINK_FLOOR: f64 = 0.1;
/// Step-control exponent `1 / (p_low + 1)` for classical step-doubling with
/// semi-implicit Euler's order `p_low = 1`. See the module documentation.
const STEP_DOUBLING_ERROR_EXPONENT: f64 = 0.5;

/// Typed architecture policy for the step-doubling adaptive integrator:
/// which history backend, [`crate::MemoryLaw`], and
/// [`crate::HistoryTransport`] compose the non-uniform memory evaluation
/// around [`crate::SemiImplicitEulerStepper`].
///
/// This mirrors [`crate::NonlocalSimulationPolicy`]'s role for the
/// fixed-step architecture, narrowed to the components this module actually
/// varies: there is no stepper type parameter, because
/// [`crate::SemiImplicitEulerStepper`] is the only
/// [`crate::WorldlineStepper`] this module can soundly drive (see the module
/// documentation for why).
#[derive(Debug, Clone, PartialEq)]
pub struct AdaptiveStepperPolicy<H, L, T> {
    history_backend: H,
    memory_law: L,
    transport: T,
}

impl<H, L, T> AdaptiveStepperPolicy<H, L, T> {
    /// Construct a policy from explicit architecture components.
    #[must_use]
    pub const fn new(history_backend: H, memory_law: L, transport: T) -> Self {
        Self {
            history_backend,
            memory_law,
            transport,
        }
    }

    /// Borrow the policy history backend.
    #[must_use]
    pub const fn history_backend(&self) -> &H {
        &self.history_backend
    }

    /// Borrow the policy memory law.
    #[must_use]
    pub const fn memory_law(&self) -> &L {
        &self.memory_law
    }

    /// Borrow the policy history transport.
    #[must_use]
    pub const fn transport(&self) -> &T {
        &self.transport
    }
}

struct StepDoublingResult<const D: usize> {
    full_state: WorldlineState<D>,
    refined_state: WorldlineState<D>,
}

/// Attempt one step-doubling trial of size `step` from `state`: one full
/// [`crate::SemiImplicitEulerStepper`] step, compared against two half steps
/// with a memory-law evaluation at the midpoint (against a throwaway
/// provisional history clone). See the module documentation.
#[allow(clippy::too_many_arguments)]
fn step_doubling_trial<B, H, L, T, const D: usize>(
    background: &B,
    state: &WorldlineState<D>,
    accepted_acceleration: &[f64; D],
    history: &H,
    memory_law: &L,
    transport: &T,
    current_parameter: f64,
    step: f64,
    initial_metric_norm: f64,
    config_template: NonlocalConfig,
    step_index: usize,
) -> NonlocalResult<StepDoublingResult<D>>
where
    B: Metric<D> + Connection<D>,
    H: HistoryBackend<D>,
    L: MemoryLaw<D>,
    T: HistoryTransport<D>,
{
    let full_step_config = NonlocalConfig::from_fractional_order(
        config_template.fractional_order(),
        config_template.coupling(),
        step,
        1,
        config_template.metric_norm_floor(),
    )?;

    let full_state = SemiImplicitEulerStepper.advance(StepperContext {
        background,
        state,
        accepted_acceleration,
        history,
        memory_law,
        transport,
        initial_metric_norm,
        step_index,
        config: full_step_config,
    })?;

    let half_step = 0.5 * step;
    let half_step_config = NonlocalConfig::from_fractional_order(
        config_template.fractional_order(),
        config_template.coupling(),
        half_step,
        1,
        config_template.metric_norm_floor(),
    )?;

    let half_state_1 = SemiImplicitEulerStepper.advance(StepperContext {
        background,
        state,
        accepted_acceleration,
        history,
        memory_law,
        transport,
        initial_metric_norm,
        step_index,
        config: half_step_config,
    })?;

    let midpoint_parameter = current_parameter + half_step;
    let mut provisional_history = history.clone();
    provisional_history.push_entry(
        background,
        transport,
        HistoryEntry::new(
            half_state_1.coordinates,
            half_state_1.velocity,
            midpoint_parameter,
        ),
    )?;

    let midpoint_evaluation = evaluate_step_with_policy(StepEvaluationInput {
        background,
        state: &half_state_1,
        history: &provisional_history,
        memory_law,
        transport,
        initial_metric_norm,
        affine_parameter: midpoint_parameter,
        step_index,
        config: half_step_config,
    })?;

    let half_state_2 = SemiImplicitEulerStepper.advance(StepperContext {
        background,
        state: &half_state_1,
        accepted_acceleration: &midpoint_evaluation.acceleration,
        history: &provisional_history,
        memory_law,
        transport,
        initial_metric_norm,
        step_index,
        config: half_step_config,
    })?;

    Ok(StepDoublingResult {
        full_state,
        refined_state: half_state_2,
    })
}

/// Simulate the experimental fractional-memory worldline model with an
/// adaptive affine-parameter step, genuinely composed with
/// [`crate::MemoryLaw`] and [`crate::SemiImplicitEulerStepper`] via
/// classical step-doubling. See the module documentation for the mechanism,
/// why [`crate::SemiImplicitEulerStepper`] specifically is reused, and why
/// [`crate::HeunPeceStepper`] is not offered here.
///
/// As with [`crate::simulate_nonlocal_worldline_adaptive_with_policy`], the
/// returned [`crate::NonlocalTrajectory`] samples a generally non-uniform
/// affine-parameter axis, integration stops once the accumulated affine
/// parameter reaches `config.target_affine_parameter()` (the final accepted
/// step is clamped so the trajectory does not overshoot it), and neither an
/// exhausted step budget
/// ([`crate::NonlocalRelativityError::AdaptiveStepBudgetExhausted`]) nor an
/// exhausted rejection budget
/// ([`crate::NonlocalRelativityError::AdaptiveRejectionBudgetExhausted`])
/// is silently swallowed into a truncated or out-of-tolerance trajectory.
pub fn simulate_nonlocal_worldline_adaptive_with_stepper_policy<B, H, L, T, const D: usize>(
    background: &B,
    initial_state: WorldlineState<D>,
    config: AdaptiveNonlocalConfig,
    policy: AdaptiveStepperPolicy<H, L, T>,
) -> NonlocalResult<NonlocalTrajectory<D>>
where
    B: Metric<D> + Connection<D>,
    H: HistoryBackend<D>,
    L: MemoryLaw<D>,
    T: HistoryTransport<D>,
{
    validate_initial_state(&initial_state)?;

    let initial_metric = validated_metric(background, &initial_state.coordinates, 0)?;
    let initial_metric_norm = validated_metric_norm(
        &initial_metric,
        &initial_state.velocity,
        config.metric_norm_floor(),
        0,
    )?;

    let mut history = policy.history_backend;
    let memory_law = policy.memory_law;
    let transport = policy.transport;

    history.push_entry(
        background,
        &transport,
        HistoryEntry::new(initial_state.coordinates, initial_state.velocity, 0.0),
    )?;

    let mut states = vec![initial_state];
    let mut diagnostics_list = Vec::new();
    let mut history_diagnostics_list = Vec::new();

    let config_template = NonlocalConfig::from_fractional_order(
        config.fractional_order(),
        config.coupling(),
        config.initial_step(),
        1,
        config.metric_norm_floor(),
    )?;

    let initial_evaluation = evaluate_step_with_policy(StepEvaluationInput {
        background,
        state: &initial_state,
        history: &history,
        memory_law: &memory_law,
        transport: &transport,
        initial_metric_norm,
        affine_parameter: 0.0,
        step_index: 0,
        config: config_template,
    })?;
    diagnostics_list.push(initial_evaluation.diagnostics);
    history_diagnostics_list.push(history.diagnostics());

    let mut accepted_acceleration = initial_evaluation.acceleration;
    let mut current_parameter = 0.0_f64;
    let mut proposed_step = config.initial_step();
    let mut accepted_count = 0_usize;

    while current_parameter < config.target_affine_parameter()
    {
        if accepted_count >= config.max_accepted_steps()
        {
            return Err(NonlocalRelativityError::AdaptiveStepBudgetExhausted {
                accepted_steps: accepted_count,
                reached_parameter: current_parameter,
                target_affine_parameter: config.target_affine_parameter(),
            });
        }

        let remaining = config.target_affine_parameter() - current_parameter;
        let mut step = proposed_step.min(remaining);
        let step_index = accepted_count + 1;
        let mut rejection_count = 0_usize;

        let (accepted_state, used_step, next_proposed_step) = loop
        {
            let trial = step_doubling_trial(
                background,
                &states[accepted_count],
                &accepted_acceleration,
                &history,
                &memory_law,
                &transport,
                current_parameter,
                step,
                initial_metric_norm,
                config_template,
                step_index,
            )?;

            // Componentwise scaled RMS local error (shared by both adaptive
            // controllers); the one full step and the two half steps are the
            // step-doubling Richardson pair.
            let normalized_error = scaled_local_error_norm(
                &trial.full_state,
                &trial.refined_state,
                config.tolerance(),
            )?;

            if normalized_error <= 1.0
            {
                let growth = if normalized_error > 0.0
                {
                    STEP_SAFETY_FACTOR * normalized_error.powf(-STEP_DOUBLING_ERROR_EXPONENT)
                }
                else
                {
                    STEP_GROWTH_CAP
                };
                let next_step = (step * growth.min(STEP_GROWTH_CAP))
                    .clamp(config.min_step(), config.max_step());
                break (trial.refined_state, step, next_step);
            }

            rejection_count += 1;
            if rejection_count >= config.max_rejections_per_step()
            {
                return Err(NonlocalRelativityError::AdaptiveRejectionBudgetExhausted {
                    accepted_step: accepted_count,
                    attempted_step: step,
                    error_estimate: normalized_error,
                });
            }

            let shrink = STEP_SAFETY_FACTOR * normalized_error.powf(-STEP_DOUBLING_ERROR_EXPONENT);
            let shrunk_step = step * shrink.max(STEP_SHRINK_FLOOR);

            if shrunk_step < config.min_step()
            {
                return Err(NonlocalRelativityError::AdaptiveRejectionBudgetExhausted {
                    accepted_step: accepted_count,
                    attempted_step: step,
                    error_estimate: normalized_error,
                });
            }

            step = shrunk_step;
        };

        current_parameter += used_step;
        history.push_entry(
            background,
            &transport,
            HistoryEntry::new(
                accepted_state.coordinates,
                accepted_state.velocity,
                current_parameter,
            ),
        )?;
        states.push(accepted_state);

        let final_evaluation = evaluate_step_with_policy(StepEvaluationInput {
            background,
            state: &accepted_state,
            history: &history,
            memory_law: &memory_law,
            transport: &transport,
            initial_metric_norm,
            affine_parameter: current_parameter,
            step_index,
            config: config_template,
        })?;
        diagnostics_list.push(final_evaluation.diagnostics);
        history_diagnostics_list.push(history.diagnostics());

        accepted_acceleration = final_evaluation.acceleration;
        proposed_step = next_proposed_step;
        accepted_count += 1;
    }

    Ok(NonlocalTrajectory::new(
        states,
        diagnostics_list,
        history_diagnostics_list,
    ))
}

/// Simulate the experimental fractional-memory worldline model with an
/// adaptive, step-doubling-controlled affine-parameter step, using plain
/// non-uniform coordinate memory (complete history, no geometric transport).
///
/// This is exactly
/// [`simulate_nonlocal_worldline_adaptive_with_stepper_policy`] with
/// [`AdaptiveStepperPolicy::new`]`(CompleteUniformHistory::new(),
/// NonuniformCaputoCoordinateMemory, IdentityHistoryTransport)`; use the
/// `_with_stepper_policy` entry point directly to compose with
/// [`crate::NonuniformModulatedCaputoCoordinateMemory`],
/// [`crate::DiscreteConnectionTransport`], or
/// [`crate::BoundedShortMemoryHistory`].
pub fn simulate_nonlocal_worldline_adaptive_with_stepper<B, const D: usize>(
    background: &B,
    initial_state: WorldlineState<D>,
    config: AdaptiveNonlocalConfig,
) -> NonlocalResult<NonlocalTrajectory<D>>
where
    B: Metric<D> + Connection<D>,
{
    simulate_nonlocal_worldline_adaptive_with_stepper_policy(
        background,
        initial_state,
        config,
        AdaptiveStepperPolicy::new(
            CompleteUniformHistory::new(),
            NonuniformCaputoCoordinateMemory,
            IdentityHistoryTransport,
        ),
    )
}

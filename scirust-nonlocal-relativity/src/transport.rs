//! Typed history entries and discrete parallel-transport approximation.
//!
//! The Phase 2 [`crate::HistoryTransport`] contract only received a retained
//! sample's velocity components and the current worldline state: enough for
//! coordinate-identity transport, but not enough to carry a vector between
//! two distinct tangent spaces, since that requires the vector's *source*
//! point. [`HistoryEntry`] is the typed accepted sample (coordinates,
//! velocity, parameter) that supplies that source point, and
//! [`DiscreteConnectionTransport`] is a deterministic, explicitly
//! discretized approximation of parallel transport built on top of it.

use crate::{
    Connection, HistoryTransport, NonlocalRelativityError, NonlocalResult, WorldlineState,
    validate_history_velocity,
};

type TransportGenerator<const D: usize> = [[f64; D]; D];

/// One accepted worldline sample retained for history-dependent evaluation.
///
/// This is the typed replacement for a bare velocity component array: it
/// keeps the coordinates and parameter value where the sample was accepted,
/// which a geometric [`HistoryTransport`] needs in order to carry the
/// sample's vector from the tangent space where it was recorded into the
/// tangent space at a later worldline state. Nothing in this crate transports
/// a vector between distinct tangent spaces using components alone; when only
/// components are available (for example through the legacy
/// [`HistoryBackend::push_velocity`] path), transport is limited to the
/// coordinate-identity contract.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HistoryEntry<const D: usize> {
    /// Accepted coordinates `x^rho` at this sample.
    pub coordinates: [f64; D],
    /// Accepted contravariant velocity `u^rho` at this sample.
    pub velocity: [f64; D],
    /// Accepted affine or proper-time parameter value at this sample.
    pub parameter: f64,
}

impl<const D: usize> HistoryEntry<D> {
    /// Construct a history entry from coordinates, velocity, and parameter.
    #[must_use]
    pub const fn new(coordinates: [f64; D], velocity: [f64; D], parameter: f64) -> Self {
        Self {
            coordinates,
            velocity,
            parameter,
        }
    }

    /// Construct a velocity-only entry for callers that do not track source
    /// coordinates or a parameter value.
    ///
    /// Coordinates are set to the origin and the parameter to zero. These
    /// default values are never read by [`crate::IdentityHistoryTransport`] or
    /// by the coordinate-memory pipeline; a geometric transport must not be
    /// combined with entries constructed this way, since it would silently
    /// transport from the wrong point.
    #[must_use]
    pub const fn from_velocity_only(velocity: [f64; D]) -> Self {
        Self {
            coordinates: [0.0; D],
            velocity,
            parameter: 0.0,
        }
    }
}

/// Transport all entries in place from their current frame (the last
/// retained entry, when one exists) to `destination`, then leave `entries`
/// ready for the caller to push `destination` itself.
///
/// This is called once per accepted (or provisional) segment by
/// [`HistoryBackend::push_entry`] implementations. The generic batch contract is
/// transactional. [`DiscreteConnectionTransport`] specializes it by constructing
/// one segment operator in `O(D^3)` and applying it to `H` retained vectors in
/// `O(H * D^2)`, rather than recomputing Christoffel contractions `H` times.
pub(crate) fn transport_retained_entries<const D: usize, B, T>(
    entries: &mut [HistoryEntry<D>],
    background: &B,
    transport: &T,
    destination: &HistoryEntry<D>,
) -> NonlocalResult<()>
where
    B: Connection<D>,
    T: HistoryTransport<D>,
{
    let Some(last) = entries.last().copied()
    else
    {
        return Ok(());
    };

    let from_state = WorldlineState::new(last.coordinates, last.velocity);
    let to_state = WorldlineState::new(destination.coordinates, destination.velocity);
    let segment_step = destination.parameter - last.parameter;

    let mut velocities: Vec<[f64; D]> = entries.iter().map(|entry| entry.velocity).collect();
    transport.transport_batch(
        background,
        &mut velocities,
        &from_state,
        &to_state,
        segment_step,
    )?;
    for (retained_index, velocity) in velocities.iter().enumerate()
    {
        validate_transported_vector(velocity, retained_index)?;
    }
    for (entry, velocity) in entries.iter_mut().zip(velocities)
    {
        entry.velocity = velocity;
    }

    Ok(())
}

/// Transport `vector`, initially expressed at `waypoints[0]`'s point, along
/// the polyline described by `waypoints` (ordered oldest to newest), one
/// discrete segment at a time via `transport`. Returns the vector expressed
/// at `waypoints.last()`'s point.
///
/// This is the same per-segment mechanism [`HistoryBackend::push_entry`]
/// applies internally, exposed directly over an explicit path for research
/// and validation use: it lets a caller measure how a [`HistoryTransport`]'s
/// accumulated numerical error behaves under path refinement, independent of
/// the full simulation/backend machinery (see
/// `examples/exact_transport_convergence.rs`, which compares
/// [`DiscreteConnectionTransport`] against
/// [`crate::exact_cylindrical_minkowski_transport`]).
///
/// `waypoints` must contain at least one entry. A single-waypoint polyline
/// has no segments to transport across and returns `vector` unchanged.
pub fn transport_vector_along_polyline<const D: usize, B, T>(
    background: &B,
    transport: &T,
    vector: [f64; D],
    waypoints: &[HistoryEntry<D>],
) -> NonlocalResult<[f64; D]>
where
    B: Connection<D>,
    T: HistoryTransport<D>,
{
    if waypoints.is_empty()
    {
        return Err(NonlocalRelativityError::EmptyTransportPolyline);
    }

    let mut current = vector;

    for (segment_index, window) in waypoints.windows(2).enumerate()
    {
        let from_state = WorldlineState::new(window[0].coordinates, window[0].velocity);
        let to_state = WorldlineState::new(window[1].coordinates, window[1].velocity);
        let segment_step = window[1].parameter - window[0].parameter;

        current = transport.transport_segment(
            segment_index,
            background,
            current,
            &from_state,
            &to_state,
            segment_step,
        )?;
        validate_transported_vector(&current, segment_index)?;
    }

    Ok(current)
}

/// Deterministic discrete parallel-transport approximation for retained
/// history vectors.
///
/// Each accepted segment is transported with a single Heun
/// predict-evaluate-correct-evaluate step of the linear transport equation
///
/// `dV^mu / dlambda = - Gamma^mu_(alpha beta) u^alpha V^beta`:
///
/// 1. evaluate the transport derivative at the segment start;
/// 2. predict the vector at the segment end;
/// 3. evaluate the connection and velocity at the segment end;
/// 4. correct with the average of the two derivatives.
///
/// [`HistoryBackend::push_entry`] dispatches one batch per accepted segment.
/// This implementation reuses one segment operator across all currently
/// retained vectors, so transport still accumulates along the actual accepted
/// worldline polyline rather than jumping directly between a sample's original
/// recorded point and the current point.
///
/// This is a discrete numerical approximation to parallel transport along a
/// polyline. It is **not** an exact analytic bitensor propagator, **not** a
/// proof of covariance, and discretization error accumulates with the
/// segment step and the number of transported segments. The legacy
/// [`HistoryTransport::transport_velocity`] method has no source point to
/// transport from, so this type implements it as an identity passthrough;
/// its real work happens in
/// [`transport_segment`](HistoryTransport::transport_segment), which is
/// called once per accepted segment rather than once per memory evaluation.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct DiscreteConnectionTransport;

impl DiscreteConnectionTransport {
    fn segment_generators<B, const D: usize>(
        background: &B,
        from_state: &WorldlineState<D>,
        to_state: &WorldlineState<D>,
    ) -> NonlocalResult<(TransportGenerator<D>, TransportGenerator<D>)>
    where
        B: Connection<D>,
    {
        let start_symbols = background.christoffel(&from_state.coordinates);
        validate_transport_christoffel(&start_symbols)?;
        let end_symbols = background.christoffel(&to_state.coordinates);
        validate_transport_christoffel(&end_symbols)?;
        Ok((
            transport_generator(&start_symbols, &from_state.velocity),
            transport_generator(&end_symbols, &to_state.velocity),
        ))
    }

    /// Construct the linear Heun transport operator for one accepted segment.
    ///
    /// Let `A_0^mu_beta = Gamma^mu_(alpha beta)(x_0) u_0^alpha` and
    /// analogously `A_1` at the segment end. With `H_0 = h A_0` and
    /// `H_1 = h A_1`, the existing scalar Heun update is algebraically
    ///
    /// `V_1 = [I - 1/2 H_0 - 1/2 H_1 + 1/2 H_1 H_0] V_0`.
    ///
    /// The generators are scaled before multiplication. This avoids forming
    /// `A_1 A_0` first and only then multiplying by `h^2`, which can overflow
    /// even when the final Heun update is finite. Floating-point regrouping
    /// means results are numerically equivalent to, but not generally
    /// bit-identical with, [`HistoryTransport::transport_segment`].
    pub fn segment_operator<B, const D: usize>(
        background: &B,
        from_state: &WorldlineState<D>,
        to_state: &WorldlineState<D>,
        segment_step: f64,
    ) -> NonlocalResult<[[f64; D]; D]>
    where
        B: Connection<D>,
    {
        if !segment_step.is_finite()
        {
            return Err(NonlocalRelativityError::InvalidTransportSegmentStep(
                segment_step,
            ));
        }

        let (start_generator, end_generator) =
            Self::segment_generators(background, from_state, to_state)?;
        let mut start_scaled = start_generator;
        let mut end_scaled = end_generator;
        scale_transport_generator(&mut start_scaled, segment_step)?;
        scale_transport_generator(&mut end_scaled, segment_step)?;

        let mut operator = [[0.0_f64; D]; D];
        for (row_index, row) in operator.iter_mut().enumerate()
        {
            for (column_index, coefficient) in row.iter_mut().enumerate()
            {
                let mut product = 0.0_f64;
                for (end_value, start_row) in end_scaled[row_index].iter().zip(start_scaled.iter())
                {
                    product += *end_value * start_row[column_index];
                }
                let identity = if row_index == column_index { 1.0 } else { 0.0 };
                *coefficient = identity
                    - 0.5 * start_scaled[row_index][column_index]
                    - 0.5 * end_scaled[row_index][column_index]
                    + 0.5 * product;
                if !coefficient.is_finite()
                {
                    return Err(NonlocalRelativityError::NonFiniteTransportOperator {
                        row: row_index,
                        column: column_index,
                        value: *coefficient,
                    });
                }
            }
        }
        Ok(operator)
    }

    /// Apply a precomputed segment operator to one retained vector.
    pub fn apply_segment_operator<const D: usize>(
        operator: &[[f64; D]; D],
        vector: [f64; D],
        retained_index: usize,
    ) -> NonlocalResult<[f64; D]> {
        validate_history_velocity(&vector, retained_index)?;
        let mut transported = [0.0_f64; D];
        for (row, output) in operator.iter().zip(transported.iter_mut())
        {
            for (coefficient, value) in row.iter().zip(vector.iter())
            {
                *output += coefficient * value;
            }
        }
        validate_transported_vector(&transported, retained_index)?;
        Ok(transported)
    }
}

impl<const D: usize> HistoryTransport<D> for DiscreteConnectionTransport {
    fn transport_velocity(
        &self,
        retained_index: usize,
        velocity: [f64; D],
        _current_state: &WorldlineState<D>,
    ) -> NonlocalResult<[f64; D]> {
        validate_history_velocity(&velocity, retained_index)?;
        Ok(velocity)
    }

    fn transport_segment<B>(
        &self,
        retained_index: usize,
        background: &B,
        vector: [f64; D],
        from_state: &WorldlineState<D>,
        to_state: &WorldlineState<D>,
        segment_step: f64,
    ) -> NonlocalResult<[f64; D]>
    where
        B: Connection<D>,
    {
        heun_discrete_parallel_transport(
            retained_index,
            background,
            vector,
            from_state,
            to_state,
            segment_step,
        )
    }

    fn transport_batch<B>(
        &self,
        background: &B,
        vectors: &mut [[f64; D]],
        from_state: &WorldlineState<D>,
        to_state: &WorldlineState<D>,
        segment_step: f64,
    ) -> NonlocalResult<()>
    where
        B: Connection<D>,
    {
        if !segment_step.is_finite()
        {
            return Err(NonlocalRelativityError::InvalidTransportSegmentStep(
                segment_step,
            ));
        }
        let (start_generator, end_generator) =
            Self::segment_generators(background, from_state, to_state)?;
        let mut outputs = Vec::with_capacity(vectors.len());
        for (retained_index, vector) in vectors.iter().copied().enumerate()
        {
            outputs.push(heun_discrete_parallel_transport_from_generators(
                retained_index,
                vector,
                &start_generator,
                &end_generator,
                segment_step,
            )?);
        }
        vectors.copy_from_slice(&outputs);
        Ok(())
    }
}

fn heun_discrete_parallel_transport_from_generators<const D: usize>(
    retained_index: usize,
    vector: [f64; D],
    start_generator: &[[f64; D]; D],
    end_generator: &[[f64; D]; D],
    segment_step: f64,
) -> NonlocalResult<[f64; D]> {
    validate_history_velocity(&vector, retained_index)?;
    let start_derivative = generator_transport_derivative(start_generator, &vector);
    validate_transported_vector(&start_derivative, retained_index)?;

    let mut predicted = [0.0_f64; D];
    for mu in 0..D
    {
        predicted[mu] = vector[mu] + segment_step * start_derivative[mu];
    }
    validate_transported_vector(&predicted, retained_index)?;

    let end_derivative = generator_transport_derivative(end_generator, &predicted);
    validate_transported_vector(&end_derivative, retained_index)?;

    let mut corrected = [0.0_f64; D];
    for mu in 0..D
    {
        corrected[mu] =
            vector[mu] + 0.5 * segment_step * (start_derivative[mu] + end_derivative[mu]);
    }
    validate_transported_vector(&corrected, retained_index)?;
    Ok(corrected)
}

fn heun_discrete_parallel_transport<B, const D: usize>(
    retained_index: usize,
    background: &B,
    vector: [f64; D],
    from_state: &WorldlineState<D>,
    to_state: &WorldlineState<D>,
    segment_step: f64,
) -> NonlocalResult<[f64; D]>
where
    B: Connection<D>,
{
    if !segment_step.is_finite()
    {
        return Err(NonlocalRelativityError::InvalidTransportSegmentStep(
            segment_step,
        ));
    }

    let start_symbols = background.christoffel(&from_state.coordinates);
    validate_transport_christoffel(&start_symbols)?;
    let start_derivative =
        parallel_transport_derivative(&start_symbols, &from_state.velocity, &vector);
    validate_transported_vector(&start_derivative, retained_index)?;

    let mut predicted = [0.0_f64; D];
    for mu in 0..D
    {
        predicted[mu] = vector[mu] + segment_step * start_derivative[mu];
    }
    validate_transported_vector(&predicted, retained_index)?;

    let end_symbols = background.christoffel(&to_state.coordinates);
    validate_transport_christoffel(&end_symbols)?;
    let end_derivative =
        parallel_transport_derivative(&end_symbols, &to_state.velocity, &predicted);
    validate_transported_vector(&end_derivative, retained_index)?;

    let mut corrected = [0.0_f64; D];
    for mu in 0..D
    {
        corrected[mu] =
            vector[mu] + 0.5 * segment_step * (start_derivative[mu] + end_derivative[mu]);
    }
    validate_transported_vector(&corrected, retained_index)?;

    Ok(corrected)
}

/// Build `A^mu_beta = Gamma^mu_(alpha beta) u^alpha`, the linear generator
/// whose negative multiplies a vector in the parallel-transport equation.
fn transport_generator<const D: usize>(
    christoffel: &[[[f64; D]; D]; D],
    local_velocity: &[f64; D],
) -> [[f64; D]; D] {
    let mut generator = [[0.0_f64; D]; D];
    for (row_index, row) in generator.iter_mut().enumerate()
    {
        for (column_index, coefficient) in row.iter_mut().enumerate()
        {
            for (alpha, velocity) in local_velocity.iter().copied().enumerate()
            {
                *coefficient += christoffel[row_index][alpha][column_index] * velocity;
            }
        }
    }
    generator
}

fn scale_transport_generator<const D: usize>(
    generator: &mut [[f64; D]; D],
    segment_step: f64,
) -> NonlocalResult<()> {
    for (row_index, row) in generator.iter_mut().enumerate()
    {
        for (column_index, coefficient) in row.iter_mut().enumerate()
        {
            *coefficient *= segment_step;
            if !coefficient.is_finite()
            {
                return Err(NonlocalRelativityError::NonFiniteTransportOperator {
                    row: row_index,
                    column: column_index,
                    value: *coefficient,
                });
            }
        }
    }
    Ok(())
}

/// Evaluate `-A^mu_beta V^beta` after the connection/velocity contraction.
fn generator_transport_derivative<const D: usize>(
    generator: &[[f64; D]; D],
    vector: &[f64; D],
) -> [f64; D] {
    let mut derivative = [0.0_f64; D];
    for (row, output) in generator.iter().zip(derivative.iter_mut())
    {
        for (coefficient, value) in row.iter().zip(vector.iter())
        {
            *output -= coefficient * value;
        }
    }
    derivative
}

/// Evaluate `-Gamma^mu_(alpha beta) u^alpha V^beta`, the linear parallel
/// transport derivative of `vector` along local direction `local_velocity`.
fn parallel_transport_derivative<const D: usize>(
    christoffel: &[[[f64; D]; D]; D],
    local_velocity: &[f64; D],
    vector: &[f64; D],
) -> [f64; D] {
    let mut derivative = [0.0_f64; D];

    for mu in 0..D
    {
        for alpha in 0..D
        {
            for beta in 0..D
            {
                derivative[mu] -=
                    christoffel[mu][alpha][beta] * local_velocity[alpha] * vector[beta];
            }
        }
    }

    derivative
}

fn validate_transport_christoffel<const D: usize>(
    symbols: &[[[f64; D]; D]; D],
) -> NonlocalResult<()> {
    for (rho, rho_values) in symbols.iter().enumerate()
    {
        for (mu, mu_values) in rho_values.iter().enumerate()
        {
            for (nu, value) in mu_values.iter().copied().enumerate()
            {
                if !value.is_finite()
                {
                    return Err(NonlocalRelativityError::NonFiniteTransportChristoffel {
                        rho,
                        mu,
                        nu,
                        value,
                    });
                }
            }
        }
    }

    Ok(())
}

fn validate_transported_vector<const D: usize>(
    vector: &[f64; D],
    retained_index: usize,
) -> NonlocalResult<()> {
    for (component, value) in vector.iter().copied().enumerate()
    {
        if !value.is_finite()
        {
            return Err(NonlocalRelativityError::NonFiniteTransportedVector {
                retained_index,
                component,
                value,
            });
        }
    }

    Ok(())
}

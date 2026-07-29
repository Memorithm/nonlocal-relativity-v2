//! BSSN on a uniform periodic one-dimensional grid (Layer 3.4).
//!
//! Layer 3.3 delivered the BSSN algebra at a *point*, where every
//! spatial-derivative term vanished identically because every field was
//! spatially constant. This module supplies the grid, making BSSN a partial
//! differential equation for the first time.
//!
//! # 1D3V, not 3D
//!
//! Fields vary only along `x`, but every grid point stores complete
//! three-dimensional tensors: `gammatilde_ij` and `Atilde_ij` are full symmetric
//! 3x3 objects and `Gammatilde^i` is a full 3-vector. Derivatives along `y` and
//! `z` are exactly zero **by construction** — the provider is a function of `x`
//! alone, so a difference in `y` subtracts a value from itself.
//!
//! This is **not** a general three-dimensional numerical-relativity solver.
//!
//! # How this reuses Layer 3.3
//!
//! Every Layer 3.1/3.3 evaluator samples its fields *by coordinate* and takes
//! all spatial derivatives internally by nested central differences. Rather than
//! rewrite those equations, this module supplies fields that **read from the
//! grid**: [`GridMetric`] and [`GridCurvature`] answer a query at coordinate `x`
//! by reconstructing the ADM pair at the nearest grid point.
//!
//! For that to be a *grid* difference rather than nonsense, the internal step
//! must land exactly on grid points, so this module always evaluates with
//! `spatial_step = metric_step = dx`. Callers do not choose these; the grid
//! does.
//!
//! No BSSN equation, no Ricci engine, no constraint evaluator, and no
//! integrator is duplicated here.
//!
//! # Genuine BSSN form, and why it matters
//!
//! An earlier revision of this module reused Layer 3.3's [`conformal_ricci`]
//! unchanged. That computes `Rtilde_ij` as the *generic* Ricci tensor of
//! `gammatilde_ij` — correct as a tensor identity, but **not the BSSN system**.
//! It also inherited Layer 3.3's `d_t Gammatilde^i = 0`, which was right there
//! (every term carries a spatial gradient, and it had none) and wrong here.
//!
//! The result was BSSN *variables* carrying ADM's *principal part*, and it
//! behaved exactly like weakly hyperbolic ADM: `N >= 64` blew up before `t = 1`
//! at every Courant factor, worse as the grid was refined, and Kreiss-Oliger
//! dissipation could not cure it.
//!
//! Two things fix it, and both are needed:
//!
//! - [`crate::bssn::conformal_ricci_from_derivatives`] writes `Rtilde_ij` in
//!   BSSN form, so the term `gammatilde_{k(i} d_{j)} Gammatilde^k` carries the
//!   **evolved** connection and removes the mixed second derivatives from the
//!   principal part, leaving the manifestly elliptic
//!   `-1/2 gammatilde^{lm} d_l d_m`.
//! - [`crate::bssn::bssn_connection_rhs`] supplies the `d_t Gammatilde^i`
//!   equation that Layer 3.3 explicitly deferred.
//!
//! Evolving `Gammatilde^i` without the first change does not help, because the
//! Ricci tensor never reads it.
//!
//! Measured afterwards: every resolution from `N = 32` to `N = 256` reaches
//! `t = 4` with the error *decreasing* under refinement; the connection
//! constraint grows only secularly and is resolution-independent; and the
//! Courant boundary is a genuine, resolution-independent `C <= 1` with `C = 2`
//! rejected at both `N = 64` and `N = 128`. Dissipation becomes very nearly a
//! no-op, which is the correct behaviour when there is no high-frequency growth
//! to damp.
//!
//! **None of this proves strong hyperbolicity.** That is an analytic property of
//! the continuum system; these are measurements on one discretisation of a
//! one-dimensional reduction.
//!
//! # Flat state layout
//!
//! The method-of-lines state is `17 N` values in structure-of-arrays order;
//! component `c` occupies the contiguous slice `[c*N, (c+1)*N)`:
//!
//! ```text
//!  0      phi
//!  1..7   gammatilde_xx, _xy, _xz, _yy, _yz, _zz
//!  7      K
//!  8..14  Atilde_xx, _xy, _xz, _yy, _yz, _zz
//! 14..17  Gammatilde^x, ^y, ^z
//! ```
//!
//! Symmetric tensors store only their six independent components, so symmetry
//! is preserved *structurally* — a round trip cannot break it, because the
//! redundant entries do not exist.
//!
//! # Gauge and constraints
//!
//! Prescribed gauge only (`alpha = 1`, `beta^i = 0`); gauge is **not evolved**.
//! Evolution is free: constraints are monitored and never enforced. Projection
//! is explicit and opt-in, and is never applied between RK4 stages.
//!
//! # Example
//!
//! ```
//! use scirust_relativity::adm_evolution::AdmSources;
//! use scirust_relativity::bssn_grid::{
//!     BssnGridState, BssnGridSystem, bssn_grid_constraints, evolve_bssn_grid,
//! };
//! use scirust_relativity::grid1d::UniformGrid1d;
//!
//! let grid = UniformGrid1d::new(16, 0.0, 1.0).expect("valid grid");
//! let system = BssnGridSystem::vacuum(grid);
//! let initial = BssnGridState::minkowski(grid);
//!
//! let samples = evolve_bssn_grid(&system, &initial, 0.0, 1.0, 0.25 * grid.spacing())
//!     .expect("Minkowski evolves");
//!
//! // Minkowski is exactly stationary: the right-hand side is exactly zero, so
//! // RK4 adds exactly zero and there is no rounding to accumulate.
//! let final_state = &samples.last().expect("non-empty").state;
//! assert_eq!(final_state.as_slice(), initial.as_slice());
//!
//! let constraints =
//!     bssn_grid_constraints(&final_state.view(), &AdmSources::VACUUM).expect("constraints");
//! assert_eq!(constraints.hamiltonian.max_abs, 0.0);
//! ```

use core::fmt;

use scirust_sim::{SimError, System, simulate};

use crate::Metric;
use crate::adm_evolution::{
    AdmEvolutionError, AdmEvolutionSettings, AdmSources, SpatialTensorField,
    hamiltonian_constraint, momentum_constraint,
};
use crate::bssn::{
    BssnConnectionRhs, BssnError, BssnGauge, BssnLapseDerivatives, BssnProjection,
    BssnSecondDerivatives, BssnSpatialDerivatives, BssnState, BssnSuppliedTerms, ConformalRicci,
    adm_to_bssn, bssn_algebraic_constraints, bssn_connection_rhs, bssn_evolution_rhs_with_ricci,
    bssn_to_adm, conformal_ricci, conformal_ricci_from_derivatives, one_plus_log_lapse_rhs,
    project_trace_free, project_unit_determinant,
};
use crate::grid1d::{Grid1dError, GridReduction, UniformGrid1d};

/// Evolved scalar component arrays per grid point.
pub const COMPONENTS_PER_POINT: usize = 18;

/// Slot of `phi` in the flat layout.
pub const SLOT_CONFORMAL_FACTOR: usize = 0;
/// First slot of `gammatilde_ij` (six independent components).
pub const SLOT_CONFORMAL_METRIC: usize = 1;
/// Slot of `K`.
pub const SLOT_MEAN_CURVATURE: usize = 7;
/// First slot of `Atilde_ij` (six independent components).
pub const SLOT_CONFORMAL_CURVATURE: usize = 8;
/// First slot of `Gammatilde^i` (three components).
pub const SLOT_CONFORMAL_CONNECTION: usize = 14;
/// Slot of the lapse `alpha`.
///
/// Stored at every point whether or not it is evolved: a prescribed lapse is
/// simply one whose right-hand side is zero, which keeps the state layout
/// independent of the slicing choice.
pub const SLOT_LAPSE: usize = 17;

/// The independent components of a symmetric 3x3 tensor, in storage order.
const SYMMETRIC_PAIRS: [(usize, usize); 6] = [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)];

/// A failure in the grid BSSN pipeline, carrying the context needed to locate it.
#[derive(Debug, Clone, PartialEq)]
pub enum BssnGridError {
    /// The grid itself is invalid.
    Grid(Grid1dError),
    /// The flat state length does not match `17 N`.
    FlatStateLength {
        /// The length `17 N` the grid requires.
        expected: usize,
        /// The length supplied.
        actual: usize,
    },
    /// A stored field failed validation at a specific point.
    InvalidState {
        /// The coordinate time at which the failure was detected.
        time: f64,
        /// The grid index.
        index: usize,
        /// The field name.
        field: &'static str,
        /// The offending value.
        value: f64,
        /// What was wrong with it.
        category: StateFailure,
    },
    /// A pointwise BSSN evaluation failed.
    Pointwise {
        /// The coordinate time.
        time: f64,
        /// The grid index.
        index: usize,
        /// What was being evaluated.
        stage: &'static str,
        /// The underlying Layer 3.3 error.
        source: BssnError,
    },
    /// A constraint reconstruction failed.
    Constraint {
        /// The grid index.
        index: usize,
        /// Which constraint.
        constraint: &'static str,
        /// The underlying Layer 3.1 error.
        source: AdmEvolutionError,
    },
    /// An explicit projection failed.
    Projection {
        /// The grid index.
        index: usize,
        /// The underlying Layer 3.3 error.
        source: BssnError,
    },
    /// The integrator rejected the run or stopped on a non-finite state.
    Integration(SimError),
    /// A requested evolution parameter is invalid.
    InvalidRequest {
        /// Which parameter.
        parameter: &'static str,
        /// The offending value.
        value: f64,
    },
}

/// Why a stored field failed validation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StateFailure {
    /// The value is `NaN` or infinite.
    NonFinite,
    /// A quantity required to be strictly positive was not.
    NonPositive,
    /// The conformal metric could not be inverted.
    SingularConformalMetric,
    /// The reconstructed physical metric has a non-positive determinant.
    SingularPhysicalMetric,
}

impl fmt::Display for StateFailure {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self
        {
            Self::NonFinite => write!(f, "not finite"),
            Self::NonPositive => write!(f, "not strictly positive"),
            Self::SingularConformalMetric => write!(f, "singular conformal metric"),
            Self::SingularPhysicalMetric => write!(f, "singular physical metric"),
        }
    }
}

impl fmt::Display for BssnGridError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self
        {
            Self::Grid(error) => write!(f, "grid: {error}"),
            Self::FlatStateLength { expected, actual } =>
            {
                write!(f, "flat state must have {expected} values, got {actual}")
            },
            Self::InvalidState {
                time,
                index,
                field,
                value,
                category,
            } => write!(
                f,
                "at t = {time}, point {index}: field '{field}' = {value} is {category}"
            ),
            Self::Pointwise {
                time,
                index,
                stage,
                source,
            } => write!(f, "at t = {time}, point {index}: {stage} failed: {source}"),
            Self::Constraint {
                index,
                constraint,
                source,
            } => write!(f, "point {index}: {constraint} constraint failed: {source}"),
            Self::Projection { index, source } =>
            {
                write!(f, "point {index}: projection failed: {source}")
            },
            Self::Integration(error) => write!(f, "integration failed: {error}"),
            Self::InvalidRequest { parameter, value } =>
            {
                write!(f, "parameter '{parameter}' is invalid: {value}")
            },
        }
    }
}

impl std::error::Error for BssnGridError {}

impl From<Grid1dError> for BssnGridError {
    fn from(error: Grid1dError) -> Self {
        Self::Grid(error)
    }
}

// ---------------------------------------------------------------------------
// Storage
// ---------------------------------------------------------------------------

/// Read the BSSN state stored at `index` of a flat structure-of-arrays slice.
fn state_from_flat(flat: &[f64], points: usize, index: usize) -> BssnState {
    let at = |slot: usize| flat[slot * points + index];

    let mut conformal_metric = [[0.0_f64; 3]; 3];
    let mut conformal_curvature = [[0.0_f64; 3]; 3];
    for (slot, &(i, j)) in SYMMETRIC_PAIRS.iter().enumerate()
    {
        let metric_value = at(SLOT_CONFORMAL_METRIC + slot);
        let curvature_value = at(SLOT_CONFORMAL_CURVATURE + slot);
        conformal_metric[i][j] = metric_value;
        conformal_metric[j][i] = metric_value;
        conformal_curvature[i][j] = curvature_value;
        conformal_curvature[j][i] = curvature_value;
    }

    BssnState {
        conformal_factor: at(SLOT_CONFORMAL_FACTOR),
        conformal_metric,
        mean_curvature: at(SLOT_MEAN_CURVATURE),
        conformal_curvature,
        conformal_connection: [
            at(SLOT_CONFORMAL_CONNECTION),
            at(SLOT_CONFORMAL_CONNECTION + 1),
            at(SLOT_CONFORMAL_CONNECTION + 2),
        ],
    }
}

/// Write a BSSN state into a flat structure-of-arrays slice at `index`.
///
/// Only the six independent components of each symmetric tensor are stored, so
/// the symmetry cannot be lost.
fn state_into_flat(flat: &mut [f64], points: usize, index: usize, state: &BssnState) {
    flat[SLOT_CONFORMAL_FACTOR * points + index] = state.conformal_factor;
    flat[SLOT_MEAN_CURVATURE * points + index] = state.mean_curvature;
    for (slot, &(i, j)) in SYMMETRIC_PAIRS.iter().enumerate()
    {
        flat[(SLOT_CONFORMAL_METRIC + slot) * points + index] = state.conformal_metric[i][j];
        flat[(SLOT_CONFORMAL_CURVATURE + slot) * points + index] = state.conformal_curvature[i][j];
    }
    for component in 0..3
    {
        flat[(SLOT_CONFORMAL_CONNECTION + component) * points + index] =
            state.conformal_connection[component];
    }
}

/// A borrowed, read-only view of a flat BSSN grid state.
///
/// Constructing one allocates nothing: the flat method-of-lines array *is* the
/// storage, so the right-hand-side path never copies the grid.
#[derive(Debug, Clone, Copy)]
pub struct BssnGridView<'a> {
    grid: UniformGrid1d,
    components: &'a [f64],
}

impl<'a> BssnGridView<'a> {
    /// Borrow `components` as a grid state, checking only the length.
    pub fn new(grid: UniformGrid1d, components: &'a [f64]) -> Result<Self, BssnGridError> {
        let expected = COMPONENTS_PER_POINT * grid.points();
        if components.len() != expected
        {
            return Err(BssnGridError::FlatStateLength {
                expected,
                actual: components.len(),
            });
        }
        Ok(Self { grid, components })
    }

    /// The underlying grid.
    #[must_use]
    pub const fn grid(&self) -> &UniformGrid1d {
        &self.grid
    }

    /// The raw flat slice.
    #[must_use]
    pub const fn as_slice(&self) -> &'a [f64] {
        self.components
    }

    /// The BSSN state stored at `index` (wrapped periodically).
    #[must_use]
    pub fn state_at(&self, index: usize) -> BssnState {
        state_from_flat(
            self.components,
            self.grid.points(),
            self.grid.wrap_usize(index),
        )
    }

    /// The stored lapse at `index`.
    #[must_use]
    pub fn lapse_at(&self, index: usize) -> f64 {
        self.components[SLOT_LAPSE * self.grid.points() + self.grid.wrap_usize(index)]
    }

    /// The prescribed-or-evolved gauge at `index` (zero shift throughout).
    #[must_use]
    pub fn gauge_at(&self, index: usize) -> BssnGauge {
        BssnGauge {
            lapse: self.lapse_at(index),
            shift: [0.0; 3],
        }
    }

    /// The lapse gradient and coordinate Hessian at `index`, from compact
    /// stencils on the stored lapse array.
    #[must_use]
    pub fn lapse_derivatives_at(&self, index: usize) -> BssnLapseDerivatives {
        let mut hessian = [[0.0_f64; 3]; 3];
        hessian[0][0] = component_curvature(self, SLOT_LAPSE, index);
        BssnLapseDerivatives {
            gradient: [component_gradient(self, SLOT_LAPSE, index), 0.0, 0.0],
            hessian,
        }
    }

    /// The finite-difference settings this grid requires.
    ///
    /// Both steps equal `dx` so that every internally differenced sample lands
    /// exactly on a grid point. This is not a tuning parameter.
    #[must_use]
    pub fn settings(&self) -> AdmEvolutionSettings {
        AdmEvolutionSettings {
            spatial_step: self.grid.spacing(),
            metric_step: self.grid.spacing(),
        }
    }

    /// The metric field backed by this grid.
    #[must_use]
    pub const fn metric(&self) -> GridMetric<'_, 'a> {
        GridMetric { view: self }
    }

    /// The extrinsic-curvature field backed by this grid.
    #[must_use]
    pub const fn curvature(&self) -> GridCurvature<'_, 'a> {
        GridCurvature { view: self }
    }
}

/// An owned BSSN grid state.
#[derive(Debug, Clone, PartialEq)]
pub struct BssnGridState {
    grid: UniformGrid1d,
    components: Vec<f64>,
}

impl BssnGridState {
    /// Allocate a zeroed state on `grid`.
    ///
    /// A zeroed state is *not* physically valid (its conformal metric is
    /// singular); use [`Self::minkowski`] or [`Self::from_adm_fields`] to build
    /// usable initial data.
    #[must_use]
    pub fn zeroed(grid: UniformGrid1d) -> Self {
        Self {
            components: vec![0.0; COMPONENTS_PER_POINT * grid.points()],
            grid,
        }
    }

    /// Adopt an existing flat array as a grid state.
    pub fn from_flat(grid: UniformGrid1d, components: Vec<f64>) -> Result<Self, BssnGridError> {
        let expected = COMPONENTS_PER_POINT * grid.points();
        if components.len() != expected
        {
            return Err(BssnGridError::FlatStateLength {
                expected,
                actual: components.len(),
            });
        }
        Ok(Self { grid, components })
    }

    /// Stationary Minkowski: identity conformal metric, everything else zero.
    #[must_use]
    pub fn minkowski(grid: UniformGrid1d) -> Self {
        let mut state = Self::zeroed(grid);
        let identity = BssnState {
            conformal_factor: 0.0,
            conformal_metric: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            mean_curvature: 0.0,
            conformal_curvature: [[0.0; 3]; 3],
            conformal_connection: [0.0; 3],
        };
        for index in 0..state.grid.points()
        {
            state.set_state_at(index, &identity);
            state.set_lapse_at(index, 1.0);
        }
        state
    }

    /// Build initial data by converting analytic ADM fields at every grid point
    /// through the **production** [`adm_to_bssn`] path.
    ///
    /// `Gammatilde^i` is therefore seeded by exactly the reconstruction the
    /// constraint monitor later uses, so the connection constraint starts at
    /// machine zero by construction rather than at truncation error.
    ///
    /// The supplied fields are sampled only at grid points and at offsets of
    /// whole multiples of `dx`, so an analytic field and its grid restriction
    /// give identical results.
    pub fn from_adm_fields<G: Metric<3>>(
        grid: UniformGrid1d,
        spatial_metric: &G,
        extrinsic_curvature: &impl SpatialTensorField,
    ) -> Result<Self, BssnGridError> {
        let settings = AdmEvolutionSettings {
            spatial_step: grid.spacing(),
            metric_step: grid.spacing(),
        };
        let mut state = Self::zeroed(grid);
        for index in 0..grid.points()
        {
            let at = [grid.coordinate(index), 0.0, 0.0];
            let local = adm_to_bssn(spatial_metric, extrinsic_curvature, &at, &settings).map_err(
                |source| BssnGridError::Pointwise {
                    time: 0.0,
                    index,
                    stage: "adm_to_bssn",
                    source,
                },
            )?;
            state.set_state_at(index, &local);
            // Unit lapse unless a caller overwrites it: the same synchronous
            // gauge every earlier increment used.
            state.set_lapse_at(index, 1.0);
        }
        Ok(state)
    }

    /// The underlying grid.
    #[must_use]
    pub const fn grid(&self) -> &UniformGrid1d {
        &self.grid
    }

    /// The flat method-of-lines array.
    #[must_use]
    pub fn as_slice(&self) -> &[f64] {
        &self.components
    }

    /// The flat method-of-lines array, mutably.
    pub fn as_mut_slice(&mut self) -> &mut [f64] {
        &mut self.components
    }

    /// Borrow as a read-only view.
    #[must_use]
    pub fn view(&self) -> BssnGridView<'_> {
        BssnGridView {
            grid: self.grid,
            components: &self.components,
        }
    }

    /// The BSSN state stored at `index`.
    #[must_use]
    pub fn state_at(&self, index: usize) -> BssnState {
        state_from_flat(
            &self.components,
            self.grid.points(),
            self.grid.wrap_usize(index),
        )
    }

    /// The stored lapse at `index`.
    #[must_use]
    pub fn lapse_at(&self, index: usize) -> f64 {
        self.components[SLOT_LAPSE * self.grid.points() + self.grid.wrap_usize(index)]
    }

    /// Overwrite the stored lapse at `index`.
    pub fn set_lapse_at(&mut self, index: usize, lapse: f64) {
        let points = self.grid.points();
        let wrapped = self.grid.wrap_usize(index);
        self.components[SLOT_LAPSE * points + wrapped] = lapse;
    }

    /// Store `state` at `index`.
    pub fn set_state_at(&mut self, index: usize, state: &BssnState) {
        let points = self.grid.points();
        let wrapped = self.grid.wrap_usize(index);
        state_into_flat(&mut self.components, points, wrapped, state);
    }

    /// Validate every stored field, returning the first failure in ascending
    /// index order.
    ///
    /// Nothing is clamped, replaced, or repaired.
    pub fn validate(&self, time: f64) -> Result<(), BssnGridError> {
        for index in 0..self.grid.points()
        {
            let state = self.state_at(index);
            let fail = |field: &'static str, value: f64, category: StateFailure| {
                BssnGridError::InvalidState {
                    time,
                    index,
                    field,
                    value,
                    category,
                }
            };

            if !state.conformal_factor.is_finite()
            {
                return Err(fail("phi", state.conformal_factor, StateFailure::NonFinite));
            }
            // `phi` is canonical, so `chi = e^{-4 phi}` is derived and can only
            // fail by overflowing to zero or infinity.
            let chi = state.chi();
            if !chi.is_finite() || chi <= 0.0
            {
                return Err(fail("chi", chi, StateFailure::NonPositive));
            }
            if !state.mean_curvature.is_finite()
            {
                return Err(fail("K", state.mean_curvature, StateFailure::NonFinite));
            }
            let lapse = self.lapse_at(index);
            if !lapse.is_finite()
            {
                return Err(fail("alpha", lapse, StateFailure::NonFinite));
            }
            if lapse <= 0.0
            {
                return Err(fail("alpha", lapse, StateFailure::NonPositive));
            }
            for &(i, j) in &SYMMETRIC_PAIRS
            {
                if !state.conformal_metric[i][j].is_finite()
                {
                    return Err(fail(
                        "gammatilde",
                        state.conformal_metric[i][j],
                        StateFailure::NonFinite,
                    ));
                }
                if !state.conformal_curvature[i][j].is_finite()
                {
                    return Err(fail(
                        "Atilde",
                        state.conformal_curvature[i][j],
                        StateFailure::NonFinite,
                    ));
                }
            }
            for component in 0..3
            {
                if !state.conformal_connection[component].is_finite()
                {
                    return Err(fail(
                        "Gammatilde",
                        state.conformal_connection[component],
                        StateFailure::NonFinite,
                    ));
                }
            }
            if state.inverse_conformal_metric().is_err()
            {
                return Err(fail(
                    "gammatilde",
                    f64::NAN,
                    StateFailure::SingularConformalMetric,
                ));
            }
            match bssn_to_adm(&state)
            {
                Ok(adm) =>
                {
                    let determinant = crate::determinant::<3>(&adm.spatial_metric);
                    match determinant
                    {
                        Ok(value) if value.is_finite() && value > 0.0 =>
                        {},
                        Ok(value) =>
                        {
                            return Err(fail("gamma", value, StateFailure::SingularPhysicalMetric));
                        },
                        Err(_) =>
                        {
                            return Err(fail(
                                "gamma",
                                f64::NAN,
                                StateFailure::SingularPhysicalMetric,
                            ));
                        },
                    }
                },
                Err(_) =>
                {
                    return Err(fail(
                        "gamma",
                        f64::NAN,
                        StateFailure::SingularPhysicalMetric,
                    ));
                },
            }
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Grid-backed fields
// ---------------------------------------------------------------------------

/// The physical spatial metric `gamma_ij` reconstructed from the grid.
///
/// A query at coordinate `x` is answered from the **nearest grid point**. The
/// caller is responsible for only ever asking about coordinates that are whole
/// multiples of `dx` away from a grid point; [`BssnGridView::settings`]
/// guarantees exactly that for every Layer 3.1/3.3 evaluator.
///
/// A reconstruction failure emits `NaN` rather than a silent substitute, so the
/// downstream evaluator fails loudly instead of returning a wrong number.
#[derive(Debug, Clone, Copy)]
pub struct GridMetric<'v, 'a> {
    view: &'v BssnGridView<'a>,
}

impl Metric<3> for GridMetric<'_, '_> {
    fn components(&self, coordinates: &[f64; 3]) -> [[f64; 3]; 3] {
        let index = self.view.grid.nearest_index(coordinates[0]);
        match bssn_to_adm(&self.view.state_at(index))
        {
            Ok(adm) => adm.spatial_metric,
            Err(_) => [[f64::NAN; 3]; 3],
        }
    }
}

/// The extrinsic curvature `K_ij` reconstructed from the grid.
///
/// Same nearest-point contract and same loud-failure policy as [`GridMetric`].
#[derive(Debug, Clone, Copy)]
pub struct GridCurvature<'v, 'a> {
    view: &'v BssnGridView<'a>,
}

impl SpatialTensorField for GridCurvature<'_, '_> {
    fn components(&self, coordinates: &[f64; 3]) -> [[f64; 3]; 3] {
        let index = self.view.grid.nearest_index(coordinates[0]);
        match bssn_to_adm(&self.view.state_at(index))
        {
            Ok(adm) => adm.extrinsic_curvature,
            Err(_) => [[f64::NAN; 3]; 3],
        }
    }
}

// ---------------------------------------------------------------------------
// Right-hand side
// ---------------------------------------------------------------------------

/// The centred periodic first difference of one component array at `index`.
///
/// Compact: it samples `i-1` and `i+1` only, so the stencil is one grid spacing
/// wide. This is what makes the connection evolution well-behaved where the
/// nested coordinate-sampled path is not.
fn component_gradient(view: &BssnGridView<'_>, slot: usize, index: usize) -> f64 {
    let points = view.grid.points();
    let base = slot * points;
    let left = view.components[base + view.grid.offset(index, -1)];
    let right = view.components[base + view.grid.offset(index, 1)];
    (right - left) / (2.0 * view.grid.spacing())
}

/// The spatial gradients at `index`, taken with compact stencils on the stored
/// component arrays.
///
/// **1D3V**: only the `x` gradient can be non-zero. The `y` and `z` entries are
/// exactly zero because the grid has no extent in those directions — a
/// structural fact, not an approximation.
#[must_use]
pub fn grid_spatial_derivatives(view: &BssnGridView<'_>, index: usize) -> BssnSpatialDerivatives {
    let mut conformal_metric_gradient = [[[0.0_f64; 3]; 3]; 3];
    for (slot, &(i, j)) in SYMMETRIC_PAIRS.iter().enumerate()
    {
        let value = component_gradient(view, SLOT_CONFORMAL_METRIC + slot, index);
        conformal_metric_gradient[0][i][j] = value;
        conformal_metric_gradient[0][j][i] = value;
    }

    BssnSpatialDerivatives {
        conformal_factor_gradient: [
            component_gradient(view, SLOT_CONFORMAL_FACTOR, index),
            0.0,
            0.0,
        ],
        mean_curvature_gradient: [
            component_gradient(view, SLOT_MEAN_CURVATURE, index),
            0.0,
            0.0,
        ],
        conformal_metric_gradient,
    }
}

/// The centred periodic second difference of one component array at `index`.
///
/// Compact: `f_{i+1} - 2 f_i + f_{i-1}` over `dx^2`. Its Fourier symbol
/// `2 cos(theta) - 2` is maximal at the Nyquist mode rather than vanishing
/// there, unlike the `2 dx`-spaced difference the coordinate-sampled path
/// produces.
fn component_curvature(view: &BssnGridView<'_>, slot: usize, index: usize) -> f64 {
    let points = view.grid.points();
    let base = slot * points;
    let centre = view.components[base + view.grid.wrap_usize(index)];
    let left = view.components[base + view.grid.offset(index, -1)];
    let right = view.components[base + view.grid.offset(index, 1)];
    (right - 2.0 * centre + left) / (view.grid.spacing() * view.grid.spacing())
}

/// The second spatial derivatives at `index`, taken with compact stencils.
///
/// **1D3V**: only the `xx` entries can be non-zero.
#[must_use]
pub fn grid_second_derivatives(view: &BssnGridView<'_>, index: usize) -> BssnSecondDerivatives {
    let mut conformal_metric_hessian = [[[[0.0_f64; 3]; 3]; 3]; 3];
    for (slot, &(i, j)) in SYMMETRIC_PAIRS.iter().enumerate()
    {
        let value = component_curvature(view, SLOT_CONFORMAL_METRIC + slot, index);
        conformal_metric_hessian[0][0][i][j] = value;
        conformal_metric_hessian[0][0][j][i] = value;
    }

    let mut conformal_factor_hessian = [[0.0_f64; 3]; 3];
    conformal_factor_hessian[0][0] = component_curvature(view, SLOT_CONFORMAL_FACTOR, index);

    let mut conformal_connection_gradient = [[0.0_f64; 3]; 3];
    // Indexed rather than iterated: the loop walks tensor components and the
    // matching storage slots together, which the index expresses.
    #[allow(clippy::needless_range_loop)]
    for component in 0..3
    {
        conformal_connection_gradient[component][0] =
            component_gradient(view, SLOT_CONFORMAL_CONNECTION + component, index);
    }

    BssnSecondDerivatives {
        conformal_factor_hessian,
        conformal_metric_hessian,
        conformal_connection_gradient,
    }
}

/// The conformal Ricci tensor at one grid point, in genuine BSSN form.
///
/// Uses the **evolved** `Gammatilde^i` and compact stencils, so the principal
/// part is the manifestly elliptic `-1/2 gammatilde^{lm} d_l d_m`.
pub fn grid_conformal_ricci(
    view: &BssnGridView<'_>,
    index: usize,
) -> Result<ConformalRicci, BssnGridError> {
    let first = grid_spatial_derivatives(view, index);
    let second = grid_second_derivatives(view, index);
    conformal_ricci_from_derivatives(&view.state_at(index), &first, &second).map_err(|source| {
        BssnGridError::Pointwise {
            time: 0.0,
            index,
            stage: "conformal_ricci_from_derivatives",
            source,
        }
    })
}

/// Evaluate the conformal connection right-hand side at one grid point.
pub fn grid_connection_rhs(
    view: &BssnGridView<'_>,
    index: usize,
    gauge: &BssnGauge,
    sources: &AdmSources,
) -> Result<BssnConnectionRhs, BssnGridError> {
    let derivatives = grid_spatial_derivatives(view, index);
    let lapse = view.lapse_derivatives_at(index);
    bssn_connection_rhs(&view.state_at(index), &derivatives, &lapse, gauge, sources).map_err(
        |source| BssnGridError::Pointwise {
            time: 0.0,
            index,
            stage: "bssn_connection_rhs",
            source,
        },
    )
}

/// Evaluate the BSSN right-hand side over the whole grid into `out`.
///
/// Points are visited in ascending index order. Nothing is allocated per point:
/// the input is borrowed and the output is written in place.
pub fn bssn_grid_rhs(
    view: &BssnGridView<'_>,
    gauge: &BssnGauge,
    sources: &AdmSources,
    time: f64,
    out: &mut [f64],
) -> Result<(), BssnGridError> {
    bssn_grid_rhs_with_slicing(view, gauge, sources, BssnSlicing::Prescribed, time, out)
}

/// Evaluate the BSSN right-hand side with an explicit slicing condition.
///
/// The lapse used at each point is the **stored** one, not the `gauge`
/// argument's; `gauge` supplies only the shift, which is zero throughout this
/// increment.
pub fn bssn_grid_rhs_with_slicing(
    view: &BssnGridView<'_>,
    gauge: &BssnGauge,
    sources: &AdmSources,
    slicing: BssnSlicing,
    time: f64,
    out: &mut [f64],
) -> Result<(), BssnGridError> {
    let points = view.grid.points();
    let expected = COMPONENTS_PER_POINT * points;
    if out.len() != expected
    {
        return Err(BssnGridError::FlatStateLength {
            expected,
            actual: out.len(),
        });
    }

    let settings = view.settings();
    let metric = view.metric();
    let curvature = view.curvature();

    for index in 0..points
    {
        let at = [view.grid.coordinate(index), 0.0, 0.0];
        // The conformal Ricci comes from the genuine BSSN form -- the evolved
        // Gammatilde^k with compact stencils -- not from a generic metric Ricci.
        // That substitution is what removes the mixed second derivatives from
        // the principal part.
        let ricci = grid_conformal_ricci(view, index)?;
        // The lapse is a stored field, so the gauge is per-point. The shift
        // comes from the caller and is zero throughout this increment.
        let local_gauge = BssnGauge {
            lapse: view.lapse_at(index),
            shift: gauge.shift,
        };
        let first = grid_spatial_derivatives(view, index);
        let lapse_derivatives = view.lapse_derivatives_at(index);
        let (_, rhs) = bssn_evolution_rhs_with_ricci(
            &metric,
            &curvature,
            &at,
            &local_gauge,
            sources,
            &settings,
            &BssnSuppliedTerms {
                ricci: &ricci,
                derivatives: &first,
                lapse: &lapse_derivatives,
            },
        )
        .map_err(|source| BssnGridError::Pointwise {
            time,
            index,
            stage: "bssn_evolution_rhs_with_ricci",
            source,
        })?;

        // Layer 3.3 reports `d_t Gammatilde^i = 0`: every term of that equation
        // carries a spatial gradient, and it had none. On a grid those gradients
        // are non-zero, so the equation is supplied here from compact centred
        // differences of the stored arrays. Leaving it at zero freezes
        // `Gammatilde^i` while `gammatilde_ij` evolves, which is what reduces
        // BSSN back to a weakly hyperbolic ADM-like system.
        let connection = grid_connection_rhs(view, index, &local_gauge, sources)?;

        // The right-hand side has exactly the shape of a state, so it is stored
        // through the same layout: no second flattening convention exists.
        let as_state = BssnState {
            conformal_factor: rhs.conformal_factor,
            conformal_metric: rhs.conformal_metric,
            mean_curvature: rhs.mean_curvature,
            conformal_curvature: rhs.conformal_curvature,
            conformal_connection: connection.total,
        };
        state_into_flat(out, points, index, &as_state);

        // The slicing condition. `Prescribed` keeps the lapse fixed, so its
        // right-hand side is exactly zero.
        out[SLOT_LAPSE * points + index] = match slicing
        {
            BssnSlicing::Prescribed => 0.0,
            BssnSlicing::OnePlusLog => one_plus_log_lapse_rhs(&view.state_at(index), &local_gauge),
        };
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Method of lines
// ---------------------------------------------------------------------------

/// Which slicing condition drives the lapse.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum BssnSlicing {
    /// The lapse is **prescribed**: `d_t alpha = 0`, and it keeps whatever value
    /// the initial data gave it. The default, so enabling a live gauge is always
    /// an explicit choice.
    #[default]
    Prescribed,
    /// 1+log slicing, `d_t alpha = -2 alpha K`, with zero shift.
    ///
    /// Linearised about flat space the lapse then obeys a wave equation with
    /// characteristic speed `sqrt(2)`.
    OnePlusLog,
}

/// The semidiscrete BSSN system `dY/dt = F(t, Y)` on a periodic grid.
///
/// Implements [`System`] so that time integration reuses `scirust_sim`'s
/// existing fixed-step RK4. **No integrator is written here.**
#[derive(Debug, Clone, Copy)]
pub struct BssnGridSystem {
    grid: UniformGrid1d,
    gauge: BssnGauge,
    sources: AdmSources,
    dissipation: f64,
    slicing: BssnSlicing,
}

impl BssnGridSystem {
    /// Build a system on `grid` with a prescribed `gauge` and matter `sources`.
    #[must_use]
    pub const fn new(grid: UniformGrid1d, gauge: BssnGauge, sources: AdmSources) -> Self {
        Self {
            grid,
            gauge,
            sources,
            dissipation: 0.0,
            slicing: BssnSlicing::Prescribed,
        }
    }

    /// Choose the slicing condition. Defaults to [`BssnSlicing::Prescribed`], so
    /// a live gauge is never enabled implicitly.
    #[must_use]
    pub const fn with_slicing(mut self, slicing: BssnSlicing) -> Self {
        self.slicing = slicing;
        self
    }

    /// The slicing condition in force.
    #[must_use]
    pub const fn slicing(&self) -> BssnSlicing {
        self.slicing
    }

    /// Enable explicit Kreiss-Oliger dissipation with coefficient `sigma`.
    ///
    /// **Disabled by default** (`sigma = 0`), and never enabled implicitly. See
    /// [`kreiss_oliger_term`] for the formula and
    /// `docs/LAYER_3_BSSN_PERIODIC_1D.md` for the measured effect on accuracy
    /// and on the physical wave amplitude.
    ///
    /// This is **not** constraint damping: it acts on the evolved variables'
    /// high-frequency grid content, not on the constraint residuals, and it does
    /// not drive a solution toward the constraint surface.
    #[must_use]
    pub const fn with_dissipation(mut self, sigma: f64) -> Self {
        self.dissipation = sigma;
        self
    }

    /// The Kreiss-Oliger coefficient (`0` when dissipation is off).
    #[must_use]
    pub const fn dissipation(&self) -> f64 {
        self.dissipation
    }

    /// Vacuum with the synchronous gauge (`alpha = 1`, `beta^i = 0`) — the
    /// default validation configuration.
    #[must_use]
    pub fn vacuum(grid: UniformGrid1d) -> Self {
        Self::new(grid, BssnGauge::SYNCHRONOUS, AdmSources::VACUUM)
    }

    /// The grid.
    #[must_use]
    pub const fn grid(&self) -> &UniformGrid1d {
        &self.grid
    }

    /// The prescribed gauge.
    #[must_use]
    pub const fn gauge(&self) -> &BssnGauge {
        &self.gauge
    }

    /// Evaluate the right-hand side, surfacing a typed error.
    pub fn right_hand_side(
        &self,
        time: f64,
        state: &[f64],
        out: &mut [f64],
    ) -> Result<(), BssnGridError> {
        let view = BssnGridView::new(self.grid, state)?;
        bssn_grid_rhs_with_slicing(&view, &self.gauge, &self.sources, self.slicing, time, out)?;
        if self.dissipation != 0.0
        {
            apply_kreiss_oliger(&self.grid, state, self.dissipation, out);
        }
        Ok(())
    }
}

impl System for BssnGridSystem {
    fn dim(&self) -> usize {
        COMPONENTS_PER_POINT * self.grid.points()
    }

    fn derivatives(&self, time: f64, state: &[f64], derivative: &mut [f64]) {
        if self.right_hand_side(time, state, derivative).is_err()
        {
            // `System::derivatives` cannot report an error. Emitting a
            // non-finite derivative makes the integrator's own non-finite guard
            // stop the run, which `evolve_bssn_grid` then surfaces as a typed
            // `Integration` error rather than a silent wrong answer. Nothing is
            // replaced by zero.
            for value in derivative.iter_mut()
            {
                *value = f64::NAN;
            }
        }
    }
}

/// One sample of an evolution.
#[derive(Debug, Clone, PartialEq)]
pub struct BssnGridSample {
    /// The coordinate time.
    pub time: f64,
    /// The grid state at that time.
    pub state: BssnGridState,
}

/// Integrate a BSSN grid state from `t0` to `t_end` with fixed step `step`.
///
/// Reuses `scirust_sim::simulate` (fixed-step RK4). Free evolution: no
/// constraint is enforced, damped, or projected at any stage.
pub fn evolve_bssn_grid(
    system: &BssnGridSystem,
    initial: &BssnGridState,
    t0: f64,
    t_end: f64,
    step: f64,
) -> Result<Vec<BssnGridSample>, BssnGridError> {
    if !step.is_finite() || step <= 0.0
    {
        return Err(BssnGridError::InvalidRequest {
            parameter: "step",
            value: step,
        });
    }
    if !t0.is_finite()
    {
        return Err(BssnGridError::InvalidRequest {
            parameter: "t0",
            value: t0,
        });
    }
    if !t_end.is_finite() || t_end < t0
    {
        return Err(BssnGridError::InvalidRequest {
            parameter: "t_end",
            value: t_end,
        });
    }
    if initial.grid != system.grid
    {
        return Err(BssnGridError::FlatStateLength {
            expected: COMPONENTS_PER_POINT * system.grid.points(),
            actual: initial.components.len(),
        });
    }
    initial.validate(t0)?;

    let trajectory = simulate(system, initial.as_slice(), t0, t_end, step)
        .map_err(BssnGridError::Integration)?;

    let mut samples = Vec::with_capacity(trajectory.t.len());
    for (time, flat) in trajectory.t.iter().zip(trajectory.y.iter())
    {
        let state = BssnGridState::from_flat(system.grid, flat.clone())?;
        samples.push(BssnGridSample { time: *time, state });
    }
    Ok(samples)
}

// ---------------------------------------------------------------------------
// Constraint monitoring
// ---------------------------------------------------------------------------

/// Deterministic grid reductions of every monitored constraint.
///
/// Each constraint keeps its own reduction. Nothing is blended into a single
/// score, and nothing here enforces anything.
#[derive(Debug, Clone, PartialEq)]
pub struct BssnGridConstraints {
    /// `det gammatilde - 1`.
    pub determinant: GridReduction,
    /// `gammatilde^{ij} Atilde_ij`.
    pub trace_free: GridReduction,
    /// The largest `|Gammatilde^i_stored - Gammatilde^i_reconstructed|`.
    pub connection: GridReduction,
    /// The Hamiltonian constraint residual.
    pub hamiltonian: GridReduction,
    /// The three momentum-constraint components.
    pub momentum: [GridReduction; 3],
}

/// Evaluate every monitored constraint over the grid.
///
/// The BSSN algebraic constraints are read from the stored state; the ADM
/// physical constraints are evaluated from the reconstructed fields through the
/// Layer 3.1 evaluators, so no constraint formula is duplicated.
pub fn bssn_grid_constraints(
    view: &BssnGridView<'_>,
    sources: &AdmSources,
) -> Result<BssnGridConstraints, BssnGridError> {
    let points = view.grid.points();
    let settings = view.settings();
    let metric = view.metric();
    let curvature = view.curvature();

    let mut determinant = vec![0.0_f64; points];
    let mut trace_free = vec![0.0_f64; points];
    let mut connection = vec![0.0_f64; points];
    let mut hamiltonian = vec![0.0_f64; points];
    let mut momentum = [
        vec![0.0_f64; points],
        vec![0.0_f64; points],
        vec![0.0_f64; points],
    ];

    for index in 0..points
    {
        let at = [view.grid.coordinate(index), 0.0, 0.0];

        // The algebraic constraints are properties of the **stored** state, so
        // they must be read from it. Re-converting the reconstructed ADM pair
        // would silently satisfy them by construction and report zero for a
        // genuinely violated state -- the exact failure this monitor exists to
        // catch.
        //
        // `adm_to_bssn` is used only to obtain the `Gammatilde^i` implied by the
        // metric, which is what it seeds its own connection with; that is the
        // reconstruction the stored connection is compared against.
        let reconstructed = adm_to_bssn(&metric, &curvature, &at, &settings)
            .map_err(|source| BssnGridError::Pointwise {
                time: 0.0,
                index,
                stage: "adm_to_bssn",
                source,
            })?
            .conformal_connection;
        let algebraic = bssn_algebraic_constraints(&view.state_at(index), &reconstructed);
        determinant[index] = algebraic.determinant_residual;
        trace_free[index] = algebraic.trace_residual;
        // The reduction takes the worst of the three connection components so
        // that no direction can hide behind the others.
        connection[index] = algebraic
            .connection_residual
            .iter()
            .fold(0.0_f64, |worst, residual| worst.max(residual.abs()));

        let hamiltonian_value = hamiltonian_constraint(
            &metric, &curvature, &at, sources, &settings,
        )
        .map_err(|source| BssnGridError::Constraint {
            index,
            constraint: "hamiltonian",
            source,
        })?;
        hamiltonian[index] = hamiltonian_value.signed_residual;

        let momentum_value = momentum_constraint(&metric, &curvature, &at, sources, &settings)
            .map_err(|source| BssnGridError::Constraint {
                index,
                constraint: "momentum",
                source,
            })?;
        // Indexed rather than iterated: `momentum` and `residual` are distinct
        // collections walked together by tensor component, which is exactly the
        // pairing the index expresses. Matches the crate's existing convention.
        #[allow(clippy::needless_range_loop)]
        for component in 0..3
        {
            momentum[component][index] = momentum_value.residual[component];
        }
    }

    Ok(BssnGridConstraints {
        determinant: GridReduction::of(&determinant),
        trace_free: GridReduction::of(&trace_free),
        connection: GridReduction::of(&connection),
        hamiltonian: GridReduction::of(&hamiltonian),
        momentum: [
            GridReduction::of(&momentum[0]),
            GridReduction::of(&momentum[1]),
            GridReduction::of(&momentum[2]),
        ],
    })
}

/// The pointwise mismatch between the conformal Ricci decomposition and the
/// physical Ricci tensor computed independently from the reconstructed metric.
#[derive(Debug, Clone, PartialEq)]
pub struct RicciDecompositionReport {
    /// The reduction of `max_ij |Rtilde_ij + Rphi_ij - R_ij|` over the grid.
    pub mismatch: GridReduction,
    /// The largest `|Rtilde_ij|` anywhere on the grid.
    pub conformal_scale: f64,
    /// The largest `|Rphi_ij|` anywhere on the grid.
    pub conformal_factor_scale: f64,
    /// The largest `|R_ij|` anywhere on the grid.
    pub physical_scale: f64,
}

/// Cross-check the conformal Ricci decomposition against Layer 1's independently
/// computed physical Ricci tensor at every grid point.
///
/// The scales are reported alongside the mismatch so that agreement cannot be
/// claimed on the strength of both sides being near zero.
pub fn bssn_grid_ricci_report(
    view: &BssnGridView<'_>,
) -> Result<RicciDecompositionReport, BssnGridError> {
    let points = view.grid.points();
    let settings = view.settings();
    let metric = view.metric();

    let mut mismatch = vec![0.0_f64; points];
    let mut conformal_scale = 0.0_f64;
    let mut conformal_factor_scale = 0.0_f64;
    let mut physical_scale = 0.0_f64;

    // Indexed rather than iterated: the loop walks grid points and writes into
    // `mismatch` by the same index, matching the crate's existing convention.
    #[allow(clippy::needless_range_loop)]
    for index in 0..points
    {
        let at = [view.grid.coordinate(index), 0.0, 0.0];
        let decomposition = conformal_ricci(&metric, &at, &settings).map_err(|source| {
            BssnGridError::Pointwise {
                time: 0.0,
                index,
                stage: "conformal_ricci",
                source,
            }
        })?;
        let physical = crate::ricci_tensor_from_metric::<_, 3>(
            &metric,
            &at,
            settings.spatial_step,
            settings.metric_step,
        )
        .map_err(|_| BssnGridError::Pointwise {
            time: 0.0,
            index,
            stage: "ricci_tensor_from_metric",
            source: BssnError::UnavailableDerivative {
                quantity: "physical_ricci",
            },
        })?;

        let mut worst = 0.0_f64;
        // Tensor-component loops index several distinct rank-2 arrays together;
        // iterating one of them would obscure that pairing.
        #[allow(clippy::needless_range_loop)]
        for i in 0..3
        {
            for j in 0..3
            {
                worst = worst.max((decomposition.total[i][j] - physical[i][j]).abs());
                conformal_scale = conformal_scale.max(decomposition.conformal[i][j].abs());
                conformal_factor_scale =
                    conformal_factor_scale.max(decomposition.conformal_factor_part[i][j].abs());
                physical_scale = physical_scale.max(physical[i][j].abs());
            }
        }
        mismatch[index] = worst;
    }

    Ok(RicciDecompositionReport {
        mismatch: GridReduction::of(&mismatch),
        conformal_scale,
        conformal_factor_scale,
        physical_scale,
    })
}

// ---------------------------------------------------------------------------
// Explicit, opt-in projection
// ---------------------------------------------------------------------------

/// The result of projecting every point of a grid onto the algebraic constraints.
#[derive(Debug, Clone, PartialEq)]
pub struct GridProjectionReport {
    /// The worst residual before projection.
    pub residual_before: f64,
    /// The worst residual after projection.
    pub residual_after: f64,
    /// The largest correction applied at any point.
    pub correction_magnitude: f64,
}

fn accumulate(report: &mut Option<GridProjectionReport>, projection: &BssnProjection) {
    match report
    {
        Some(existing) =>
        {
            existing.residual_before = existing
                .residual_before
                .max(projection.residual_before.abs());
            existing.residual_after = existing.residual_after.max(projection.residual_after.abs());
            existing.correction_magnitude = existing
                .correction_magnitude
                .max(projection.correction_magnitude);
        },
        None =>
        {
            *report = Some(GridProjectionReport {
                residual_before: projection.residual_before.abs(),
                residual_after: projection.residual_after.abs(),
                correction_magnitude: projection.correction_magnitude,
            });
        },
    }
}

/// Project every point onto `det gammatilde = 1`.
///
/// **Explicit and opt-in.** Never invoked by [`evolve_bssn_grid`], and never
/// applied between RK4 stages.
pub fn project_grid_unit_determinant(
    state: &mut BssnGridState,
) -> Result<GridProjectionReport, BssnGridError> {
    let mut report = None;
    for index in 0..state.grid.points()
    {
        let mut local = state.state_at(index);
        let projection = project_unit_determinant(&mut local)
            .map_err(|source| BssnGridError::Projection { index, source })?;
        state.set_state_at(index, &local);
        accumulate(&mut report, &projection);
    }
    report.ok_or(BssnGridError::FlatStateLength {
        expected: COMPONENTS_PER_POINT * state.grid.points(),
        actual: 0,
    })
}

/// Project every point onto `gammatilde^{ij} Atilde_ij = 0`.
///
/// **Explicit and opt-in.** Never invoked by [`evolve_bssn_grid`], and never
/// applied between RK4 stages.
pub fn project_grid_trace_free(
    state: &mut BssnGridState,
) -> Result<GridProjectionReport, BssnGridError> {
    let mut report = None;
    for index in 0..state.grid.points()
    {
        let mut local = state.state_at(index);
        let projection = project_trace_free(&mut local)
            .map_err(|source| BssnGridError::Projection { index, source })?;
        state.set_state_at(index, &local);
        accumulate(&mut report, &projection);
    }
    report.ok_or(BssnGridError::FlatStateLength {
        expected: COMPONENTS_PER_POINT * state.grid.points(),
        actual: 0,
    })
}

// ---------------------------------------------------------------------------
// Optional explicit dissipation (off by default)
// ---------------------------------------------------------------------------

/// The fourth-order Kreiss-Oliger dissipation term for one component array.
///
/// ```text
/// Q f_i = -(sigma / 16 dx) ( f_{i+2} - 4 f_{i+1} + 6 f_i - 4 f_{i-1} + f_{i-2} )
/// ```
///
/// The bracket's Fourier symbol is `4 (1 - cos theta)^2`, which is zero for a
/// constant field and maximal (`16`) at the Nyquist mode `theta = pi`, so the
/// term damps the highest grid frequency at rate `sigma / dx` while leaving
/// smooth content almost untouched: for smooth `f` the bracket is `O(dx^4)`, so
/// `Q f = O(dx^3)` and the scheme's second-order accuracy is preserved.
///
/// This is **numerical dissipation on the evolved variables**. It is *not*
/// constraint damping — it neither targets nor reduces the constraint residuals
/// by construction, and any effect it has on them is indirect.
#[must_use]
pub fn kreiss_oliger_term(grid: &UniformGrid1d, samples: &[f64], index: usize, sigma: f64) -> f64 {
    let bracket = samples[grid.offset(index, 2)] - 4.0 * samples[grid.offset(index, 1)]
        + 6.0 * samples[grid.wrap_usize(index)]
        - 4.0 * samples[grid.offset(index, -1)]
        + samples[grid.offset(index, -2)];
    -sigma * bracket / (16.0 * grid.spacing())
}

/// Add Kreiss-Oliger dissipation to every evolved component of `out`.
fn apply_kreiss_oliger(grid: &UniformGrid1d, state: &[f64], sigma: f64, out: &mut [f64]) {
    let points = grid.points();
    for component in 0..COMPONENTS_PER_POINT
    {
        let base = component * points;
        let samples = &state[base..base + points];
        for index in 0..points
        {
            out[base + index] += kreiss_oliger_term(grid, samples, index, sigma);
        }
    }
}

// ---------------------------------------------------------------------------
// Linearized transverse-traceless wave (validation oracle)
// ---------------------------------------------------------------------------

/// A weak transverse-traceless gravitational wave propagating along `x`.
///
/// ```text
/// h_yy = +A sin(k (x - t)),    h_zz = -A sin(k (x - t)),
/// gamma_ij = diag(1, 1 + h_yy, 1 + h_zz),
/// K_ij     = -(1/2) d_t gamma_ij      (unit lapse, zero shift)
/// ```
///
/// The sign of `K_ij` follows this repository's convention
/// `K_ij = -1/(2 N) ( d_t gamma_ij - D_i N_j - D_j N_i )`, not a textbook's.
///
/// This is an exact solution of **linearized** general relativity in this gauge.
/// The evolution code solves the **full nonlinear** BSSN system, so a measured
/// difference contains both `O(A^2)` nonlinearity and `O(dx^2)` truncation. The
/// amplitude must be chosen small enough that truncation dominates before any
/// convergence claim is made. This is **not** a nonlinear-wave oracle.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TransverseTracelessWave {
    /// The dimensionless strain amplitude `A`.
    pub amplitude: f64,
    /// The wave number `k`.
    pub wave_number: f64,
    /// The coordinate time at which the fields are evaluated.
    pub time: f64,
}

impl TransverseTracelessWave {
    /// Build a wave sampled at coordinate time `time`.
    #[must_use]
    pub const fn new(amplitude: f64, wave_number: f64, time: f64) -> Self {
        Self {
            amplitude,
            wave_number,
            time,
        }
    }

    /// The phase `k (x - t)`.
    #[must_use]
    pub fn phase(&self, x: f64) -> f64 {
        self.wave_number * (x - self.time)
    }

    /// The strain `h_yy = A sin(k(x - t))` (and `h_zz = -h_yy`).
    #[must_use]
    pub fn strain(&self, x: f64) -> f64 {
        self.amplitude * self.phase(x).sin()
    }

    /// The spatial metric `gamma_ij`.
    #[must_use]
    pub fn spatial_metric(&self, x: f64) -> [[f64; 3]; 3] {
        let h = self.strain(x);
        [[1.0, 0.0, 0.0], [0.0, 1.0 + h, 0.0], [0.0, 0.0, 1.0 - h]]
    }

    /// The extrinsic curvature `K_ij = -(1/2) d_t gamma_ij`.
    #[must_use]
    pub fn extrinsic_curvature(&self, x: f64) -> [[f64; 3]; 3] {
        // d_t h = -A k cos(k(x - t)), so K_yy = +(A k / 2) cos(k(x - t)).
        let half_rate = 0.5 * self.amplitude * self.wave_number * self.phase(x).cos();
        [
            [0.0, 0.0, 0.0],
            [0.0, half_rate, 0.0],
            [0.0, 0.0, -half_rate],
        ]
    }

    /// The metric as a coordinate-sampled field.
    #[must_use]
    pub const fn metric_field(&self) -> WaveMetric<'_> {
        WaveMetric { wave: self }
    }

    /// The extrinsic curvature as a coordinate-sampled field.
    #[must_use]
    pub const fn curvature_field(&self) -> WaveCurvature<'_> {
        WaveCurvature { wave: self }
    }
}

/// The wave's spatial metric as a [`Metric<3>`].
#[derive(Debug, Clone, Copy)]
pub struct WaveMetric<'a> {
    wave: &'a TransverseTracelessWave,
}

impl Metric<3> for WaveMetric<'_> {
    fn components(&self, coordinates: &[f64; 3]) -> [[f64; 3]; 3] {
        self.wave.spatial_metric(coordinates[0])
    }
}

/// The wave's extrinsic curvature as a [`SpatialTensorField`].
#[derive(Debug, Clone, Copy)]
pub struct WaveCurvature<'a> {
    wave: &'a TransverseTracelessWave,
}

impl SpatialTensorField for WaveCurvature<'_> {
    fn components(&self, coordinates: &[f64; 3]) -> [[f64; 3]; 3] {
        self.wave.extrinsic_curvature(coordinates[0])
    }
}

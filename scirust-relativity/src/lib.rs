//! # `scirust-relativity` — deterministic differential geometry
//!
//! This crate provides a minimal, auditable foundation for spacetime geometry
//! and geodesic simulation:
//!
//! - covariant metric tensors;
//! - deterministic metric inversion;
//! - numerical Levi-Civita Christoffel symbols;
//! - four-dimensional Minkowski, Schwarzschild, Reissner-Nordström, Kerr,
//!   de Sitter, anti-de Sitter, and spatially flat FLRW spacetimes;
//! - Riemann, Ricci, Einstein, and Kretschmann curvature tensors from any
//!   metric-and-connection background, validated against exact analytic
//!   oracles (see [`CurvatureTensors`]);
//! - a reusable parallel-transport engine for vectors ([`transport_along_segment`],
//!   [`transport_along_polyline`], [`holonomy_defect`]), covectors, and rank-2
//!   covariant tensors ([`transport_covector_along_segment`],
//!   [`transport_covariant_tensor_along_segment`]), validated by metric
//!   compatibility and the holonomy/curvature identity;
//! - geodesic equations compatible with `scirust-sim`, geodesic-deviation
//!   (Jacobi) fields ([`integrate_geodesic_deviation`]) validated against the
//!   separation of nearby geodesics, and geodesic exponential/logarithm maps
//!   ([`geodesic_exponential`], [`geodesic_logarithm`]);
//! - local orthonormal frames (tetrads) for timelike observers
//!   ([`orthonormal_tetrad`]), built by metric Gram-Schmidt and validated by
//!   orthonormality, completeness, and preservation under parallel transport;
//! - Synge's world function and its bitensors — the biscalar and its
//!   first-derivative gradients ([`world_function`],
//!   [`world_function_with_gradients`]) and the van Vleck-Morette determinant
//!   ([`van_vleck_determinant`]) — built on the geodesic logarithm map and a
//!   deterministic [`determinant`], validated by flat exactness, the fundamental
//!   identity, base/field symmetry, and the known maximally-symmetric
//!   coincidence expansion;
//! - linearized gravity ([`LinearizedField`]) — the field equations to first
//!   order in a metric perturbation about Minkowski (the opening slice of the
//!   Layer 2 Covariant Gravity Workbench), validated by weak-field-Schwarzschild
//!   vacuum, the Newtonian Poisson limit, gauge invariance, and an `O(h^2)`
//!   cross-check against the nonlinear curvature.
//!
//! The crate does not assume that fractional calculus modifies general
//! relativity. Such models, if added later, must be exposed explicitly as
//! experimental constitutive or non-local extensions.
//!
//! ## Example
//!
//! ```
//! use scirust_relativity::{GeodesicSystem, Minkowski};
//! use scirust_sim::simulate;
//!
//! let system = GeodesicSystem::<_, 4>::new(Minkowski);
//! let initial = [0.0, 0.0, 0.0, 0.0, 1.0, 0.25, 0.0, 0.0];
//!
//! let trajectory = simulate(&system, &initial, 0.0, 2.0, 0.01)
//!     .expect("Minkowski geodesic integrates");
//!
//! let final_state = trajectory.last_state().expect("non-empty trajectory");
//! assert!((final_state[1] - 0.5).abs() < 1.0e-12);
//! ```

#![forbid(unsafe_code)]
#![deny(missing_docs)]

mod connection;
mod covariant_transport;
mod curvature;
mod de_sitter;
mod error;
mod exponential_map;
mod flrw;
mod geodesic;
mod geodesic_deviation;
mod isotropic_schwarzschild;
mod kerr;
mod linearized;
mod metric;
mod minkowski;
mod minkowski_spherical;
mod parallel_transport;
mod reissner_nordstrom;
mod schwarzschild;
mod static_spherical;
mod synge;
mod tetrad;

pub use connection::{Connection, numerical_christoffel};
pub use covariant_transport::{
    transport_covariant_tensor_along_polyline, transport_covariant_tensor_along_segment,
    transport_covector_along_polyline, transport_covector_along_segment,
};
pub use curvature::CurvatureTensors;
pub use de_sitter::{AntiDeSitter, DeSitter};
pub use error::RelativityError;
pub use exponential_map::{geodesic_exponential, geodesic_logarithm};
pub use flrw::{ExponentialScaleFactor, Flrw, PowerLawScaleFactor, ScaleFactor};
pub use geodesic::GeodesicSystem;
pub use geodesic_deviation::{JacobiSample, integrate_geodesic_deviation};
pub use isotropic_schwarzschild::IsotropicSchwarzschild;
pub use kerr::Kerr;
pub use linearized::LinearizedField;
pub use metric::{Metric, determinant, invert_metric, metric_norm};
pub use minkowski::Minkowski;
pub use minkowski_spherical::MinkowskiSpherical;
pub use parallel_transport::{holonomy_defect, transport_along_polyline, transport_along_segment};
pub use reissner_nordstrom::ReissnerNordstrom;
pub use schwarzschild::Schwarzschild;
pub use synge::{
    WorldFunction, WorldFunctionSettings, van_vleck_determinant, world_function,
    world_function_with_gradients,
};
pub use tetrad::{OrthonormalTetrad, orthonormal_tetrad};

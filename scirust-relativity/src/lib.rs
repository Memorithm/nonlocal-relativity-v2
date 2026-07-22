//! # `scirust-relativity` — deterministic differential geometry
//!
//! This crate provides a minimal, auditable foundation for spacetime geometry
//! and geodesic simulation:
//!
//! - covariant metric tensors;
//! - deterministic metric inversion;
//! - numerical Levi-Civita Christoffel symbols;
//! - four-dimensional Minkowski, Schwarzschild, Reissner-Nordström, Kerr,
//!   de Sitter, and anti-de Sitter spacetimes;
//! - Riemann, Ricci, Einstein, and Kretschmann curvature tensors from any
//!   metric-and-connection background, validated against exact analytic
//!   oracles (see [`CurvatureTensors`]);
//! - geodesic equations compatible with `scirust-sim`.
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
mod curvature;
mod de_sitter;
mod error;
mod geodesic;
mod isotropic_schwarzschild;
mod kerr;
mod metric;
mod minkowski;
mod minkowski_spherical;
mod reissner_nordstrom;
mod schwarzschild;
mod static_spherical;

pub use connection::{Connection, numerical_christoffel};
pub use curvature::CurvatureTensors;
pub use de_sitter::{AntiDeSitter, DeSitter};
pub use error::RelativityError;
pub use geodesic::GeodesicSystem;
pub use isotropic_schwarzschild::IsotropicSchwarzschild;
pub use kerr::Kerr;
pub use metric::{Metric, invert_metric, metric_norm};
pub use minkowski::Minkowski;
pub use minkowski_spherical::MinkowskiSpherical;
pub use reissner_nordstrom::ReissnerNordstrom;
pub use schwarzschild::Schwarzschild;

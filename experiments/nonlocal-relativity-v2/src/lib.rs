//! Shared helpers for the deterministic nonlocal-relativity v2 experiments.
//!
//! Every experiment binary in this crate is deterministic (no RNG, no wall
//! clock), prints a `#`-prefixed metadata/units header followed by CSV rows,
//! validates that every emitted number is finite, and closes with a short,
//! non-overclaiming interpretation. These are numerical experiments on a fixed
//! phenomenological model; none of them is a physical validation.

#![forbid(unsafe_code)]

use scirust_nonlocal_relativity::WorldlineState;
use std::f64::consts::FRAC_PI_2;

/// Circular equatorial geodesic four-velocity initial state for Schwarzschild
/// mass `mass` at areal radius `radius` (requires `radius > 3 * mass`).
#[must_use]
pub fn circular_schwarzschild_state(mass: f64, radius: f64) -> WorldlineState<4> {
    let denominator = (1.0 - 3.0 * mass / radius).sqrt();
    let u_t = 1.0 / denominator;
    let u_phi = (mass / (radius * radius * radius)).sqrt() / denominator;
    WorldlineState::new([0.0, radius, FRAC_PI_2, 0.0], [u_t, 0.0, 0.0, u_phi])
}

/// Deterministic Euclidean distance between two 4-vectors (a coordinate-chart
/// diagnostic, not a spacetime interval).
#[must_use]
pub fn euclidean_distance(left: &[f64; 4], right: &[f64; 4]) -> f64 {
    (0..4)
        .map(|component| (left[component] - right[component]).powi(2))
        .sum::<f64>()
        .sqrt()
}

/// Commit hash for provenance, read from the `NLR_EXPERIMENT_COMMIT`
/// environment variable when set (e.g. `git rev-parse HEAD`), otherwise
/// `"unset"`. This is metadata only and never enters a CSV data row, so the
/// numeric output stays reproducible regardless of it.
#[must_use]
pub fn commit_hash() -> String {
    std::env::var("NLR_EXPERIMENT_COMMIT").unwrap_or_else(|_| "unset".to_string())
}

/// Assert every value in `values` is finite, returning an error string naming
/// the offending label otherwise. Experiment binaries treat a non-finite
/// emitted value as a hard failure.
pub fn require_finite(values: &[(&str, f64)]) -> Result<(), String> {
    for (label, value) in values
    {
        if !value.is_finite()
        {
            return Err(format!("non-finite value for '{label}': {value}"));
        }
    }
    Ok(())
}

/// Print the shared metadata header lines (as `#` comments) common to every
/// experiment: the experiment name, geometric-unit convention, and commit.
pub fn print_common_header(experiment: &str) {
    println!("# experiment: {experiment}");
    println!("# layer: scirust-nonlocal-relativity (experimental, phenomenological)");
    println!(
        "# units: geometric G = c = 1; lengths, times, and the affine parameter in mass units M"
    );
    println!("# determinism: no RNG, no wall clock; identical inputs give identical output");
    println!("# commit: {}", commit_hash());
    println!("# NOTE: numerical experiment on a fixed model; not a physical validation.");
}

//! PPN (gamma, beta) extraction checks: synthetic oracles, controlled
//! contamination, isotropic Schwarzschild convergence, and areal-coordinate
//! rejection (the second Layer 2 slice; see `docs/LAYER_2_PPN.md`).
//!
//! For each metric and each (window, sample count, fit order) the experiment
//! runs the deterministic asymptotic extractor and reports the estimate, its
//! absolute error against the known oracle, the *estimated* numerical
//! uncertainty (a sensitivity, not a bound), the fit residuals, and the fit
//! conditioning. The areal-Schwarzschild rows show the coordinate-validity
//! rejection.
//!
//! Established general relativity only. The extraction is a numerical
//! approximation; the isotropic-Schwarzschild oracle (`gamma = beta = 1`)
//! validates the implementation, not an alternative theory.

#![forbid(unsafe_code)]

use nonlocal_relativity_experiments::{print_experiment_header, require_finite};
use scirust_relativity::ppn::{
    IsotropicChartAdapter, PpnDomain, PpnError, StaticIsotropicMetric, SyntheticPpnMetric,
    extract_ppn,
};
use scirust_relativity::{IsotropicSchwarzschild, Schwarzschild};

const MASS: f64 = 1.0;

/// Short, comma-free status token for a rejected extraction.
fn rejection_token(error: &PpnError) -> &'static str {
    match error
    {
        PpnError::NonIsotropicCoordinates { .. } => "rejected_non_isotropic",
        PpnError::CompactnessOutOfRange { .. } => "rejected_compactness",
        PpnError::NonAsymptoticallyFlat { .. } => "rejected_not_flat",
        _ => "rejected_other",
    }
}

#[allow(clippy::too_many_arguments)]
fn emit_row<M: StaticIsotropicMetric>(
    name: &str,
    oracle_gamma: f64,
    oracle_beta: f64,
    metric: &M,
    compactness_max: f64,
    sample_count: usize,
    fit_order: usize,
) -> Result<(), String> {
    let compactness_min = compactness_max * 0.1;
    let radius_min = MASS / compactness_max;
    let radius_max = MASS / compactness_min;
    let domain = PpnDomain::uniform_compactness(compactness_min, compactness_max, sample_count);

    match extract_ppn(metric, &domain, fit_order)
    {
        Ok(estimate) =>
        {
            let gamma_error = (estimate.gamma.asymptotic_value - oracle_gamma).abs();
            let beta_error = (estimate.beta.asymptotic_value - oracle_beta).abs();
            require_finite(&[
                ("gamma", estimate.gamma.asymptotic_value),
                ("beta", estimate.beta.asymptotic_value),
                ("gamma_uncertainty", estimate.gamma.estimated_uncertainty),
                ("beta_uncertainty", estimate.beta.estimated_uncertainty),
                ("gamma_residual", estimate.gamma.fit_residual),
                ("beta_residual", estimate.beta.fit_residual),
                ("conditioning", estimate.gamma.conditioning),
            ])?;
            println!(
                "{name},{oracle_gamma},{oracle_beta},{radius_min:.3},{radius_max:.3},\
                 {:.4e},{:.4e},{sample_count},{fit_order},{:.10},{gamma_error:.3e},{:.3e},\
                 {:.10},{beta_error:.3e},{:.3e},{:.3e},{:.3e},{:.3e},ok",
                estimate.compactness_min,
                estimate.compactness_max,
                estimate.gamma.asymptotic_value,
                estimate.gamma.estimated_uncertainty,
                estimate.beta.asymptotic_value,
                estimate.beta.estimated_uncertainty,
                estimate.gamma.fit_residual,
                estimate.beta.fit_residual,
                estimate.gamma.conditioning,
            );
        },
        Err(error) =>
        {
            println!(
                "{name},{oracle_gamma},{oracle_beta},{radius_min:.3},{radius_max:.3},\
                 {compactness_min:.4e},{compactness_max:.4e},{sample_count},{fit_order},\
                 nan,nan,nan,nan,nan,nan,nan,nan,nan,{}",
                rejection_token(&error)
            );
        },
    }
    Ok(())
}

fn main() -> Result<(), String> {
    print_experiment_header(
        "PPN parameter (gamma, beta) extraction",
        "scirust-relativity geometry core (established general relativity)",
        "asymptotic extraction of gamma, beta from isotropic weak-field metrics; numerical, not a bound.",
    );

    println!(
        "metric,oracle_gamma,oracle_beta,radius_min,radius_max,compactness_min,compactness_max,\
         sample_count,fit_order,gamma_estimate,gamma_abs_error,gamma_estimated_uncertainty,\
         beta_estimate,beta_abs_error,beta_estimated_uncertainty,gamma_fit_residual,\
         beta_fit_residual,conditioning_indicator,status"
    );

    let exact = SyntheticPpnMetric::exact(MASS, 1.0, 1.0);
    let contaminated = SyntheticPpnMetric::contaminated(MASS, 1.0, 1.0, 3.0, -5.0, 2.0, 4.0);
    let isotropic =
        IsotropicSchwarzschild::try_new(MASS).ok_or("invalid isotropic Schwarzschild")?;
    let isotropic_adapter =
        IsotropicChartAdapter::new(&isotropic, MASS).map_err(|e| e.to_string())?;
    let areal = Schwarzschild::try_new(MASS).ok_or("invalid Schwarzschild")?;
    let areal_adapter = IsotropicChartAdapter::new(&areal, MASS).map_err(|e| e.to_string())?;

    for &compactness_max in &[0.1, 0.05]
    {
        for &sample_count in &[12, 24]
        {
            for &fit_order in &[2, 3]
            {
                emit_row(
                    "synthetic_exact",
                    1.0,
                    1.0,
                    &exact,
                    compactness_max,
                    sample_count,
                    fit_order,
                )?;
                emit_row(
                    "synthetic_contaminated",
                    1.0,
                    1.0,
                    &contaminated,
                    compactness_max,
                    sample_count,
                    fit_order,
                )?;
                emit_row(
                    "isotropic_schwarzschild",
                    1.0,
                    1.0,
                    &isotropic_adapter,
                    compactness_max,
                    sample_count,
                    fit_order,
                )?;
            }
        }
    }
    // Coordinate-validity: areal Schwarzschild is rejected, not silently trusted.
    emit_row("areal_schwarzschild", 1.0, 1.0, &areal_adapter, 0.05, 16, 3)?;

    println!("# interpretation: the exact synthetic metric recovers the injected pair to near");
    println!("# machine precision; contamination biases the finite-radius estimators while the");
    println!("# extrapolation recovers the injected values (and a weaker-field window or higher");
    println!(
        "# order improves it); isotropic Schwarzschild converges to gamma = beta = 1; and areal"
    );
    println!("# Schwarzschild is rejected as non-isotropic. Numerical extraction, not a bound.");
    Ok(())
}

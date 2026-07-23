//! Validation of PPN (gamma, beta) extraction against the design-note oracles
//! (`docs/LAYER_2_PPN.md`):
//!
//! 1. exact synthetic metrics recover the injected (gamma, beta);
//! 2. controlled higher-order contamination biases the finite-radius estimators,
//!    the extrapolation recovers the injected values, a weaker-field window
//!    improves the result, and the reported sensitivity tracks the error;
//! 3. exact isotropic Schwarzschild converges to gamma = beta = 1;
//! 4. areal-coordinate Schwarzschild is rejected, never silently misused.
//!
//! Plus the domain / fit / metric error paths and determinism.

use scirust_relativity::ppn::{
    IsotropicChartAdapter, PpnDomain, PpnError, SyntheticPpnMetric, extract_ppn,
};
use scirust_relativity::{IsotropicSchwarzschild, Schwarzschild};

const MASS: f64 = 1.0;

// -------------------------------------------------------------------------
// Oracle 1 — exact synthetic
// -------------------------------------------------------------------------

#[test]
fn exact_synthetic_recovers_injected_pairs() {
    for (gamma_star, beta_star) in [(1.0, 1.0), (0.5, 1.2), (1.3, 0.8), (0.9, 1.05)]
    {
        let metric = SyntheticPpnMetric::exact(MASS, gamma_star, beta_star);
        let domain = PpnDomain::uniform_compactness(0.01, 0.1, 12);
        let estimate = extract_ppn(&metric, &domain, 2).expect("exact extraction");
        assert!(
            (estimate.gamma.asymptotic_value - gamma_star).abs() < 1.0e-10,
            "gamma {} != {gamma_star}",
            estimate.gamma.asymptotic_value
        );
        assert!(
            (estimate.beta.asymptotic_value - beta_star).abs() < 1.0e-10,
            "beta {} != {beta_star}",
            estimate.beta.asymptotic_value
        );
    }
}

// -------------------------------------------------------------------------
// Oracle 2 — controlled contamination
// -------------------------------------------------------------------------

fn contaminated() -> SyntheticPpnMetric {
    // Higher-order terms a3, a4 (in g_00) and b2, b3 (in A): the effective
    // estimators become degree-2 polynomials in U, so a degree-2 fit is exact
    // while the finite-radius values are visibly biased.
    SyntheticPpnMetric::contaminated(MASS, 1.0, 1.0, 3.0, -5.0, 2.0, 4.0)
}

#[test]
fn contaminated_finite_radius_biased_extrapolation_recovers() {
    let metric = contaminated();
    let domain = PpnDomain::uniform_compactness(0.01, 0.1, 16);
    let estimate = extract_ppn(&metric, &domain, 2).expect("contaminated extraction");

    // The strongest-field effective beta is visibly off from the injected 1.
    let strongest = estimate
        .beta
        .finite_radius_values
        .iter()
        .max_by(|a, b| a.compactness.total_cmp(&b.compactness))
        .expect("samples");
    assert!(
        (strongest.value - 1.0).abs() > 0.05,
        "finite-radius beta {} should be biased",
        strongest.value
    );

    // The degree-2 extrapolation recovers the injected values to near machine
    // precision (the contamination is exactly degree-2 in U).
    assert!((estimate.gamma.asymptotic_value - 1.0).abs() < 1.0e-9);
    assert!((estimate.beta.asymptotic_value - 1.0).abs() < 1.0e-9);
}

#[test]
fn contaminated_weaker_window_reduces_error() {
    // Degree-1 leaves a residual O(U) bias, so shrinking the window reduces it.
    let metric = contaminated();
    let wide = extract_ppn(&metric, &PpnDomain::uniform_compactness(0.012, 0.12, 12), 1)
        .expect("wide")
        .gamma
        .asymptotic_value;
    let narrow = extract_ppn(&metric, &PpnDomain::uniform_compactness(0.004, 0.04, 12), 1)
        .expect("narrow")
        .gamma
        .asymptotic_value;
    assert!(
        (narrow - 1.0).abs() < (wide - 1.0).abs(),
        "narrow error {} not < wide error {}",
        (narrow - 1.0).abs(),
        (wide - 1.0).abs()
    );
}

#[test]
fn contaminated_uncertainty_tracks_error() {
    // The estimated uncertainty must not under-report the true error.
    let metric = contaminated();
    let estimate = extract_ppn(&metric, &PpnDomain::uniform_compactness(0.012, 0.12, 12), 1)
        .expect("extraction");
    let gamma_error = (estimate.gamma.asymptotic_value - 1.0).abs();
    assert!(
        gamma_error <= 1.5 * estimate.gamma.estimated_uncertainty,
        "error {gamma_error} exceeds 1.5x uncertainty {}",
        estimate.gamma.estimated_uncertainty
    );
    // And it is not absurdly loose.
    assert!(estimate.gamma.estimated_uncertainty <= 100.0 * gamma_error.max(1.0e-15));
}

// -------------------------------------------------------------------------
// Oracle 3 — isotropic Schwarzschild
// -------------------------------------------------------------------------

#[test]
fn isotropic_schwarzschild_converges_to_general_relativity() {
    let isotropic = IsotropicSchwarzschild::try_new(MASS).unwrap();
    let adapter = IsotropicChartAdapter::new(&isotropic, MASS).unwrap();
    let domain = PpnDomain::uniform_compactness(0.005, 0.05, 20);

    let low = extract_ppn(&adapter, &domain, 2).expect("degree 2");
    let high = extract_ppn(&adapter, &domain, 3).expect("degree 3");

    // Converges to gamma = beta = 1.
    assert!((high.gamma.asymptotic_value - 1.0).abs() < 1.0e-9);
    assert!((high.beta.asymptotic_value - 1.0).abs() < 1.0e-5);

    // Higher degree is at least as accurate.
    assert!(
        (high.gamma.asymptotic_value - 1.0).abs() <= (low.gamma.asymptotic_value - 1.0).abs(),
        "degree 3 gamma not better than degree 2"
    );
    assert!(
        (high.beta.asymptotic_value - 1.0).abs() <= (low.beta.asymptotic_value - 1.0).abs(),
        "degree 3 beta not better than degree 2"
    );
}

#[test]
fn isotropic_schwarzschild_stable_under_more_samples() {
    let isotropic = IsotropicSchwarzschild::try_new(MASS).unwrap();
    let adapter = IsotropicChartAdapter::new(&isotropic, MASS).unwrap();
    let coarse = extract_ppn(
        &adapter,
        &PpnDomain::uniform_compactness(0.005, 0.05, 16),
        3,
    )
    .expect("coarse");
    let fine = extract_ppn(
        &adapter,
        &PpnDomain::uniform_compactness(0.005, 0.05, 40),
        3,
    )
    .expect("fine");
    assert!((coarse.gamma.asymptotic_value - fine.gamma.asymptotic_value).abs() < 1.0e-6);
    assert!((coarse.beta.asymptotic_value - fine.beta.asymptotic_value).abs() < 1.0e-4);
}

// -------------------------------------------------------------------------
// Oracle 4 — invalid coordinates (mandatory negative test)
// -------------------------------------------------------------------------

#[test]
fn areal_schwarzschild_is_rejected() {
    let areal = Schwarzschild::try_new(MASS).unwrap();
    let adapter = IsotropicChartAdapter::new(&areal, MASS).unwrap();
    let domain = PpnDomain::uniform_compactness(0.005, 0.05, 12);
    assert!(matches!(
        extract_ppn(&adapter, &domain, 2),
        Err(PpnError::NonIsotropicCoordinates { .. })
    ));
}

// -------------------------------------------------------------------------
// Error paths
// -------------------------------------------------------------------------

#[test]
fn rejects_invalid_inputs() {
    let good = SyntheticPpnMetric::exact(MASS, 1.0, 1.0);
    let good_domain = PpnDomain::uniform_compactness(0.01, 0.1, 12);

    // Unsupported degree.
    assert!(matches!(
        extract_ppn(&good, &good_domain, 0),
        Err(PpnError::UnsupportedExtrapolationOrder { .. })
    ));
    assert!(matches!(
        extract_ppn(&good, &good_domain, 99),
        Err(PpnError::UnsupportedExtrapolationOrder { .. })
    ));

    // Invalid mass scale.
    let bad_mass = SyntheticPpnMetric::exact(-1.0, 1.0, 1.0);
    assert!(matches!(
        extract_ppn(&bad_mass, &good_domain, 2),
        Err(PpnError::InvalidMassScale(_))
    ));

    // Compactness outside the weak-field window.
    assert!(matches!(
        extract_ppn(&good, &PpnDomain::uniform_compactness(0.05, 0.5, 12), 2),
        Err(PpnError::CompactnessOutOfRange { .. })
    ));

    // Malformed radial domain.
    assert!(matches!(
        extract_ppn(&good, &PpnDomain::logarithmic_radius(50.0, 10.0, 12), 2),
        Err(PpnError::InvalidRadialDomain { .. })
    ));

    // Too few samples for the degree.
    assert!(matches!(
        extract_ppn(&good, &PpnDomain::uniform_compactness(0.01, 0.1, 3), 3),
        Err(PpnError::InsufficientSamples { .. })
    ));

    // Non-finite metric value.
    let nan_metric = SyntheticPpnMetric::exact(MASS, f64::NAN, 1.0);
    assert!(matches!(
        extract_ppn(&nan_metric, &good_domain, 2),
        Err(PpnError::NonFiniteMetricValue { .. })
    ));

    // Not a weak-field perturbation (huge gamma makes A - 1 order one).
    let strong = SyntheticPpnMetric::exact(MASS, 100.0, 1.0);
    assert!(matches!(
        extract_ppn(&strong, &good_domain, 2),
        Err(PpnError::NonAsymptoticallyFlat { .. })
    ));

    // Degenerate (all-equal) radii give a singular / ill-conditioned fit.
    let degenerate = PpnDomain::explicit_radii(vec![20.0, 20.0, 20.0, 20.0]);
    assert!(matches!(
        extract_ppn(&good, &degenerate, 2),
        Err(PpnError::SingularFit) | Err(PpnError::IllConditionedFit { .. })
    ));
}

// -------------------------------------------------------------------------
// Determinism
// -------------------------------------------------------------------------

#[test]
fn extraction_is_deterministic() {
    let isotropic = IsotropicSchwarzschild::try_new(MASS).unwrap();
    let adapter = IsotropicChartAdapter::new(&isotropic, MASS).unwrap();
    let domain = PpnDomain::uniform_compactness(0.005, 0.05, 24);
    let first = extract_ppn(&adapter, &domain, 3).unwrap();
    let second = extract_ppn(&adapter, &domain, 3).unwrap();
    assert_eq!(first, second);
}

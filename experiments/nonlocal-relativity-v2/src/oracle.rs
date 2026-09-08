//! Reusable evidence and oracle helpers for deterministic scientific tests.
//!
//! This module does not decide whether a model is physically valid. It records
//! what kind of evidence a numerical check represents and provides small,
//! deterministic primitives for comparing measured values with exact or
//! numerical references. The distinction between exact analytic evidence,
//! independent numerical references, self-convergence, regression checks, and
//! phenomenological sensitivity is explicit in the type system instead of
//! being left to free-form prose.

use core::fmt;

/// Scientific evidence category carried by an oracle check.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EvidenceClass {
    /// Exact analytic value or identity.
    ExactAnalytic,
    /// Exact closed-form reference for a deliberately limited family.
    ClosedForm,
    /// Independent numerical implementation or independently converged path.
    IndependentNumerical,
    /// Fine-grid numerical reference, not an exact solution.
    FineGridReference,
    /// Self-convergence under a controlled refinement sweep.
    SelfConvergence,
    /// Regression/compatibility check that preserves already-established behavior.
    Regression,
    /// Sensitivity measurement inside an explicitly phenomenological model.
    PhenomenologicalSensitivity,
}

impl EvidenceClass {
    /// Stable lowercase identifier for CSV/metadata output.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self
        {
            Self::ExactAnalytic => "exact_analytic",
            Self::ClosedForm => "closed_form",
            Self::IndependentNumerical => "independent_numerical",
            Self::FineGridReference => "fine_grid_reference",
            Self::SelfConvergence => "self_convergence",
            Self::Regression => "regression",
            Self::PhenomenologicalSensitivity => "phenomenological_sensitivity",
        }
    }
}

/// Typed failure constructing or evaluating an oracle diagnostic.
#[derive(Debug, Clone, PartialEq)]
pub enum OracleError {
    /// A required numeric input was NaN or infinite.
    NonFiniteInput {
        /// Name of the offending input.
        name: &'static str,
        /// Value supplied.
        value: f64,
    },
    /// A tolerance or scale required to be strictly positive was not.
    NonPositiveScale {
        /// Name of the offending scale.
        name: &'static str,
        /// Value supplied.
        value: f64,
    },
    /// An observed-order calculation needs strictly positive errors.
    NonPositiveError {
        /// Name of the offending error term.
        name: &'static str,
        /// Value supplied.
        value: f64,
    },
    /// A refinement ratio must be strictly greater than one.
    InvalidRefinementRatio(f64),
    /// A monotonicity check needs at least two samples.
    InsufficientSamples(usize),
}

impl fmt::Display for OracleError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self
        {
            Self::NonFiniteInput { name, value } =>
            {
                write!(
                    formatter,
                    "oracle input '{name}' must be finite; got {value}"
                )
            },
            Self::NonPositiveScale { name, value } => write!(
                formatter,
                "oracle scale '{name}' must be finite and strictly positive; got {value}"
            ),
            Self::NonPositiveError { name, value } => write!(
                formatter,
                "oracle error '{name}' must be finite and strictly positive; got {value}"
            ),
            Self::InvalidRefinementRatio(value) => write!(
                formatter,
                "refinement ratio must be finite and strictly greater than one; got {value}"
            ),
            Self::InsufficientSamples(count) => write!(
                formatter,
                "oracle monotonicity check needs at least two samples; got {count}"
            ),
        }
    }
}

impl std::error::Error for OracleError {}

/// Result of one scalar oracle comparison.
#[derive(Debug, Clone, PartialEq)]
pub struct OracleCheck {
    /// Human-readable quantity name.
    pub label: String,
    /// Evidence category of this comparison.
    pub evidence: EvidenceClass,
    /// Measured value.
    pub measured: f64,
    /// Reference value.
    pub reference: f64,
    /// Absolute error `|measured-reference|`.
    pub absolute_error: f64,
    /// Normalized error used by the selected check.
    pub normalized_error: f64,
    /// Acceptance threshold for `normalized_error`.
    pub tolerance: f64,
    /// Whether the check satisfies its threshold.
    pub passed: bool,
}

impl OracleCheck {
    /// Compare against a reference with an absolute tolerance.
    pub fn absolute(
        label: impl Into<String>,
        evidence: EvidenceClass,
        measured: f64,
        reference: f64,
        tolerance: f64,
    ) -> Result<Self, OracleError> {
        require_finite("measured", measured)?;
        require_finite("reference", reference)?;
        require_positive("tolerance", tolerance)?;

        let absolute_error = (measured - reference).abs();
        Ok(Self {
            label: label.into(),
            evidence,
            measured,
            reference,
            absolute_error,
            normalized_error: absolute_error,
            tolerance,
            passed: absolute_error <= tolerance,
        })
    }

    /// Compare against a reference with a relative tolerance.
    ///
    /// `scale_floor` prevents a nominally relative test from silently becoming
    /// arbitrarily strict near a zero reference. It is explicit because the
    /// correct floor is domain-specific and must never be invented internally.
    pub fn relative(
        label: impl Into<String>,
        evidence: EvidenceClass,
        measured: f64,
        reference: f64,
        relative_tolerance: f64,
        scale_floor: f64,
    ) -> Result<Self, OracleError> {
        require_finite("measured", measured)?;
        require_finite("reference", reference)?;
        require_positive("relative_tolerance", relative_tolerance)?;
        require_positive("scale_floor", scale_floor)?;

        let absolute_error = (measured - reference).abs();
        let scale = reference.abs().max(scale_floor);
        let normalized_error = absolute_error / scale;
        Ok(Self {
            label: label.into(),
            evidence,
            measured,
            reference,
            absolute_error,
            normalized_error,
            tolerance: relative_tolerance,
            passed: normalized_error <= relative_tolerance,
        })
    }
}

/// Result of a non-triviality/activation guard.
#[derive(Debug, Clone, PartialEq)]
pub struct ActivationCheck {
    /// Quantity that must be exercised by the test.
    pub label: String,
    /// Measured absolute scale.
    pub absolute_value: f64,
    /// Minimum scale required for the term to count as activated.
    pub minimum_absolute_value: f64,
    /// Whether the tested configuration genuinely activates the quantity.
    pub passed: bool,
}

/// Require a term to be measurably non-zero in a test configuration.
///
/// This prevents a test from claiming coverage of a tensor term that vanishes
/// structurally in the chosen symmetry reduction.
pub fn activation_guard(
    label: impl Into<String>,
    value: f64,
    minimum_absolute_value: f64,
) -> Result<ActivationCheck, OracleError> {
    require_finite("activation_value", value)?;
    require_positive("minimum_absolute_value", minimum_absolute_value)?;
    let absolute_value = value.abs();
    Ok(ActivationCheck {
        label: label.into(),
        absolute_value,
        minimum_absolute_value,
        passed: absolute_value >= minimum_absolute_value,
    })
}

/// Compute observed convergence order from two errors and a refinement ratio.
///
/// For `E_h ~ C h^p` and `refinement_ratio = h_coarse/h_fine`, this returns
/// `p = ln(E_coarse/E_fine)/ln(refinement_ratio)`.
pub fn observed_order(
    coarse_error: f64,
    fine_error: f64,
    refinement_ratio: f64,
) -> Result<f64, OracleError> {
    require_positive_error("coarse_error", coarse_error)?;
    require_positive_error("fine_error", fine_error)?;
    if !refinement_ratio.is_finite() || refinement_ratio <= 1.0
    {
        return Err(OracleError::InvalidRefinementRatio(refinement_ratio));
    }
    Ok((coarse_error / fine_error).ln() / refinement_ratio.ln())
}

/// Return whether a sequence decreases strictly at every refinement step.
pub fn strictly_decreasing(values: &[f64]) -> Result<bool, OracleError> {
    if values.len() < 2
    {
        return Err(OracleError::InsufficientSamples(values.len()));
    }
    for &value in values
    {
        require_finite("monotonicity_value", value)?;
    }
    Ok(values.windows(2).all(|pair| pair[1] < pair[0]))
}

fn require_finite(name: &'static str, value: f64) -> Result<(), OracleError> {
    if !value.is_finite()
    {
        return Err(OracleError::NonFiniteInput { name, value });
    }
    Ok(())
}

fn require_positive(name: &'static str, value: f64) -> Result<(), OracleError> {
    if !value.is_finite() || value <= 0.0
    {
        return Err(OracleError::NonPositiveScale { name, value });
    }
    Ok(())
}

fn require_positive_error(name: &'static str, value: f64) -> Result<(), OracleError> {
    if !value.is_finite() || value <= 0.0
    {
        return Err(OracleError::NonPositiveError { name, value });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn evidence_identifiers_are_stable() {
        assert_eq!(EvidenceClass::ExactAnalytic.as_str(), "exact_analytic");
        assert_eq!(
            EvidenceClass::IndependentNumerical.as_str(),
            "independent_numerical"
        );
        assert_eq!(
            EvidenceClass::PhenomenologicalSensitivity.as_str(),
            "phenomenological_sensitivity"
        );
    }

    #[test]
    fn absolute_and_relative_checks_report_components_separately() {
        let absolute = OracleCheck::absolute(
            "identity",
            EvidenceClass::ExactAnalytic,
            1.0 + 1.0e-10,
            1.0,
            2.0e-10,
        )
        .unwrap();
        assert!(absolute.passed);
        assert_eq!(absolute.normalized_error, absolute.absolute_error);

        let relative = OracleCheck::relative(
            "curvature",
            EvidenceClass::ClosedForm,
            2.02e-4,
            2.00e-4,
            0.02,
            1.0e-12,
        )
        .unwrap();
        assert!(relative.passed);
        assert!((relative.normalized_error - 0.01).abs() < 1.0e-12);
    }

    #[test]
    fn activation_guard_rejects_trivial_coverage() {
        let inactive = activation_guard("mixed_term", 1.0e-14, 1.0e-8).unwrap();
        let active = activation_guard("mixed_term", 2.0e-3, 1.0e-8).unwrap();
        assert!(!inactive.passed);
        assert!(active.passed);
    }

    #[test]
    fn observed_order_recovers_second_order() {
        let order = observed_order(4.0e-4, 1.0e-4, 2.0).unwrap();
        assert!((order - 2.0).abs() < 1.0e-14);
    }

    #[test]
    fn monotonicity_is_explicit_and_validated() {
        assert!(strictly_decreasing(&[1.0, 0.5, 0.2, 0.1]).unwrap());
        assert!(!strictly_decreasing(&[1.0, 0.5, 0.6, 0.1]).unwrap());
        assert!(matches!(
            strictly_decreasing(&[1.0]),
            Err(OracleError::InsufficientSamples(1))
        ));
    }
}

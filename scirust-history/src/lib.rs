//! Domain-neutral retained-history primitives for SciRust.
//!
//! This crate is a de-specialized harvest of history-management mechanisms
//! validated in `Memorithm/nonlocal-relativity-v2`. The source concepts are
//! principally the complete/bounded history work from research PR #1, the
//! true-position/non-uniform composition from PRs #6 and #7, the bounded
//! qualification experiment from PR #10, and the consolidated non-uniform
//! memory work from PR #13.
//!
//! The implementation here is adapted/reimplemented around generic `Value`
//! and `Position` types; it is not a verbatim move of the relativistic code.
//! General-relativistic coordinates, connections, transport, proper time and
//! memory laws deliberately remain outside this crate.
//!
//! # Semantic rules
//!
//! - Every retained value carries its original logical position explicitly.
//! - Accepted positions are strictly increasing. Duplicate, reversed, or
//!   incomparable positions fail closed and do not mutate history.
//! - [`CompleteHistory`] is the retained-history **reference** policy. The word
//!   reference does not imply an exact mathematical oracle for a downstream
//!   numerical method.
//! - [`BoundedHistory`] is always classified as an explicit approximation,
//!   even before its capacity has forced an eviction.
//! - Bounded eviction is observable through [`PushOutcome`]; truncation never
//!   occurs silently.
//! - No uniform-spacing assumption is made or reconstructed from an index.

#![no_std]
#![forbid(unsafe_code)]
#![deny(missing_docs)]

extern crate alloc;

mod kernel;
mod transform;
pub use kernel::HistoryKernel;
pub use transform::{
    FiniteHistoryWeight, HistoryTransform, HistoryWeight, HistoryWeightError,
    IdentityHistoryTransform, IdentityHistoryWeight, validate_history_weight,
};

use alloc::collections::VecDeque;
use alloc::vec::Vec;
use core::cmp::Ordering;
use core::fmt;
use core::num::NonZeroUsize;

/// One retained historical value at its original logical position.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct HistoryEntry<Value, Position> {
    value: Value,
    position: Position,
}

impl<Value, Position> HistoryEntry<Value, Position> {
    /// Construct a history entry without altering or reconstructing its position.
    #[must_use]
    pub const fn new(value: Value, position: Position) -> Self {
        Self { value, position }
    }

    /// Return the retained value.
    #[must_use]
    pub const fn value(&self) -> &Value {
        &self.value
    }

    /// Return the original logical position supplied at insertion.
    #[must_use]
    pub const fn position(&self) -> &Position {
        &self.position
    }

    /// Decompose this entry into its retained value and original position.
    #[must_use]
    pub fn into_parts(self) -> (Value, Position) {
        (self.value, self.position)
    }
}

/// Scientific status of a retention policy relative to complete retained history.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum HistoryFidelity {
    /// Complete retained history used as the reference behavior.
    Reference,
    /// A deliberately reduced retained history that must be evaluated as an approximation.
    Approximation,
}

/// Storage-only retention policy for a history backend.
///
/// Product policy such as tenant quotas, attention semantics, authorization or
/// retrieval ranking does not belong in this type.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum RetentionPolicy {
    /// Retain every accepted entry.
    Complete,
    /// Retain at most the given non-zero number of most recent entries.
    Bounded(NonZeroUsize),
    /// Retain a fixed-capacity age-adaptive approximation: the first sample
    /// remains an anchor, the newest tail remains exact, and older samples are
    /// thinned deterministically with a dyadic insertion-age priority.
    AdaptiveAge {
        /// Maximum number of retained entries.
        capacity: NonZeroUsize,
        /// Number of newest entries protected from age-based eviction.
        recent_samples: NonZeroUsize,
    },
}

impl RetentionPolicy {
    /// Construct a bounded policy, rejecting zero capacity explicitly.
    pub fn bounded(capacity: usize) -> Result<Self, HistoryError> {
        NonZeroUsize::new(capacity)
            .map(Self::Bounded)
            .ok_or(HistoryError::ZeroCapacity)
    }

    /// Construct an age-adaptive bounded policy.
    ///
    /// The capacity must exceed the non-zero exact recent-tail length so an
    /// older anchor/summary slot always exists.
    pub fn adaptive_age(capacity: usize, recent_samples: usize) -> Result<Self, HistoryError> {
        let Some(capacity_nonzero) = NonZeroUsize::new(capacity)
        else
        {
            return Err(HistoryError::InvalidAdaptiveAgeConfiguration {
                capacity,
                recent_samples,
            });
        };
        let Some(recent_nonzero) = NonZeroUsize::new(recent_samples)
        else
        {
            return Err(HistoryError::InvalidAdaptiveAgeConfiguration {
                capacity,
                recent_samples,
            });
        };
        if capacity <= recent_samples
        {
            return Err(HistoryError::InvalidAdaptiveAgeConfiguration {
                capacity,
                recent_samples,
            });
        }

        Ok(Self::AdaptiveAge {
            capacity: capacity_nonzero,
            recent_samples: recent_nonzero,
        })
    }

    /// Return the scientific fidelity classification implied by this policy.
    #[must_use]
    pub const fn fidelity(self) -> HistoryFidelity {
        match self
        {
            Self::Complete => HistoryFidelity::Reference,
            Self::Bounded(_) | Self::AdaptiveAge { .. } => HistoryFidelity::Approximation,
        }
    }

    /// Return the fixed maximum retained sample count, or `None` for complete retention.
    #[must_use]
    pub const fn capacity(self) -> Option<usize> {
        match self
        {
            Self::Complete => None,
            Self::Bounded(capacity) | Self::AdaptiveAge { capacity, .. } => Some(capacity.get()),
        }
    }

    /// Return the protected exact recent-tail length for an age-adaptive
    /// policy, or `None` for other policies.
    #[must_use]
    pub const fn recent_samples(self) -> Option<usize> {
        match self
        {
            Self::AdaptiveAge { recent_samples, .. } => Some(recent_samples.get()),
            Self::Complete | Self::Bounded(_) => None,
        }
    }
}

/// Failure to accept a history entry while preserving backend invariants.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HistoryError {
    /// A bounded backend was requested with zero retained capacity.
    ZeroCapacity,
    /// An age-adaptive backend has no room for both its exact recent tail and
    /// at least one older retained sample.
    InvalidAdaptiveAgeConfiguration {
        /// Requested total retained capacity.
        capacity: usize,
        /// Requested exact recent-tail length.
        recent_samples: usize,
    },
    /// The proposed position is equal to the most recently accepted position.
    DuplicatePosition {
        /// Number of samples successfully accepted before the rejected proposal.
        observed_samples: usize,
    },
    /// The proposed position precedes the most recently accepted position.
    OutOfOrderPosition {
        /// Number of samples successfully accepted before the rejected proposal.
        observed_samples: usize,
    },
    /// The proposed position cannot be ordered with itself or the latest retained position.
    IncomparablePosition {
        /// Number of samples successfully accepted before the rejected proposal.
        observed_samples: usize,
    },
}

impl fmt::Display for HistoryError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self
        {
            Self::ZeroCapacity =>
            {
                formatter.write_str("history retention capacity must be non-zero")
            },
            Self::InvalidAdaptiveAgeConfiguration {
                capacity,
                recent_samples,
            } => write!(
                formatter,
                "adaptive-age history requires capacity > recent_samples > 0; got capacity={capacity}, recent_samples={recent_samples}"
            ),
            Self::DuplicatePosition { observed_samples } => write!(
                formatter,
                "history position duplicates the latest accepted position after {observed_samples} samples"
            ),
            Self::OutOfOrderPosition { observed_samples } => write!(
                formatter,
                "history position is older than the latest accepted position after {observed_samples} samples"
            ),
            Self::IncomparablePosition { observed_samples } => write!(
                formatter,
                "history position is not comparable after {observed_samples} samples"
            ),
        }
    }
}

/// Two physical slices forming one oldest-to-newest logical history sequence.
pub type HistorySlices<'a, Value, Position> = (
    &'a [HistoryEntry<Value, Position>],
    &'a [HistoryEntry<Value, Position>],
);

/// Immutable oldest-to-newest view of entries retained by one backend.
///
/// A circular backend can expose two physical slices while preserving one
/// logical sequence. Consumers should normally use [`Self::iter`] or
/// [`Self::get`] rather than assuming contiguous storage.
#[derive(Debug, Clone, Copy)]
pub struct HistoryView<'a, Value, Position> {
    first: &'a [HistoryEntry<Value, Position>],
    second: &'a [HistoryEntry<Value, Position>],
    policy: RetentionPolicy,
    observed_samples: usize,
}

impl<'a, Value, Position> HistoryView<'a, Value, Position> {
    /// Construct a view from up to two oldest-to-newest physical slices.
    ///
    /// This is public so downstream crates can implement [`HistoryBackend`]
    /// without delegating storage to SciRust's built-in backends. The caller is
    /// responsible for supplying slices in logical order and accounting that is
    /// consistent with the backend's retention policy.
    #[must_use]
    pub const fn from_slices(
        first: &'a [HistoryEntry<Value, Position>],
        second: &'a [HistoryEntry<Value, Position>],
        policy: RetentionPolicy,
        observed_samples: usize,
    ) -> Self {
        Self {
            first,
            second,
            policy,
            observed_samples,
        }
    }

    /// Return the physical slices that form this oldest-to-newest logical view.
    #[must_use]
    pub const fn as_slices(&self) -> HistorySlices<'a, Value, Position> {
        (self.first, self.second)
    }

    /// Iterate retained entries from oldest logical position to newest.
    pub fn iter(&self) -> impl Iterator<Item = &HistoryEntry<Value, Position>> {
        self.first.iter().chain(self.second.iter())
    }

    /// Return one retained entry by oldest-to-newest logical index.
    #[must_use]
    pub fn get(&self, index: usize) -> Option<&HistoryEntry<Value, Position>> {
        if index < self.first.len()
        {
            self.first.get(index)
        }
        else
        {
            self.second.get(index - self.first.len())
        }
    }

    /// Return the number of entries currently retained.
    #[must_use]
    pub const fn retained_samples(&self) -> usize {
        self.first.len() + self.second.len()
    }

    /// Return the number of entries successfully accepted during this backend lifetime.
    #[must_use]
    pub const fn observed_samples(&self) -> usize {
        self.observed_samples
    }

    /// Return the configured storage-only retention policy.
    #[must_use]
    pub const fn retention_policy(&self) -> RetentionPolicy {
        self.policy
    }

    /// Return the reference-versus-approximation classification of the policy.
    #[must_use]
    pub const fn fidelity(&self) -> HistoryFidelity {
        self.policy.fidelity()
    }

    /// Return whether every successfully observed sample remains retained right now.
    ///
    /// This is distinct from [`Self::fidelity`]. A bounded policy remains an
    /// approximation policy even before its capacity has caused an eviction.
    #[must_use]
    pub const fn retains_all_observed(&self) -> bool {
        self.retained_samples() == self.observed_samples
    }

    /// Return whether no entries are currently retained.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.first.is_empty() && self.second.is_empty()
    }
}

/// Result of one successful history insertion.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PushOutcome<Value, Position> {
    evicted: Option<HistoryEntry<Value, Position>>,
    retained_samples: usize,
    observed_samples: usize,
}

impl<Value, Position> PushOutcome<Value, Position> {
    /// Construct a successful insertion outcome.
    ///
    /// This is public so downstream implementations of [`HistoryBackend`] can
    /// report observable eviction and accounting without delegating to a built-in
    /// backend. The implementor is responsible for supplying consistent counts.
    #[must_use]
    pub const fn new(
        evicted: Option<HistoryEntry<Value, Position>>,
        retained_samples: usize,
        observed_samples: usize,
    ) -> Self {
        Self {
            evicted,
            retained_samples,
            observed_samples,
        }
    }

    /// Return the entry evicted by this insertion, if bounded retention removed one.
    #[must_use]
    pub const fn evicted(&self) -> Option<&HistoryEntry<Value, Position>> {
        self.evicted.as_ref()
    }

    /// Consume the outcome and return the evicted entry, if any.
    #[must_use]
    pub fn into_evicted(self) -> Option<HistoryEntry<Value, Position>> {
        self.evicted
    }

    /// Return the retained sample count after the insertion.
    #[must_use]
    pub const fn retained_samples(&self) -> usize {
        self.retained_samples
    }

    /// Return the total successfully observed sample count after the insertion.
    #[must_use]
    pub const fn observed_samples(&self) -> usize {
        self.observed_samples
    }
}

/// Minimal contract implemented by deterministic history storage backends.
pub trait HistoryBackend<Value, Position>
where
    Position: PartialOrd,
{
    /// Return the storage-only retention policy.
    fn retention_policy(&self) -> RetentionPolicy;

    /// Return the total number of successfully accepted samples.
    fn observed_samples(&self) -> usize;

    /// Return an immutable ordered view of retained history.
    fn view(&self) -> HistoryView<'_, Value, Position>;

    /// Accept one entry while enforcing strict position ordering.
    ///
    /// Failed insertion must leave the backend unchanged. Any bounded eviction
    /// is returned explicitly in [`PushOutcome`].
    fn push(
        &mut self,
        entry: HistoryEntry<Value, Position>,
    ) -> Result<PushOutcome<Value, Position>, HistoryError>;
}

/// Complete retained-history reference backend.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CompleteHistory<Value, Position> {
    entries: Vec<HistoryEntry<Value, Position>>,
}

impl<Value, Position> CompleteHistory<Value, Position> {
    /// Construct an empty complete-history backend.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            entries: Vec::new(),
        }
    }

    /// Construct an empty complete-history backend with reserved capacity.
    #[must_use]
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            entries: Vec::with_capacity(capacity),
        }
    }
}

impl<Value, Position> Default for CompleteHistory<Value, Position> {
    fn default() -> Self {
        Self::new()
    }
}

impl<Value, Position> HistoryBackend<Value, Position> for CompleteHistory<Value, Position>
where
    Position: PartialOrd,
{
    fn retention_policy(&self) -> RetentionPolicy {
        RetentionPolicy::Complete
    }

    fn observed_samples(&self) -> usize {
        self.entries.len()
    }

    fn view(&self) -> HistoryView<'_, Value, Position> {
        HistoryView::from_slices(
            &self.entries,
            &[],
            RetentionPolicy::Complete,
            self.entries.len(),
        )
    }

    fn push(
        &mut self,
        entry: HistoryEntry<Value, Position>,
    ) -> Result<PushOutcome<Value, Position>, HistoryError> {
        let observed_samples = self.entries.len();
        validate_next_position(self.entries.last(), entry.position(), observed_samples)?;
        self.entries.push(entry);
        Ok(PushOutcome::new(
            None,
            self.entries.len(),
            self.entries.len(),
        ))
    }
}

/// Explicit bounded retained-history approximation backend.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BoundedHistory<Value, Position> {
    capacity: NonZeroUsize,
    entries: VecDeque<HistoryEntry<Value, Position>>,
    observed_samples: usize,
}

impl<Value, Position> BoundedHistory<Value, Position> {
    /// Construct an empty bounded-history backend.
    ///
    /// Capacity is a generic storage limit and may be one. Numerical kernels
    /// needing a larger stencil must enforce that requirement in their own layer.
    pub fn new(capacity: usize) -> Result<Self, HistoryError> {
        let capacity = NonZeroUsize::new(capacity).ok_or(HistoryError::ZeroCapacity)?;
        Ok(Self {
            capacity,
            entries: VecDeque::with_capacity(capacity.get()),
            observed_samples: 0,
        })
    }

    /// Return the maximum number of entries retained by this backend.
    #[must_use]
    pub const fn capacity(&self) -> usize {
        self.capacity.get()
    }
}

impl<Value, Position> HistoryBackend<Value, Position> for BoundedHistory<Value, Position>
where
    Position: PartialOrd,
{
    fn retention_policy(&self) -> RetentionPolicy {
        RetentionPolicy::Bounded(self.capacity)
    }

    fn observed_samples(&self) -> usize {
        self.observed_samples
    }

    fn view(&self) -> HistoryView<'_, Value, Position> {
        let (first, second) = self.entries.as_slices();
        HistoryView::from_slices(
            first,
            second,
            RetentionPolicy::Bounded(self.capacity),
            self.observed_samples,
        )
    }

    fn push(
        &mut self,
        entry: HistoryEntry<Value, Position>,
    ) -> Result<PushOutcome<Value, Position>, HistoryError> {
        validate_next_position(self.entries.back(), entry.position(), self.observed_samples)?;

        let evicted = if self.entries.len() == self.capacity.get()
        {
            self.entries.pop_front()
        }
        else
        {
            None
        };
        self.entries.push_back(entry);
        self.observed_samples += 1;

        Ok(PushOutcome::new(
            evicted,
            self.entries.len(),
            self.observed_samples,
        ))
    }
}

/// Deterministic fixed-capacity age-adaptive retained-history approximation.
///
/// This backend preserves the first accepted sample as an anchor and protects
/// the configured newest recent_samples entries exactly. When full, one older
/// interior entry is evicted before each insertion. Candidates are ranked by
/// the number of trailing zero bits in their insertion ordinal: powers of two
/// and their multiples survive longer, producing deterministic
/// dyadic/multiresolution thinning of older history without reconstructing or
/// altering logical positions.
///
/// The policy is always an approximation, even before the first eviction. It
/// changes retention only; downstream kernels remain responsible for
/// interpreting the resulting non-uniform positions and for quantifying
/// approximation error against CompleteHistory.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AdaptiveAgeHistory<Value, Position> {
    capacity: NonZeroUsize,
    recent_samples: NonZeroUsize,
    entries: Vec<HistoryEntry<Value, Position>>,
    ordinals: Vec<usize>,
    observed_samples: usize,
    total_evictions: usize,
    last_selection_checks: usize,
    total_selection_checks: usize,
}

/// Deterministic accounting for AdaptiveAgeHistory.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AdaptiveAgeDiagnostics {
    /// Entries currently retained.
    pub retained_samples: usize,
    /// Entries accepted during this backend lifetime.
    pub observed_samples: usize,
    /// Entries removed by age-adaptive thinning.
    pub total_evictions: usize,
    /// Eviction candidates examined by the most recent insertion.
    pub last_selection_checks: usize,
    /// Cumulative eviction candidates examined over the backend lifetime.
    pub total_selection_checks: usize,
}

impl<Value, Position> AdaptiveAgeHistory<Value, Position> {
    /// Construct an empty age-adaptive history.
    ///
    /// capacity must exceed the non-zero recent_samples tail length.
    pub fn new(capacity: usize, recent_samples: usize) -> Result<Self, HistoryError> {
        let policy = RetentionPolicy::adaptive_age(capacity, recent_samples)?;
        let RetentionPolicy::AdaptiveAge {
            capacity,
            recent_samples,
        } = policy
        else
        {
            unreachable!()
        };

        Ok(Self {
            capacity,
            recent_samples,
            entries: Vec::with_capacity(capacity.get()),
            ordinals: Vec::with_capacity(capacity.get()),
            observed_samples: 0,
            total_evictions: 0,
            last_selection_checks: 0,
            total_selection_checks: 0,
        })
    }

    /// Return the fixed retained-entry capacity.
    #[must_use]
    pub const fn capacity(&self) -> usize {
        self.capacity.get()
    }

    /// Return the protected exact recent-tail length.
    #[must_use]
    pub const fn recent_samples(&self) -> usize {
        self.recent_samples.get()
    }

    /// Return deterministic retention-operation accounting.
    #[must_use]
    pub const fn diagnostics(&self) -> AdaptiveAgeDiagnostics {
        AdaptiveAgeDiagnostics {
            retained_samples: self.entries.len(),
            observed_samples: self.observed_samples,
            total_evictions: self.total_evictions,
            last_selection_checks: self.last_selection_checks,
            total_selection_checks: self.total_selection_checks,
        }
    }

    fn eviction_index(&self) -> (usize, usize) {
        debug_assert_eq!(self.entries.len(), self.capacity.get());
        let protected_existing_recent = self.recent_samples.get().saturating_sub(1);
        let candidate_end = self.entries.len() - protected_existing_recent;
        debug_assert!(candidate_end > 1);

        let mut best_index = 1;
        let mut best_priority = self.ordinals[best_index].trailing_zeros();
        let mut checks = 1;
        for index in 2..candidate_end
        {
            checks += 1;
            let priority = self.ordinals[index].trailing_zeros();
            if priority < best_priority
            {
                best_index = index;
                best_priority = priority;
            }
        }
        (best_index, checks)
    }
}

impl<Value, Position> HistoryBackend<Value, Position> for AdaptiveAgeHistory<Value, Position>
where
    Position: PartialOrd,
{
    fn retention_policy(&self) -> RetentionPolicy {
        RetentionPolicy::AdaptiveAge {
            capacity: self.capacity,
            recent_samples: self.recent_samples,
        }
    }

    fn observed_samples(&self) -> usize {
        self.observed_samples
    }

    fn view(&self) -> HistoryView<'_, Value, Position> {
        HistoryView::from_slices(
            &self.entries,
            &[],
            self.retention_policy(),
            self.observed_samples,
        )
    }

    fn push(
        &mut self,
        entry: HistoryEntry<Value, Position>,
    ) -> Result<PushOutcome<Value, Position>, HistoryError> {
        validate_next_position(self.entries.last(), entry.position(), self.observed_samples)?;

        self.last_selection_checks = 0;
        let evicted = if self.entries.len() == self.capacity.get()
        {
            let (index, checks) = self.eviction_index();
            self.last_selection_checks = checks;
            self.total_selection_checks += checks;
            self.ordinals.remove(index);
            self.total_evictions += 1;
            Some(self.entries.remove(index))
        }
        else
        {
            None
        };

        self.entries.push(entry);
        self.ordinals.push(self.observed_samples);
        self.observed_samples += 1;

        Ok(PushOutcome::new(
            evicted,
            self.entries.len(),
            self.observed_samples,
        ))
    }
}

fn validate_next_position<Value, Position>(
    latest: Option<&HistoryEntry<Value, Position>>,
    position: &Position,
    observed_samples: usize,
) -> Result<(), HistoryError>
where
    Position: PartialOrd,
{
    if position.partial_cmp(position).is_none()
    {
        return Err(HistoryError::IncomparablePosition { observed_samples });
    }

    let Some(latest) = latest
    else
    {
        return Ok(());
    };

    match latest.position().partial_cmp(position)
    {
        Some(Ordering::Less) => Ok(()),
        Some(Ordering::Equal) => Err(HistoryError::DuplicatePosition { observed_samples }),
        Some(Ordering::Greater) => Err(HistoryError::OutOfOrderPosition { observed_samples }),
        None => Err(HistoryError::IncomparablePosition { observed_samples }),
    }
}

#[cfg(test)]
mod tests {
    use alloc::vec::Vec;

    use super::{
        AdaptiveAgeHistory, BoundedHistory, CompleteHistory, HistoryBackend, HistoryEntry,
        HistoryError, HistoryFidelity, HistoryView, PushOutcome, RetentionPolicy,
    };

    #[test]
    fn complete_history_is_deterministic_reference() {
        let mut history = CompleteHistory::new();
        history.push(HistoryEntry::new(10_u64, 1_u64)).unwrap();
        history.push(HistoryEntry::new(20_u64, 3_u64)).unwrap();
        history.push(HistoryEntry::new(30_u64, 9_u64)).unwrap();

        let view = history.view();
        assert_eq!(view.fidelity(), HistoryFidelity::Reference);
        assert_eq!(view.retention_policy(), RetentionPolicy::Complete);
        assert_eq!(view.retained_samples(), 3);
        assert_eq!(view.observed_samples(), 3);
        assert!(view.retains_all_observed());
        assert_eq!(*view.get(0).unwrap().position(), 1);
        assert_eq!(*view.get(1).unwrap().position(), 3);
        assert_eq!(*view.get(2).unwrap().position(), 9);
    }

    #[test]
    fn bounded_policy_is_explicit_approximation_before_eviction() {
        let mut history = BoundedHistory::new(4).unwrap();
        history.push(HistoryEntry::new(7_u8, 10_u64)).unwrap();

        let view = history.view();
        assert_eq!(view.fidelity(), HistoryFidelity::Approximation);
        assert_eq!(view.retention_policy().capacity(), Some(4));
        assert!(view.retains_all_observed());
    }

    #[test]
    fn bounded_eviction_is_observable_and_preserves_original_positions() {
        let mut history = BoundedHistory::new(2).unwrap();
        history.push(HistoryEntry::new('a', 1_u64)).unwrap();
        history.push(HistoryEntry::new('b', 8_u64)).unwrap();
        let outcome = history.push(HistoryEntry::new('c', 32_u64)).unwrap();

        let evicted = outcome.evicted().expect("third insertion must evict");
        assert_eq!(*evicted.value(), 'a');
        assert_eq!(*evicted.position(), 1);
        assert_eq!(outcome.retained_samples(), 2);
        assert_eq!(outcome.observed_samples(), 3);

        let positions: Vec<u64> = history
            .view()
            .iter()
            .map(|entry| *entry.position())
            .collect();
        assert_eq!(positions, [8, 32]);
        assert!(!history.view().retains_all_observed());
    }

    #[test]
    fn bounded_wraparound_keeps_oldest_to_newest_iteration() {
        let mut history = BoundedHistory::new(3).unwrap();
        for position in 1_u64..=8
        {
            history
                .push(HistoryEntry::new(position * 10, position))
                .unwrap();
        }

        let positions: Vec<u64> = history
            .view()
            .iter()
            .map(|entry| *entry.position())
            .collect();
        assert_eq!(positions, [6, 7, 8]);
        assert_eq!(history.view().retained_samples(), 3);
        assert_eq!(history.view().observed_samples(), 8);
    }

    #[test]
    fn bounded_covering_all_samples_matches_complete_history_exactly() {
        let samples = [(2.5_f64, 0.0_f64), (-7.0, 0.125), (11.25, 0.9)];
        let mut complete = CompleteHistory::new();
        let mut bounded = BoundedHistory::new(samples.len()).unwrap();

        for (value, position) in samples
        {
            complete.push(HistoryEntry::new(value, position)).unwrap();
            let outcome = bounded.push(HistoryEntry::new(value, position)).unwrap();
            assert!(outcome.evicted().is_none());
        }

        assert!(bounded.view().iter().eq(complete.view().iter()));
        assert!(bounded.view().retains_all_observed());
        assert_eq!(bounded.view().fidelity(), HistoryFidelity::Approximation);
    }

    #[test]
    fn public_constructors_support_external_backend_implementors() {
        let entries = [HistoryEntry::new(7_u8, 11_u64)];
        let view = HistoryView::from_slices(&entries, &[], RetentionPolicy::Complete, 1);
        let outcome = PushOutcome::new(None::<HistoryEntry<u8, u64>>, 1, 1);

        assert_eq!(view.retained_samples(), 1);
        assert_eq!(*view.get(0).unwrap().value(), 7);
        assert_eq!(outcome.retained_samples(), 1);
        assert_eq!(outcome.observed_samples(), 1);
    }

    #[test]
    fn duplicate_position_is_rejected_without_mutation() {
        let mut history = CompleteHistory::new();
        history.push(HistoryEntry::new(1_u8, 4_u64)).unwrap();
        let before = history.clone();

        assert_eq!(
            history.push(HistoryEntry::new(2_u8, 4_u64)),
            Err(HistoryError::DuplicatePosition {
                observed_samples: 1
            })
        );
        assert_eq!(history, before);
    }

    #[test]
    fn reversed_position_is_rejected_without_mutation() {
        let mut history = BoundedHistory::new(2).unwrap();
        history.push(HistoryEntry::new(1_u8, 9_u64)).unwrap();
        let before = history.clone();

        assert_eq!(
            history.push(HistoryEntry::new(2_u8, 3_u64)),
            Err(HistoryError::OutOfOrderPosition {
                observed_samples: 1
            })
        );
        assert_eq!(history, before);
    }

    #[test]
    fn incomparable_position_is_rejected_without_mutation() {
        let mut history = CompleteHistory::<u8, f64>::new();
        assert_eq!(
            history.push(HistoryEntry::new(1, f64::NAN)),
            Err(HistoryError::IncomparablePosition {
                observed_samples: 0
            })
        );
        assert!(history.view().is_empty());
    }

    #[test]
    fn bounded_zero_capacity_is_rejected() {
        assert_eq!(
            BoundedHistory::<u8, u64>::new(0),
            Err(HistoryError::ZeroCapacity)
        );
        assert_eq!(RetentionPolicy::bounded(0), Err(HistoryError::ZeroCapacity));
    }

    #[test]
    fn adaptive_age_policy_rejects_invalid_tail_geometry() {
        assert_eq!(
            AdaptiveAgeHistory::<u8, u64>::new(4, 0),
            Err(HistoryError::InvalidAdaptiveAgeConfiguration {
                capacity: 4,
                recent_samples: 0
            })
        );
        assert_eq!(
            AdaptiveAgeHistory::<u8, u64>::new(4, 4),
            Err(HistoryError::InvalidAdaptiveAgeConfiguration {
                capacity: 4,
                recent_samples: 4
            })
        );
        assert_eq!(
            RetentionPolicy::adaptive_age(0, 1),
            Err(HistoryError::InvalidAdaptiveAgeConfiguration {
                capacity: 0,
                recent_samples: 1
            })
        );
    }

    #[test]
    fn adaptive_age_retention_is_deterministic_and_preserves_recent_tail() {
        let mut history = AdaptiveAgeHistory::new(6, 2).unwrap();
        let mut evicted_positions = Vec::new();

        for position in 1_u64..=12
        {
            let outcome = history
                .push(HistoryEntry::new(position * 10, position))
                .unwrap();
            if let Some(evicted) = outcome.into_evicted()
            {
                evicted_positions.push(*evicted.position());
            }
        }

        let positions: Vec<u64> = history
            .view()
            .iter()
            .map(|entry| *entry.position())
            .collect();
        assert_eq!(positions, [1, 5, 7, 9, 11, 12]);
        assert_eq!(evicted_positions, [2, 4, 6, 8, 3, 10]);
        assert_eq!(history.view().fidelity(), HistoryFidelity::Approximation);
        assert_eq!(history.view().retention_policy().capacity(), Some(6));
        assert_eq!(history.view().retention_policy().recent_samples(), Some(2));

        let diagnostics = history.diagnostics();
        assert_eq!(diagnostics.retained_samples, 6);
        assert_eq!(diagnostics.observed_samples, 12);
        assert_eq!(diagnostics.total_evictions, 6);
        assert_eq!(diagnostics.last_selection_checks, 4);
        assert_eq!(diagnostics.total_selection_checks, 24);
    }

    #[test]
    fn adaptive_age_failed_insert_is_transactional() {
        let mut history = AdaptiveAgeHistory::new(5, 2).unwrap();
        for position in 1_u64..=7
        {
            history.push(HistoryEntry::new(position, position)).unwrap();
        }
        let before = history.clone();

        assert_eq!(
            history.push(HistoryEntry::new(99, 7)),
            Err(HistoryError::DuplicatePosition {
                observed_samples: 7
            })
        );
        assert_eq!(history, before);
    }

    #[test]
    fn bit_replay_is_deterministic_for_identical_insertions() {
        let samples = [
            (1.0_f64, 0.0_f64),
            (-0.0, 0.25),
            (f64::from_bits(0x3fe0_0000_0000_0001), 1.75),
        ];
        let mut left = CompleteHistory::new();
        let mut right = CompleteHistory::new();

        for (value, position) in samples
        {
            left.push(HistoryEntry::new(value, position)).unwrap();
            right.push(HistoryEntry::new(value, position)).unwrap();
        }

        for (left_entry, right_entry) in left.view().iter().zip(right.view().iter())
        {
            assert_eq!(left_entry.value().to_bits(), right_entry.value().to_bits());
            assert_eq!(
                left_entry.position().to_bits(),
                right_entry.position().to_bits()
            );
        }
    }
}

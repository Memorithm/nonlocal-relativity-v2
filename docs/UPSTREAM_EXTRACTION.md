# Explicit upstream route for this frozen snapshot

Audit finding: B-NR-02. Documentation-only reconciliation; no generic crate, numerical method or dependency is changed.

## Observed source

Source repository: Memorithm/nonlocal-relativity-v2  
Source base: `2c11a00b88357aae7950f40c6e6ae466713ec0f5`

The root Cargo metadata still identifies a historical SciRust workspace. `FORK.md` classifies it as a frozen fork snapshot. The previous root overview is retained unchanged at `docs/SCIRUST-SNAPSHOT-README.md` for attribution and historical provenance.

## Ownership and candidate classification

| Candidate surface | Canonical owner | Treatment |
| --- | --- | --- |
| `scirust-history/` | SciRust | Compare exact implementation/test lineage first; do not promote a stale copy |
| `scirust-nonlocal-relativity/src/adaptive_control.rs` | SciRust for reusable control primitives; specialised layer for relativistic coupling | Classify generic versus domain assumptions before any split |
| `scirust-fractional/` | SciRust | Preserve independent reference oracles and declared approximation contracts |
| `scirust-nonlocal-relativity/` domain-specific code | Specialised relativistic layer | Retain model and coordinate assumptions with each candidate |
| `experiments/nonlocal-relativity-v2/` | Research protocol/evidence producer | Preserve fixtures, negative results and scoped interpretation |
| Unrelated embedded SciRust workspace | Canonical SciRust | Frozen support only; never bulk-sync |

This table assigns responsibilities; it does not assert that extraction, scientific promotion, API convergence or canonical destination tests have already occurred.

## Gate sequence for a future extraction

1. Resolve the actual default branches of both repositories and record their full HEAD and tree SHAs.
2. Read both repositories' AGENTS, mandatory roadmaps, licences and notices; confirm the candidate is permitted despite the fork freeze.
3. For every candidate file, compare source bytes with the exact canonical target version. Classify it as identical, older, divergent, research-only or genuinely missing. Do not copy if canonical code is equivalent/newer.
4. Publish a versioned manifest with source SHA/tree, target base SHA/tree, file paths/digests, licence/attribution, scientific category, dependencies, oracle scope and test commands. A mutable branch name or repository name is insufficient provenance.
5. Propose the smallest coherent slice **to canonical SciRust**. Generic history/adaptive-control/fractional code must not gain an implicit dependency on the specialised experiment.
6. Execute source and destination tests at exact SHAs, including reference/parity, non-finite refusal, determinism and relevant approximation/error contracts. Retain input/output hashes and CI URLs. A standalone source check is not destination qualification.
7. Recheck diff, target base, exact candidate head, review state and all required CI before a policy-permitted merge. Revalidate the integrated destination and record its merge SHA.
8. Only then consider replacing specialised adapters with the qualified canonical dependency. Never delete source research evidence or archive the fork as an implicit side effect.

## Required evidence record

A promotion record must include: finding/candidate ID, exact source and target identities, affected files and hashes, destination API/toolchain/features, scientific limitations, commands and observed results, required CI, review disposition, merge SHA and post-merge validation.

Missing evidence means **not promoted**. Failure means **rejected or blocked**, not resolved. Nonlocal numerical self-consistency never upgrades to empirical or field-equation validation.

## Status of this reconciliation

- Root landing page and historical-snapshot distinction: supplied by this documentation candidate.
- Ownership and explicit extraction procedure: supplied, without altering the frozen platform.
- Physical extraction/upstream promotion: **not performed**.
- Dependency/API convergence: **not assumed**.
- Numerical and scientific methods: **unchanged**.

# Nonlocal relativity — frozen SciRust fork snapshot

This repository is a **frozen platform snapshot with a specialised research lineage**, not the canonical SciRust distribution. Start with [FORK.md](FORK.md) before changing or consuming it.

| Surface | Responsibility |
| --- | --- |
| Canonical generic scientific platform | [Memorithm/scirust](https://github.com/Memorithm/scirust) |
| `scirust-nonlocal-relativity/` | Historical specialised test-particle memory implementations and local adapters |
| `experiments/nonlocal-relativity-v2/` | Deterministic scenarios, named numerical oracles and scoped evidence |
| This fork's other SciRust crates | Frozen historical support; no new generic platform development here |

The previous SciRust overview is preserved verbatim in [SCIRUST-SNAPSHOT-README.md](SCIRUST-SNAPSHOT-README.md). Its capability counts and measurements are historical; they do not describe the current canonical platform or independently qualify this fork.

## Scientific boundary

Established-general-relativity implementation tests and the experimental phenomenological nonlocal layer are different evidence classes. The current nonlocal model operates on a fixed background metric/connection; it does not change Einstein field equations, stress-energy, curvature sourcing or backreaction.

Determinism, finite-grid agreement and numerical self-convergence are not empirical validation, a mathematical convergence proof or evidence of new physics. Preserve each protocol's assumptions, approximation error and coordinate dependence.

Read the [specialised crate overview](scirust-nonlocal-relativity/README.md), [experiment protocol](experiments/nonlocal-relativity-v2/README.md) and [research paper](paper/nonlocal-relativity-v2.md) before interpreting outputs.

## Maintenance and upstream routing

No bulk merge, platform refresh, destructive extraction or sibling-fork creation is authorised by this overview. Documentation may restate the freeze notice. Specialised patches may be proposed to canonical SciRust only through the exact source/target and destination-validation process in [docs/UPSTREAM_EXTRACTION.md](docs/UPSTREAM_EXTRACTION.md).

Generic history, adaptive-control and fractional primitives belong to canonical SciRust. Identify equivalent/newer canonical implementations before proposing any extraction; a path dependency in this historical tree is not a published cross-repository contract.

## Reproduction and licensing

Existing Cargo locks, CI gates and experiment commands remain unchanged. This documentation does not claim a new executed campaign or successful upstream promotion. Repository-specific required CI must pass on an exact candidate before any permitted merge.

[LICENSE](LICENSE), [LICENSE.md](LICENSE.md) and [LICENSING.md](LICENSING.md) retain the existing PolyForm Noncommercial 1.0.0 terms and required attribution. No relicensing occurs.

Agents must read [AGENTS.md](AGENTS.md) and its off-default-branch roadmap before decisions. This README does not override `FORK.md`, holdouts, merge restrictions or scientific gates.

# Nonlocal Relativity Agent Bootstrap Contract

This repository has a persistent agent-oriented execution roadmap that is intentionally kept off the default branch.

## Mandatory first step

Before any autonomous coding, scientific interpretation, numerical-method change, cross-repository promotion, architectural decision, PR creation, or merge decision, read the current roadmap from:

`origin/agent/ecosystem-roadmap:.agent/NONLOCAL_RELATIVITY_ECOSYSTEM_ROADMAP.yaml`

Recommended command:

```bash
git fetch origin agent/ecosystem-roadmap && \
git show origin/agent/ecosystem-roadmap:.agent/NONLOCAL_RELATIVITY_ECOSYSTEM_ROADMAP.yaml
```

If the roadmap cannot be fetched or read, fail closed: do not make a major scientific-interpretation, architecture, cross-repository promotion, integration-contract, or merge decision. Read-only diagnosis is allowed.

## Mandatory reread points

Reread the roadmap:

1. at the start of every agent session on this repository;
2. before selecting the next major task or research phase;
3. before interpreting, promoting, or reclassifying a scientific result;
4. before any promotion or integration with another Memorithm repository;
5. after any user instruction that changes scientific scope, ecosystem role, invariants, or strategy;
6. before opening or merging a PR that changes scientific semantics or cross-repository contracts.

## Scientific boundary that must never be blurred

The repository contains both established-general-relativity implementation validation and an experimental phenomenological nonlocal test-particle layer. They are different evidence classes.

Never describe numerical self-consistency, convergence, or deterministic execution of the phenomenological layer as empirical validation, a field-equation result, or evidence of new physics. Preserve the exact scientific limitations documented by the repository.

## Ecosystem boundary

This repository is a research incubator and evidence producer for SciRust. Mature reusable work may be promoted to `Memorithm/scirust` only through explicit provenance, exact source/target SHAs, targeted diffs, destination tests, and green destination CI.

Do not silently turn this repository into a second incompatible SciRust distribution. Do not absorb responsibilities owned by SciRust Hub, SciCapsule, SciRust-Verify, Forge, ElasticXxx, NNIS, FLAT-ATTENTION, or SLHAv2.

## Mandatory roadmap maintenance

Update the off-main roadmap when:

- a research candidate is promoted, rejected, or reclassified;
- a cross-repository contract is published, changed, or rejected;
- a roadmap phase changes status;
- a new exact oracle or scientific boundary is introduced;
- the promotion target or ecosystem strategy changes.

Do not merge the roadmap itself into the default branch unless the user explicitly requests it.

## Core constraints

- correctness and declared semantics dominate optimization;
- no fabricated performance or scientific novelty;
- no empirical-validation claim without empirical evidence;
- numerical self-convergence is not physical validation;
- cross-repository contracts are never assumed before they are published and tested;
- required CI must be green on the exact PR head before merge;
- missing roadmap or missing evidence causes fail-closed behavior for major decisions.

This file is only the bootstrap pointer. The off-main roadmap is the persistent source of current strategy, ecosystem interfaces, and research state.

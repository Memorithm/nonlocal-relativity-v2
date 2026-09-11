# This repository is a fork snapshot, not a science bench

Effective 2026-09-12 (Memorithm audit v3 / plan V0-2).

## What this tree actually is

The root package is named `scirust` version `0.14.0` with `repository = "https://github.com/Memorithm/scirust"`. The workspace member list is a historical SciRust platform graph, plus `experiments/nonlocal-relativity-v2` and `scirust-nonlocal-relativity`.

Size at audit: 128 564 KB. That is a platform fork, not an experiment crate.

## Canon

| What you want | Where it lives |
|---|---|
| SciRust platform | https://github.com/Memorithm/scirust |
| Relativity / nonlocal experiment | `experiments/nonlocal-relativity-v2` **inside** Memorithm/scirust |
| This repository | Frozen snapshot. Do not evolve the platform copy. |

## Freeze rules

Allowed:
- documentation that restates this notice
- extracting *only* `experiments/nonlocal-relativity-v2` + `scirust-nonlocal-relativity` patches back to canon SciRust
- archiving the repository after the extract

Forbidden:
- new SciRust platform crates here
- treating this default branch as the living SciRust
- AUTO_MERGE of platform work into this tree
- adding a 31st sibling fork

Orchestrator: classify this repository as `FORK_SNAPSHOT`. Do not schedule feature work against it.

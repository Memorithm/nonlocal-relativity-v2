# Nonlocal Relativity repository agent instructions

Before making repository changes, fetch and read the persistent off-main agent roadmap:

```bash
git fetch origin agent/ecosystem-roadmap && \
git show origin/agent/ecosystem-roadmap:.agent/NONLOCAL_RELATIVITY_ECOSYSTEM_ROADMAP.yaml
```

Treat `AGENTS.md` at the repository root as mandatory bootstrap policy.

The roadmap must be reread at every session start, before a new major task or research phase, before interpreting or promoting a scientific result, before cross-repository promotion/integration, after strategy or scientific-scope changes, and before PR or merge decisions that affect scientific semantics or ecosystem contracts.

If the roadmap is unavailable, fail closed for major scientific-interpretation, architecture, cross-repository promotion, integration-contract, or merge decisions. Do not substitute guesses for missing roadmap state.

Preserve the repository's explicit distinction between established-GR implementation validation and the experimental phenomenological nonlocal layer. Numerical convergence or reproducibility of the latter is not empirical validation or evidence of new physics.

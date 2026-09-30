# Adaptive refinement evidence, version 2

## Rejected strict-trend hypothesis

The original five-row tolerance sweep does not demonstrate strictly decreasing endpoint error. The unmodified numerical experiment was measured on the Dell in RemoteOps run 36764843501, job 110056296427, control revision 49adef81ad225e6ab285eb7e5ccc326f0c2269cf, source candidate 19f3feede69429e4e6cca2a1837cdb0250e3465e. That candidate contains current main's guards plus the explicit f64 literal needed to compile its new generic evidence test.

Embedded Heun-Euler: accepted-step counts 17,17,19,56,175. Errors 5.47906579504352642e-9, 5.47906579504352642e-9, 4.67511164899225218e-9, 1.89076272235724917e-9, 3.52346248938996648e-10.

Step doubling: counts 17,17,17,40,124. Errors 2.92432738471547213e-6, 2.92432738471547213e-6, 2.92432738471547213e-6, 1.19665282414632863e-6, 3.38715887400592547e-7.

All repeated-error rows also had bit-identical complete final states and identical accepted-step counts. This does not prove identical trajectories or identify a causal mechanism. Both strict-trend predicates were false. These negative findings are retained, not promoted to a strict convergence result.

## Explicit regression contract

The CSV protocol, five tolerances, both fine-grid references, integrators, step bounds and error computation are unchanged. The generic strictly_decreasing function is unchanged and its actual result is still emitted. The experiment additionally reports the explicitly distinct adaptive-refinement/v2 contract: aligned finite evidence, non-increasing endpoint error, non-decreasing positive accepted-step counts, at least one overall error improvement and overall count increase. Equal-error rows require bit-identical coordinate and velocity endpoints and equal counts. A plateau is counted explicitly, never called a strict improvement. No epsilon is used to hide increases.

Tests reject a one-ULP error increase, reduced step counts, entirely flat/no-progress evidence, missing/misaligned/nonfinite data, and equal scalar errors with changed endpoint bits or counts. The original strict-convergence assertion was invalid for its unchanged input protocol; this version explicitly replaces that claim, rather than reporting it as passed.

This is a bounded fine-grid numerical self-consistency contract, not a theorem of convergence, an exact analytical oracle, coordinate-invariance evidence or empirical validation of the phenomenological model.

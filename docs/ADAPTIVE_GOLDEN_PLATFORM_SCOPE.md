# Adaptive golden regression scope

The existing Linux bit patterns remain unchanged. A separate exact fixture is selected only for aarch64-apple-darwin. Every comparison remains an equality of f64::to_bits(); no tolerance, ignored test or production numerical change is introduced. The 176 accepted samples and exact affine endpoint 0.8 remain required.

## Measurement provenance

Source: f312ce4bbaf87607d269775ca3fa84aa457d42fc. Toolchain: nightly-2026-07-02, rustc 1.98.0-nightly (4c9d2bfe4ad7a65669098754964aaebe0ec1ced2), LLVM 22.1.8, aarch64-apple-darwin on the Mac mini. RemoteOps control commit 559a026b36fe7830515bfd8df06eb7d7b7b9e969, run 36756792064, job 110028969360, 2026-09-30.

Both unoptimized and opt-level 2 builds produced identical Darwin fixtures. Two repetitions per profile had bit-identical endpoints. The original Linux fixture failed under both profiles, so the discrepancy was not introduced by enabling test optimization. The existing fine fixed-step Heun-PECE comparison and identity-policy/plain-entry-point bit equality both passed. The specific elementary operation responsible for cross-platform rounding has not been isolated.

The Darwin coordinate bits are 3fee99d33872904b, 4023fbe788e8f7f6, 3ff921fb54442d18, 3f9efcc3603d09e2. Velocity bits are 3ff3209f6d6aa13b, bf84797d4428fc0d, 3bf46a69ccd61f69, 3fa361e17609ee99. Affine bits are 3fe999999999999a and memory-norm bits are 3f34510b8cc118e0.

Maximum absolute endpoint difference from the Linux anchor is 6.661338147750939e-16. The memory-norm difference is 2.2837092460148423e-15. These measurements do not establish bitwise portability, coordinate invariance or empirical model validation. Other targets continue to use the unchanged strict legacy fixture; this is not a claim that all such environments have been qualified. Future toolchain/OS changes must be investigated rather than automatically updating fixtures.

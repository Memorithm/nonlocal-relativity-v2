# Self-hosted CI execution

Linux validation runs on repository-scoped Dell x86_64 runners; native ARM64 scientific tests and the native macOS workspace check run on the Mac mini. RemoteOps owns registration and host operations.

The Windows job checks every workspace target for x86_64-pc-windows-msvc from Linux. It is explicitly a target compilation check, not native Windows execution, packaging or build-host qualification. Native Windows evidence is not claimed.

The GnuCOBOL reference and CUDA-unavailable jobs use the pinned Ubuntu 22.04 image ubuntu@sha256:b8b6ee6aa931ecd9d0d952abc34dc0e5f7c6a30c6bb71b079fe399fde0329c02. The latter also checks that no CUDA runtime is available. All original Linux validation commands and informational/gating classifications are retained. Fork PR code is excluded from internal runners and checkout credentials are not persisted.

Tests use opt-level 2 with debug assertions and integer-overflow checks explicitly retained. No test, oracle, tolerance or production profile is removed or relaxed.

Self-hosted process execution is not a security sandbox. Scientific self-consistency remains distinct from empirical validation. Successful execution does not broaden the model's scientific claims.

## Runner-specific dependency execution

The pinned Ubuntu checks now execute through `.github/scripts/nonlocal-ubuntu-check.sh`. Docker bridge DNS/outbound connectivity was unavailable on the measured Dell host, while host HTTP/HTTPS succeeded. Each container therefore explicitly uses host networking; network isolation is not claimed. Only a read-only source archive at the exact HEAD and the reviewed script are mounted. No host credentials, Docker socket or GPU devices are mounted. Both unchanged GnuCOBOL 3.1.2 baseline checks and both CUDA-no-runtime tests passed under this execution adapter in RemoteOps run 36757227221 before it was integrated.

The security audit runs the same cargo-deny 0.19.8 version and all eight dependency graphs on the host rather than building a Docker action through the unavailable bridge. License/advisory policy is unchanged. Linux package installs wait for the dpkg lock instead of deleting it or racing other workers. Native ARM64 exact regression anchors have their observed scope documented in ADAPTIVE_GOLDEN_PLATFORM_SCOPE.md.

## TLS security remediation and sanitizer execution

RUSTSEC-2026-0285 / GHSA-2mjx-qc3c-rqvc (published 2026-09-14) affects rustls 0.23.13 through 0.23.44. The lockfile now resolves rustls 0.23.45, the minimum patched version. The resolver update is confined to rustls 0.23.45, rustls-webpki 0.103.15, aws-lc-rs 1.18.1 and aws-lc-sys 0.45.0. The upstream rustls v/0.23.45 workspace requires aws-lc-rs >=1.18 and rustls-webpki >=0.103.14; no unrelated package, manifest feature, license policy or advisory exception is changed. Sources: https://rustsec.org/advisories/RUSTSEC-2026-0285 and https://github.com/rustls/rustls/security/advisories/GHSA-2mjx-qc3c-rqvc .

All eight unchanged cargo-deny 0.19.8 dependency graphs and all three optional-network compilation checks passed before publication in https://github.com/Memorithm/RemoteOps/actions/runs/36762493842. Final-head repository CI remains required before merge. The TLS transcript remains authenticated; this advisory is not represented as arbitrary handshake alteration.

The Dell host was observed with vm.overcommit_memory=2 and unlimited process address-space limits in RemoteOps run 36759858042. Linux ASan failed during shadow-memory reservation before processing fuzz input. The shared host policy is unchanged. CI retains Linux x86_64 sanitizer-instrumented builds and runs both unchanged parser harnesses on native Mac ARM64 with AddressSanitizer, cargo-fuzz 0.13.2 and nightly-2026-07-02. Both campaigns retain -max_total_time=45 and -runs=500000. Native execution passed in RemoteOps run 36760682917, job 110042147778, control revision ad1fddb0276c1b819275d519ec569e6d0879c69c.

Mac ASan execution is not Linux runtime sanitizer evidence or exhaustive parser verification. Build and execution evidence are named separately; no sanitizer is disabled. The pre-existing informational fuzz classification is unchanged.

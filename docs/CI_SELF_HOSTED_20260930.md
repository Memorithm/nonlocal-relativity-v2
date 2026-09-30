# Self-hosted CI execution

Linux validation runs on repository-scoped Dell x86_64 runners; native ARM64 scientific tests and the native macOS workspace check run on the Mac mini. RemoteOps owns registration and host operations.

The Windows job checks every workspace target for x86_64-pc-windows-msvc from Linux. It is explicitly a target compilation check, not native Windows execution, packaging or build-host qualification. Native Windows evidence is not claimed.

The GnuCOBOL reference and CUDA-unavailable jobs use the pinned Ubuntu 22.04 image ubuntu@sha256:b8b6ee6aa931ecd9d0d952abc34dc0e5f7c6a30c6bb71b079fe399fde0329c02. The latter also checks that no CUDA runtime is available. All original Linux validation commands and informational/gating classifications are retained. Fork PR code is excluded from internal runners and checkout credentials are not persisted.

Tests use opt-level 2 with debug assertions and integer-overflow checks explicitly retained. No test, oracle, tolerance or production profile is removed or relaxed.

Self-hosted process execution is not a security sandbox. Scientific self-consistency remains distinct from empirical validation. Successful execution does not broaden the model's scientific claims.

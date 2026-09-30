#!/usr/bin/env bash
# Pinned Ubuntu userspace; host networking is explicit, not network isolation.
# Only a read-only archive and this reviewed script enter the container.
set -euo pipefail

if [[ ${1:-} == --inside ]]; then
  mode=${2:?missing check mode}
  case "$mode" in gnucobol|cuda) ;; *) exit 2 ;; esac
  export DEBIAN_FRONTEND=noninteractive
  apt-get -o APT::Update::Error-Mode=any update -qq
  apt-get install -y -qq ca-certificates git python3 sudo curl build-essential pkg-config
  mkdir -p /src
  tar -xf /source.tar -C /src
  cd /src
  if [[ "$mode" == gnucobol ]]; then
    cd scirust-finmigrate/tests/compiler-derived/gnucobol-3.1.2
    sudo apt-get update -qq
    sudo apt-get install -y -qq gnucobol
    cobc --version
    cobc --version | sed -n '1p' | grep -F '3.1.2'
    sha256sum -c SHA256SUMS
    python3 tools/run_baselines.py check
    python3 tools/run_baselines.py verify
  else
    test ! -e /dev/nvidiactl
    python3 - <<'PY'
import ctypes
try:
    ctypes.CDLL("libcuda.so.1")
except OSError:
    pass
else:
    raise SystemExit("Unexpected CUDA runtime: fallback evidence would be invalid")
PY
    export CARGO_HOME=/tmp/cargo RUSTUP_HOME=/tmp/rustup
    mkdir -p "$CARGO_HOME" "$RUSTUP_HOME"
    curl --proto '=https' --tlsv1.2 --retry 3 --location --silent --show-error --fail https://sh.rustup.rs \
      | sh -s -- --default-toolchain stable --profile minimal -y
    export PATH="$CARGO_HOME/bin:$PATH"
    rustc +stable --version --verbose
    cargo +stable test -p scirust-cuda --features cuda --locked
    cargo +stable test -p scirust-gpu --features cuda --locked
  fi
  printf 'UBUNTU_CHECK_PASSED mode=%s\n' "$mode"
  exit 0
fi

mode=${1:?usage: nonlocal-ubuntu-check.sh gnucobol|cuda}
case "$mode" in gnucobol|cuda) ;; *) echo 'Unsupported check mode' >&2; exit 2 ;; esac
image=ubuntu@sha256:b8b6ee6aa931ecd9d0d952abc34dc0e5f7c6a30c6bb71b079fe399fde0329c02
source_sha=$(git rev-parse HEAD)
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
script="$script_dir/$(basename -- "${BASH_SOURCE[0]}")"
temporary=$(mktemp -d "${RUNNER_TEMP:-${TMPDIR:-/tmp}}/nonlocal-ubuntu.XXXXXX")
trap 'rm -f -- "$temporary/source.tar"; rmdir -- "$temporary"' EXIT
git archive --format=tar "$source_sha" > "$temporary/source.tar"
printf 'UBUNTU_CHECK_SOURCE=%s MODE=%s IMAGE=%s NETWORK=host\n' "$source_sha" "$mode" "$image"
docker run --rm --network host \
  --volume "$temporary/source.tar:/source.tar:ro" \
  --volume "$script:/ci-check.sh:ro" \
  --env RUSTFLAGS='-D warnings' \
  --env CARGO_INCREMENTAL=0 \
  --env CARGO_BUILD_JOBS=4 \
  --env CARGO_PROFILE_TEST_OPT_LEVEL=2 \
  --env CARGO_PROFILE_TEST_DEBUG_ASSERTIONS=true \
  --env CARGO_PROFILE_TEST_OVERFLOW_CHECKS=true \
  "$image" bash /ci-check.sh --inside "$mode"

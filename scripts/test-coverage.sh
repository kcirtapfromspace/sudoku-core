#!/usr/bin/env bash
# Measure the entire production library, including the complete soundness suite.
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"
if (( $# )); then
  printf 'Usage: %s\n' "$0" >&2
  exit 2
fi

for program in cargo cargo-llvm-cov python3; do
  command -v "$program" >/dev/null || { printf 'Required tool missing: %s\n' "$program" >&2; exit 1; }
done
python3 -c 'import tomllib' # Python 3.11 or newer.

# Keep instrumentation separate from developer builds and baseline measurements.
export CARGO_TARGET_DIR="$repo_root/target/coverage-build"
mkdir -p target/coverage
final_status=0
check() {
  if "$@"; then
    return 0
  else
    printf 'Check failed: %s\n' "$*" >&2
    final_status=1
  fi
}

check python3 -m unittest discover -s scripts -p test_coverage.py -v
# Shared coverage counters make concurrent solver searches substantially slower.
# Run every test sequentially; soundness remains required on every acceptance run.
check cargo llvm-cov --workspace --all-targets --all-features --profile coverage --locked --no-fail-fast --no-report -- --test-threads=1
# Preserve reports on test failures so CI retains useful diagnostics.
check cargo llvm-cov report --profile coverage --lcov --output-path target/coverage/lcov.info
check cargo llvm-cov report --profile coverage --html --output-dir target/coverage/rust
check python3 scripts/coverage.py
exit "$final_status"

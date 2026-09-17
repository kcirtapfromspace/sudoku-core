#!/usr/bin/env bash
# Run Verus formal verification over the verified kernel in verification/.
#
# Requires a Verus release (https://github.com/verus-lang/verus/releases) plus
# its pinned Rust toolchain; point VERUS_BIN at the `verus` binary if it is not
# installed at ~/verus/verus. Running `verus` once prints the exact
# `rustup install <toolchain>` command if the pinned toolchain is missing.
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

VERUS_BIN="${VERUS_BIN:-$HOME/verus/verus}"
if ! command -v "$VERUS_BIN" >/dev/null; then
  printf 'Verus binary not found at %s (set VERUS_BIN or install to ~/verus)\n' "$VERUS_BIN" >&2
  exit 1
fi

status=0
for file in verification/*.rs; do
  printf '==> verifying %s\n' "$file"
  "$VERUS_BIN" --crate-type=lib "$file" || status=1
done
exit "$status"

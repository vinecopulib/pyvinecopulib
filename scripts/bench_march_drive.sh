#!/usr/bin/env bash
# Drive scripts/bench_march.py for both -march=native and -march=x86-64-v3
# builds. Uninstalls + clears the build dir between variants so neither
# benchmark reads stale .so bytes.
#
# Usage:  bash scripts/bench_march_drive.sh
# Output: /tmp/bench_march/{native,v3}.json + a printed comparison table.
set -euo pipefail

cd "$(dirname "$0")/.."

UV="${UV:-$HOME/.local/bin/uv}"
RESULTS=/tmp/bench_march
mkdir -p "$RESULTS"

build_and_bench() {
  local label=$1

  echo
  echo "================================================================"
  echo "[$label] uninstall + clean build cache"
  echo "================================================================"
  $UV pip uninstall pyvinecopulib 2>&1 | tail -2 || true
  rm -rf build

  echo
  echo "[$label] editable install"
  if [ "$label" = "v3" ]; then
    # The env var CI sets; adds -march=x86-64-v3 to the release flags.
    CIBUILDWHEEL=1 $UV pip install -e . --no-build-isolation 2>&1 | tail -3
  else
    # -march=native is opt-in; a plain build gets the redistributable baseline.
    $UV pip install -e . --no-build-isolation \
      -C cmake.define.PYVINECOPULIB_NATIVE_ARCH=ON 2>&1 | tail -3
  fi

  echo
  echo "[$label] effective release flags (from CMakeCache):"
  grep "^CMAKE_CXX_FLAGS_RELEASE" build/*/CMakeCache.txt | head -2

  echo
  echo "[$label] verifying instruction set in compiled .so:"
  if command -v objdump >/dev/null 2>&1; then
    local so
    so=$(find build -name "pyvinecopulib_ext*.so" -print -quit)
    if [ -n "$so" ]; then
      echo "    sample of AVX-512 instructions (count):"
      objdump -d --no-show-raw-insn "$so" 2>/dev/null \
        | grep -oE 'vpermq|vfmadd[0-9]{3}pd|vpaddq|vpcmp|vmovdqu64|zmm[0-9]+' \
        | sort -u | head -10
      printf "    AVX-512 zmm-register uses: "
      objdump -d "$so" 2>/dev/null | grep -cE 'zmm[0-9]+' || true
    fi
  fi

  echo
  echo "[$label] running benchmark"
  $UV run python scripts/bench_march.py \
    --label "$label" --output "$RESULTS/$label.json"
}

build_and_bench native
build_and_bench v3

echo
echo "================================================================"
echo "Side-by-side"
echo "================================================================"
$UV run python - <<'PYEOF'
import json
from pathlib import Path

native = json.loads(Path("/tmp/bench_march/native.json").read_text())
v3 = json.loads(Path("/tmp/bench_march/v3.json").read_text())

def mean_ms(xs):
  return sum(xs) / len(xs) * 1000

print(f"{'metric':<20s} {'native (ms)':>14s} {'v3 (ms)':>14s} {'v3/native':>12s}")
print("-" * 62)
for kind in ("bicop", "vine"):
  for fn in ("itau", "itau_par", "tll"):
    n_ms = mean_ms(native["results"][kind][fn])
    v_ms = mean_ms(v3["results"][kind][fn])
    ratio = v_ms / n_ms
    print(f"{kind+'/'+fn:<20s} {n_ms:>14.2f} {v_ms:>14.2f} {ratio:>11.2f}x")
PYEOF

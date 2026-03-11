#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -lt 1 ]; then
  echo "usage: $0 <foldcomp-binary>" >&2
  exit 2
fi

FOLDCOMP_BIN=$1
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
TEST_DIR=$SCRIPT_DIR

assert_close() {
  local observed=$1
  local expected=$2
  local tolerance=$3
  awk -v check="$observed" -v target="$expected" -v tol="$tolerance" \
    'BEGIN { diff = check - target; if (diff < 0) diff = -diff; if (diff > tol) { print check "!=" target; exit 1 } }'
}

assert_rmsd_exact() {
  local source=$1
  local roundtrip=$2
  local expected_bb=$3
  local expected_aa=$4
  local tolerance=${5:-0.001}
  local line bb aa
  line=$("$FOLDCOMP_BIN" rmsd "$source" "$roundtrip")
  bb=$(printf '%s\n' "$line" | cut -f5)
  aa=$(printf '%s\n' "$line" | cut -f6)
  assert_close "$bb" "$expected_bb" "$tolerance"
  assert_close "$aa" "$expected_aa" "$tolerance"
}

run_reference_cases() {
  rm -f "$TEST_DIR/test.fcz" "$TEST_DIR/test_fcz.cif" "$TEST_DIR/test.cif.fcz" "$TEST_DIR/test.cif_fcz.cif"

  "$FOLDCOMP_BIN" compress "$TEST_DIR/test.pdb"
  "$FOLDCOMP_BIN" decompress "$TEST_DIR/test.fcz" "$TEST_DIR/test_fcz.cif"
  assert_rmsd_exact "$TEST_DIR/test.pdb" "$TEST_DIR/test_fcz.cif" "0.0322624" "0.0769879"

  "$FOLDCOMP_BIN" compress "$TEST_DIR/test.cif.gz"
  "$FOLDCOMP_BIN" decompress -a "$TEST_DIR/test.cif.fcz" "$TEST_DIR/test.cif_fcz.cif"
  assert_rmsd_exact "$TEST_DIR/test.cif.gz" "$TEST_DIR/test.cif_fcz.cif" "0.0380966" "0.122364"
}

run_regression_case() {
  local name=$1
  local expected_default_bb=$2
  local expected_default_aa=$3
  local expected_max_bb=$4
  local expected_max_aa=$5
  local source="$TEST_DIR/${name}.cif.gz"
  local default_base="/tmp/foldcomp-smoke-${name}-default"
  local max_base="/tmp/foldcomp-smoke-${name}-max"

  rm -f "${default_base}" "${default_base}." "${default_base}.cif"
  "$FOLDCOMP_BIN" compress "$source" "$default_base"
  "$FOLDCOMP_BIN" decompress "${default_base}." "${default_base}.cif"
  assert_rmsd_exact "$source" "${default_base}.cif" "$expected_default_bb" "$expected_default_aa"

  rm -f "${max_base}" "${max_base}." "${max_base}.cif"
  "$FOLDCOMP_BIN" compress --max-backbone-rmsd 0.2 "$source" "$max_base"
  "$FOLDCOMP_BIN" decompress "${max_base}." "${max_base}.cif"
  assert_rmsd_exact "$source" "${max_base}.cif" "$expected_max_bb" "$expected_max_aa"
}

run_reference_cases
run_regression_case "1cfg" "0" "0" "0" "0"
run_regression_case "1adn" "0.261305" "0.311185" "0" "0"
run_regression_case "1dlp" "0.0415883" "0.336396" "0.014391" "0.132124"
run_regression_case "1hrn" "0.0730179" "0.168113" "0" "0"
run_regression_case "1jwt" "0.0783066" "0.265589" "0" "0"
run_regression_case "1lyz" "0.0870845" "0.214444" "0" "0"
run_regression_case "2ap2" "0.0359238" "0.106224" "0" "0"
run_regression_case "2gn5" "0.0237757" "0.439115" "0.00025578" "0.0256872"
run_regression_case "4kuk" "0.0411448" "0.0972185" "0" "0"
run_regression_case "7c2s" "0" "0" "0" "0"

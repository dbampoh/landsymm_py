#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
LANDSYMM_ROOT="$(cd "${PROJECT_DIR}/.." && pwd)"

HILDA_DATA="${LANDSYMM_ROOT}/data/geodata_py/HILDA+/data"
DEFAULT_OUTPUT_DIR="${HILDA_DATA}/output"
DEFAULT_GRIDLIST="${HILDA_DATA}/input_data/gridlist_in_62892_and_climate.txt"
DEFAULT_STATES=""

show_help() {
  cat <<'EOF'
Run smoothing -> (optional mini-gridlist) -> upscaling in one command.

Required:
  --states PATH                HILDA+ LULC states NetCDF (unsmoothed)

Optional inputs:
  --fm PATH                    Forest management NetCDF (optional)
  --gridlist PATH              Full gridlist (default: data/input_data/gridlist_in_62892_and_climate.txt)
  --output-dir PATH            Output directory (default: data/output)

Mode selection:
  --smooth-mode single|parallel    (default: parallel)
  --upscale-mode single|parallel   (default: parallel)

HILDA+ -> LPJ-GUESS land-cover mapping (passed through to upscaling script):
  --mapping-config PATH           Path to a custom mapping YAML
                                  (see hildaplus/config/README.md)
  --mapping-profile NAME          Named profile (e.g. lpjg_v3_default,
                                  lpjg_legacy_v1, lpjg_treecrops_as_forest)

Benchmark options:
  --benchmark                     Run benchmark chain (mini NetCDF + mini gridlist)
  --benchmark-chunks N            Smoothing benchmark chunks (default: 2)
  --benchmark-lines N             Upscaling benchmark lines (default: 200)

Parallel options:
  --smooth-workers N              Smoothing workers (default: 8)
  --smooth-chunk-size N           Smoothing chunk size (default: 100)
  --upscale-workers N             Upscaling workers (default: 8)
  --upscale-chunk-lines N         Upscaling gridlist chunk size (default: 500)
  --print-config                 Print resolved config and exit
  --dry-run                      Print commands without running
  --inspect                      Run inspection scripts after each stage

Examples:
  # Full run (parallel smoothing + parallel upscaling)
  scripts/run_chain.sh --states /path/to/states.nc

  # Benchmark run
  scripts/run_chain.sh --states /path/to/states.nc --benchmark

  # Single-process smoothing + single-process upscaling
  scripts/run_chain.sh --states /path/to/states.nc --smooth-mode single --upscale-mode single
EOF
}

SMOOTH_MODE="parallel"
UPSCALE_MODE="parallel"
OUTPUT_DIR="${DEFAULT_OUTPUT_DIR}"
GRIDLIST="${DEFAULT_GRIDLIST}"
STATES="${DEFAULT_STATES}"
FMFILE=""
BENCHMARK="false"
BENCHMARK_CHUNKS=2
BENCHMARK_LINES=200
SMOOTH_WORKERS=8
SMOOTH_CHUNK_SIZE=100
UPSCALE_WORKERS=8
UPSCALE_CHUNK_LINES=500
PRINT_CONFIG="false"
DRY_RUN="false"
INSPECT="false"
MAPPING_CONFIG=""
MAPPING_PROFILE=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --states) STATES="$2"; shift 2 ;;
    --fm) FMFILE="$2"; shift 2 ;;
    --gridlist) GRIDLIST="$2"; shift 2 ;;
    --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
    --smooth-mode) SMOOTH_MODE="$2"; shift 2 ;;
    --upscale-mode) UPSCALE_MODE="$2"; shift 2 ;;
    --benchmark) BENCHMARK="true"; shift ;;
    --benchmark-chunks) BENCHMARK_CHUNKS="$2"; shift 2 ;;
    --benchmark-lines) BENCHMARK_LINES="$2"; shift 2 ;;
    --smooth-workers) SMOOTH_WORKERS="$2"; shift 2 ;;
    --smooth-chunk-size) SMOOTH_CHUNK_SIZE="$2"; shift 2 ;;
    --upscale-workers) UPSCALE_WORKERS="$2"; shift 2 ;;
    --upscale-chunk-lines) UPSCALE_CHUNK_LINES="$2"; shift 2 ;;
    --mapping-config) MAPPING_CONFIG="$2"; shift 2 ;;
    --mapping-profile) MAPPING_PROFILE="$2"; shift 2 ;;
    --print-config) PRINT_CONFIG="true"; shift ;;
    --dry-run) DRY_RUN="true"; shift ;;
    --inspect) INSPECT="true"; shift ;;
    -h|--help) show_help; exit 0 ;;
    *) echo "Unknown argument: $1"; show_help; exit 1 ;;
  esac
done

# Build pass-through args for the upscaling script
UPSCALE_MAPPING_ARGS=""
[[ -n "${MAPPING_CONFIG}" ]]  && UPSCALE_MAPPING_ARGS+=" --mapping-config \"${MAPPING_CONFIG}\""
[[ -n "${MAPPING_PROFILE}" ]] && UPSCALE_MAPPING_ARGS+=" --mapping-profile \"${MAPPING_PROFILE}\""

if [[ -z "${STATES}" ]]; then
  echo "Error: --states is required"
  show_help
  exit 1
fi
if [[ "${PRINT_CONFIG}" == "true" ]]; then
  echo "Resolved configuration:"
  echo "  STATES=${STATES}"
  echo "  FMFILE=${FMFILE:-<none>}"
  echo "  GRIDLIST=${GRIDLIST}"
  echo "  OUTPUT_DIR=${OUTPUT_DIR}"
  echo "  SMOOTH_MODE=${SMOOTH_MODE}"
  echo "  UPSCALE_MODE=${UPSCALE_MODE}"
  echo "  BENCHMARK=${BENCHMARK}"
  echo "  BENCHMARK_CHUNKS=${BENCHMARK_CHUNKS}"
  echo "  BENCHMARK_LINES=${BENCHMARK_LINES}"
  echo "  SMOOTH_WORKERS=${SMOOTH_WORKERS}"
  echo "  SMOOTH_CHUNK_SIZE=${SMOOTH_CHUNK_SIZE}"
  echo "  UPSCALE_WORKERS=${UPSCALE_WORKERS}"
  echo "  UPSCALE_CHUNK_LINES=${UPSCALE_CHUNK_LINES}"
  exit 0
fi

if [[ "${DRY_RUN}" == "true" ]]; then
  echo "Resolved file paths:"
  echo "  SMOOTH_OUT=${SMOOTH_OUT}"
  echo "  MINI_SMOOTH_OUT=${MINI_SMOOTH_OUT}"
  echo "  MINI_GRIDLIST=${MINI_GRIDLIST}"
  echo "  NETFRAC_OUT=${NETFRAC_OUT}"
  echo "  MINI_NETFRAC_OUT=${MINI_NETFRAC_OUT}"
  echo "  INSPECT=${INSPECT}"
fi

run_cmd() {
  if [[ "${DRY_RUN}" == "true" ]]; then
    echo "[dry-run] $*"
  else
    eval "$@"
  fi
}

mkdir -p "${OUTPUT_DIR}"

SMOOTH_OUT="${OUTPUT_DIR}/hildaplus_smoothed.nc"
MINI_SMOOTH_OUT="${OUTPUT_DIR}/mini_smoothed.nc"
MINI_GRIDLIST="${OUTPUT_DIR}/mini_gridlist.txt"
NETFRAC_OUT="${OUTPUT_DIR}/hildaplus_netfrac_1901_2020.txt"
MINI_NETFRAC_OUT="${OUTPUT_DIR}/mini_netfrac.txt"

echo "==> Smoothing (${SMOOTH_MODE})"
if [[ "${SMOOTH_MODE}" == "parallel" ]]; then
  SMOOTH_SCRIPT="${SCRIPT_DIR}/hilda_smoothing/smoothing_local_parallel.py"
  if [[ "${BENCHMARK}" == "true" ]]; then
    run_cmd "python \"${SMOOTH_SCRIPT}\" \
      --states \"${STATES}\" \
      ${FMFILE:+--fm \"${FMFILE}\"} \
      --output \"${MINI_SMOOTH_OUT}\" \
      --workers \"${SMOOTH_WORKERS}\" \
      --chunk-size \"${SMOOTH_CHUNK_SIZE}\" \
      --benchmark-chunks \"${BENCHMARK_CHUNKS}\""
  else
    run_cmd "python \"${SMOOTH_SCRIPT}\" \
      --states \"${STATES}\" \
      ${FMFILE:+--fm \"${FMFILE}\"} \
      --output \"${SMOOTH_OUT}\" \
      --workers \"${SMOOTH_WORKERS}\" \
      --chunk-size \"${SMOOTH_CHUNK_SIZE}\""
  fi
else
  SMOOTH_SCRIPT="${SCRIPT_DIR}/hilda_smoothing/smoothing_local.py"
  if [[ "${BENCHMARK}" == "true" ]]; then
    run_cmd "python \"${SMOOTH_SCRIPT}\" \
      --states \"${STATES}\" \
      ${FMFILE:+--fm \"${FMFILE}\"} \
      --output \"${MINI_SMOOTH_OUT}\" \
      --chunk-size \"${SMOOTH_CHUNK_SIZE}\" \
      --benchmark-chunks \"${BENCHMARK_CHUNKS}\""
  else
    run_cmd "python \"${SMOOTH_SCRIPT}\" \
      --states \"${STATES}\" \
      ${FMFILE:+--fm \"${FMFILE}\"} \
      --output \"${SMOOTH_OUT}\" \
      --chunk-size \"${SMOOTH_CHUNK_SIZE}\""
  fi
fi

if [[ "${BENCHMARK}" == "true" ]]; then
  echo "==> Creating mini gridlist"
  run_cmd "python \"${SCRIPT_DIR}/hilda_smoothing/make_minigridlist.py\" \
    --netcdf \"${MINI_SMOOTH_OUT}\" \
    --gridlist \"${GRIDLIST}\" \
    --output \"${MINI_GRIDLIST}\""
fi

echo "==> Upscaling (${UPSCALE_MODE})"
if [[ "${UPSCALE_MODE}" == "parallel" ]]; then
  UPSCALE_SCRIPT="${SCRIPT_DIR}/hildaplus-upscaling/hildap_tables_netfrac_v3_parallel.py"
  if [[ "${BENCHMARK}" == "true" ]]; then
    run_cmd "python \"${UPSCALE_SCRIPT}\" \
      --datafile \"${MINI_SMOOTH_OUT}\" \
      --gridlist \"${MINI_GRIDLIST}\" \
      ${FMFILE:+--fmfile \"${FMFILE}\"} \
      --output \"${MINI_NETFRAC_OUT}\" \
      --workers \"${UPSCALE_WORKERS}\" \
      --chunk-lines \"${UPSCALE_CHUNK_LINES}\" \
      --benchmark-lines \"${BENCHMARK_LINES}\"${UPSCALE_MAPPING_ARGS}"
  else
    run_cmd "python \"${UPSCALE_SCRIPT}\" \
      --datafile \"${SMOOTH_OUT}\" \
      --gridlist \"${GRIDLIST}\" \
      ${FMFILE:+--fmfile \"${FMFILE}\"} \
      --output \"${NETFRAC_OUT}\" \
      --workers \"${UPSCALE_WORKERS}\" \
      --chunk-lines \"${UPSCALE_CHUNK_LINES}\"${UPSCALE_MAPPING_ARGS}"
  fi
else
  UPSCALE_SCRIPT="${SCRIPT_DIR}/hildaplus-upscaling/hildap_tables_netfrac_v3.py"
  if [[ "${BENCHMARK}" == "true" ]]; then
    run_cmd "python \"${UPSCALE_SCRIPT}\" \
      --datafile \"${MINI_SMOOTH_OUT}\" \
      --gridlist \"${MINI_GRIDLIST}\" \
      ${FMFILE:+--fmfile \"${FMFILE}\"} \
      --output \"${MINI_NETFRAC_OUT}\" \
      --benchmark-lines \"${BENCHMARK_LINES}\"${UPSCALE_MAPPING_ARGS}"
  else
    run_cmd "python \"${UPSCALE_SCRIPT}\" \
      --datafile \"${SMOOTH_OUT}\" \
      --gridlist \"${GRIDLIST}\" \
      ${FMFILE:+--fmfile \"${FMFILE}\"} \
      --output \"${NETFRAC_OUT}\"${UPSCALE_MAPPING_ARGS}"
  fi
fi

if [[ "${INSPECT}" == "true" ]]; then
  echo "==> Inspecting smoothed output"
  SMOOTH_FILE="${SMOOTH_OUT}"
  if [[ "${BENCHMARK}" == "true" ]]; then
    SMOOTH_FILE="${MINI_SMOOTH_OUT}"
  fi
  run_cmd "python \"${SCRIPT_DIR}/hilda_smoothing/inspect_smoothed.py\" \
    --input \"${SMOOTH_FILE}\" \
    --output \"${OUTPUT_DIR}/inspect_smoothed.txt\""

  echo "==> Inspecting netfrac output"
  NETFRAC_FILE="${NETFRAC_OUT}"
  if [[ "${BENCHMARK}" == "true" ]]; then
    NETFRAC_FILE="${MINI_NETFRAC_OUT}"
  fi
  run_cmd "python \"${SCRIPT_DIR}/hildaplus-upscaling/inspect_netfrac_forestfrac.py\" \
    --input \"${NETFRAC_FILE}\" \
    --output \"${OUTPUT_DIR}/inspect_netfrac_forestfrac.txt\""
fi

echo "==> Done"
echo "Outputs written to: ${OUTPUT_DIR}"

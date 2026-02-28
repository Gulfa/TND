#!/usr/bin/env bash
# run_in_docker.sh — Execute VE simulation inside the FHI R container.
#
# Usage:
#   bash run_in_docker.sh [container_name]
#
# Default container name: fhi_r
# The script mounts nothing; it relies on the host path being visible
# inside the container (e.g. via a bind mount configured at container start).

set -euo pipefail

CONTAINER="${1:-fhi}"
SCRIPT_DIR="/fhi/TND"

echo "==> Target container: ${CONTAINER}"
echo "==> Working dir:      ${SCRIPT_DIR}"

# ---------------------------------------------------------------------------
# Optional first-time package installation.
# Run once manually if packages are missing:
#
#   docker exec "${CONTAINER}" Rscript -e \
#     "install.packages(c('odin','dplyr','tidyr','ggplot2','scales','purrr'), \
#      repos='https://cloud.r-project.org')"
# ---------------------------------------------------------------------------

echo "==> Running main.R ..."
docker exec \
  -w "${SCRIPT_DIR}" \
  "${CONTAINER}" \
  Rscript main.R

echo "==> Done. Check ${SCRIPT_DIR}/outputs/ for results."

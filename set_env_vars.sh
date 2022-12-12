#!/bin/bash
#
# Set environment variables used by other scripts
#
function set_conda_env_prefixes() {
  local ORIG_DIR="$(pwd)" || return 1

  local REL_SCRIPT_DIR="$(dirname ${BASH_SOURCE[0]})" || return 1
  cd "${REL_SCRIPT_DIR}" || return 1
  SCRIPT_DIR="$(pwd)" || return 1

  cd "${ORIG_DIR}" || return 1

  CONDA_ENV_PREFIX_2="${SCRIPT_DIR}/conda_env_2"
  CONDA_ENV_PREFIX_3="${SCRIPT_DIR}/conda_env_3"
}

function main() {
  # need to use the setup that conda init writes to .bashrc
  source "${HOME}/.bashrc" || return 1
  set_conda_env_prefixes || return 1
}

main "$@"

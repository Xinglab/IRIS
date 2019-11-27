#!/bin/bash
#
# Provide functions for working with conda
# Usage:
# * `source conda.sh`
# * `conda::create_env_with_name_and_python_version {env_name} {python_version}`
#   + example: `conda::create_env_with_name_and_python_version my-conda-env 3.6`
# * `conda::activate_env {env_name}`
#   + example: `conda::activate_env my-conda-env`
# * `conda::deactivate_env`
#
# Assumes that conda is installed
# https://docs.conda.io/en/latest/miniconda.html
# https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh

function conda::create_env_with_name_and_python_version() {
  local ERROR_PREFIX="error in conda::create_env_with_name_and_python_version()"

  local ENV_NAME="$1"
  local PYTHON_VERSION="$2"

  local FOUND_COUNT=$(conda info --envs | grep "^${ENV_NAME} .*/${ENV_NAME}$" | wc -l)
  if [[ "$?" -ne 0 ]]; then
    echo "${ERROR_PREFIX}: checking conda envs" >&2
    return 1
  fi

  if [[ "${FOUND_COUNT}" -eq 1 ]]; then
    echo "using existing ${ENV_NAME} conda environment"
    return 0
  fi

  echo "creating new conda environment: ${ENV_NAME} python=${PYTHON_VERSION}"
  conda create --name "${ENV_NAME}" python="${PYTHON_VERSION}"
  if [[ "$?" -ne 0 ]]; then
    echo "${ERROR_PREFIX}: creating env" >&2
    return 1
  fi
}
export -f conda::create_env_with_name_and_python_version

function conda::activate_env() {
  conda activate "$1"
}
export -f conda::activate_env

function conda::deactivate_env() {
  conda deactivate
}
export -f conda::deactivate_env

function main() {
  # need to use the setup that conda init writes to .bashrc
  source ${HOME}/.bashrc || return 1
}

main "$@"

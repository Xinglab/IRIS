#!/bin/bash

SCRIPT_DIR=$(pwd)
export PATH=${PATH}:${SCRIPT_DIR}/bedtools/bedtools2/bin

# need to use the setup that conda init writes to .bashrc
source ${HOME}/.bashrc || exit 1

CONDA_ENV_NAME_2='iris2'
CONDA_ENV_NAME_3='iris3'
conda activate ${CONDA_ENV_NAME_2} || exit 1

echo
echo "Formatting"

IRIS formatting ./example/SJ_matrices/matrices.txt ./example/SJ_matrices/samples.txt -s 2 -d ./IRIS_data/db/ -t SE -n Glioma_test || exit 1

echo
echo "Screening"

IRIS screening ./example/Test.para -t -o ../Glioma_test/screening || exit 1

echo
echo "prediction"

TEMP_F_NAME=$(mktemp)
# TODO --iedb-local should be required=True
IRIS prediction ../Glioma_test/screening -c 5 -m ./example/HLA_types/hla_types.list -p ./example/Test.para --iedb-local ./IEDB/mhc_i/src > ${TEMP_F_NAME} || exit 1
cat ${TEMP_F_NAME}

QSUB_CMDS=()
while read LINE
do
    GREP_RES=$(echo ${LINE} | grep '^qsub.*\.sh$')
    if [[ -n ${GREP_RES} ]]
    then
        QSUB_CMDS+=("${LINE}")
    fi
done < ${TEMP_F_NAME}

rm ${TEMP_F_NAME}

if [[ ${#QSUB_CMDS[@]} == 0 ]]
then
    echo "the prediction step did not produce any qsub commands"
    exit 1
fi

echo
echo "executing qsub commands"

# Need to use Python3 for submit_and_wait.py.
# Also need to have the Python2 conda environment to run IRIS.
# Get the Python3 path and then go back to the Python2 environment

conda deactivate || exit 1
conda activate ${CONDA_ENV_NAME_3} || exit 1

PYTHON_3_EXECUTABLE=$(which python) || exit 1

conda deactivate || exit 1
conda activate ${CONDA_ENV_NAME_2} || exit 1

SUBMIT_AND_WAIT_CMD="${PYTHON_3_EXECUTABLE} submit_and_wait.py"
TEMP_F_NAME=$(mktemp)

for QSUB_CMD in "${QSUB_CMDS[@]}"
do
    echo "${QSUB_CMD}" >> ${TEMP_F_NAME}
done

echo "execute: ${SUBMIT_AND_WAIT_CMD}"
echo "with qsub commands:"
cat ${TEMP_F_NAME}

${SUBMIT_AND_WAIT_CMD} ${TEMP_F_NAME} || exit 1

rm ${TEMP_F_NAME}

echo
echo "epitope post"

# TODO -e is actually not required?
IRIS epitope_post -p ./example/Test.para -o ../Glioma_test/screening -m ./example/HLA_types/hla_patient.tsv || exit 1

conda deactivate || exit 1

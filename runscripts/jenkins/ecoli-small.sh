set -e

module load wcEcoli/python3
pyenv local wcEcoli3

# Prevent differences when run with more than one thread using the current
# OpenBLAS library. This could be removed if a test shows the current library
# is not thread sensitive.
export OPENBLAS_NUM_THREADS=1

make clean compile

PYTHONPATH=$PWD:$PYTHONPATH pytest --cov=wholecell --cov-report xml \
    --junitxml=unittests.xml

runscripts/debug/mypy.sh

module load wcEcoli/python3
pyenv local wcEcoli3

# Prevent differences when run with more than one thread using the current
# OpenBLAS library. This could be removed if a test shows the current library
# is not thread sensitive.
export OPENBLAS_NUM_THREADS=1

make clean
make compile

set -e

# Running it this way prints all timing measurements:
python -m wholecell.tests.utils.test_library_performance

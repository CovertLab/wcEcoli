set -e

module load wcEcoli/python3

WCECOLI_PYENV=wcEcoli3
pyenv local ${WCECOLI_PYENV}

make clean compile

PYTHONPATH=$PWD:$PYTHONPATH pytest --cov=wholecell --cov-report xml \
    --junitxml=unittests.xml

(export PYENV_VERSION="mypy:${WCECOLI_PYENV}"; echo ---Running mypy---; mypy)

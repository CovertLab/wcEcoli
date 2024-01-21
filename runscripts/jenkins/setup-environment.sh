set -e

echo se1
export PYTHONPATH=$PWD
echo se2
module load wcEcoli/python3

if [ -d "${PYENV_ROOT}" ]; then
    echo se3
    export PATH="${PYENV_ROOT}/bin:${PATH}"
    echo se4
    eval "$(pyenv init -)"
    echo se5
    eval "$(pyenv virtualenv-init -)"
fi

### Edit this line to make this branch use another pyenv like wcEcoli3-staging
echo se6
pyenv local wcEcoli3
echo se7
pyenv activate

echo se8
make clean compile
echo se9

set -e

export PYTHONPATH=$PWD
module load wcEcoli/python3

### -------------------------------------------------------------------
### Edit this line to make the PR build use another pyenv like wcEcoli3-staging.
### Revert it to `wcEcoli3` before merging the PR into master.
### -------------------------------------------------------------------
pyenv local wcEcoli3

make clean compile

set -e

### ---------------------------------------------------------------------------
### Edit the pyenv in setup-environment.sh to make the PR build temporarily use
### another pyenv for testing (eg. wcEcoli3-staging). Revert the change before
### merging the PR into master to prevent changing it for other Jenkins builds.
### ---------------------------------------------------------------------------
echo eprt1
source runscripts/jenkins/setup-environment.sh

echo eprt2
echo eprt2 pyenv: "$(pyenv version)"
pip list | grep 'numpy\|scipy'
runscripts/debug/numpy_benchmark.py

echo eprt3
pytest --cov=wholecell --cov-report xml --junitxml=unittests.xml

echo eprt4
echo eprt4 pyenv: "$(pyenv version)"
runscripts/debug/numpy_benchmark.py

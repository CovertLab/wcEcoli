set -e

### ---------------------------------------------------------------------------
### Edit the pyenv in setup-environment.sh to make the PR build temporarily use
### another pyenv for testing (eg. wcEcoli3-staging). Revert the change before
### merging the PR into master to prevent changing it for other Jenkins builds.
### ---------------------------------------------------------------------------
echo eprt1
source runscripts/jenkins/setup-environment.sh

echo eprt2
echo pyenv version: $(pyenv version)

echo eprt3
pytest --cov=wholecell --cov-report xml --junitxml=unittests.xml
echo eprt4

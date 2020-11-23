set -e

sh runscripts/jenkins/setup-environment.sh
PYTHONPATH=$PWD:$PYTHONPATH pytest --cov=wholecell --cov-report xml \
    --junitxml=unittests.xml

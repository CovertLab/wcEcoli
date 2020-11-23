set -e

sh runscripts/jenkins/setup-environment.sh
runscripts/debug/mypy.sh

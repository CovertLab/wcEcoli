#!/bin/sh
# Use Google Cloud Build servers to build a personalized "${ID}-wcm-code" Docker
# Image and store it in a Google Artifact Registry.
#
# TODO: Require and use config artifacts/location, artifacts/repository?
#       Otherwise check that $REGION is not empty.
#
# COMMAND LINE ARGUMENTS:
#   ARG1: Distinguishing ID prefix for the "wcm-code" Docker Image tag;
#     defaults to "${USER}".
#   ARG2: Docker tag for the wcm-runtime Image in Artifact Registry to build FROM;
#     defaults to "${ID}-wcm-runtime".
#     The named Docker Image must already exist in the Artifact Registry.
#
# ASSUMES: The current working dir is the wcEcoli/ project root.

set -eu

ID="${1:-$USER}"

REGION="$(gcloud config get compute/region)"
WCM_RUNTIME="${2:-${ID}-wcm-runtime}"
WCM_CODE="${ID}-wcm-code"
GIT_HASH=$(git rev-parse HEAD)
GIT_BRANCH=$(git symbolic-ref --short HEAD)
TIMESTAMP=$(date '+%Y%m%d.%H%M%S')

mkdir -p source-info
git diff HEAD > source-info/git_diff.txt

echo "=== Cloud-building WCM code Docker Image ${WCM_CODE} on ${WCM_RUNTIME} ==="
echo "=== git hash ${GIT_HASH}, git branch ${GIT_BRANCH} ==="

# This needs a config file to identify the project files to upload and the
# Dockerfile to run.
# --region sets the Cloud Build region.
gcloud builds submit --timeout=15m --config config-build-2-wcm-code.json \
    --region="${REGION}" \
    --substitutions="_WCM_RUNTIME=${WCM_RUNTIME},_WCM_CODE=${WCM_CODE},\
_GIT_HASH=${GIT_HASH},_GIT_BRANCH=${GIT_BRANCH},_TIMESTAMP=${TIMESTAMP}"

rm source-info/git_diff.txt

#!/bin/sh
# Use Google Cloud Build servers to build a personalized "${ID}-wcm-runtime"
# Docker Image and store it in a Google Artifact Registry.
#
# COMMAND LINE ARGUMENTS:
#   ARG1: Distinguishing ID prefix for the "${ID}-wcm-runtime" tag for the
#     Docker Image to build; defaults to "${USER}". Must be in lower case.
#
# ASSUMES: The current working dir is the wcEcoli/ project root.

set -eu

ID="${1:-$USER}"

PROJECT="$(gcloud config get core/project)"
REGION="$(gcloud config get compute/region)"
WCM_RUNTIME="${ID}-wcm-runtime"
TAG="${REGION}-docker.pkg.dev/${PROJECT}/wcm/${WCM_RUNTIME}"

echo "=== Cloud-building WCM runtime Docker Image: ${TAG} ==="

# This needs only one payload file so copy it in rather than using a config at
# the project root which would upload the entire project.
# --region sets the Cloud Build region while TAG sets the artifact region et al.
cp requirements.txt cloud/docker/runtime/
gcloud builds submit --timeout=3h --region="${REGION}" \
    --tag="${TAG}" cloud/docker/runtime/
rm cloud/docker/runtime/requirements.txt

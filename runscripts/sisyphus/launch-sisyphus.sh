#!/usr/bin/env bash

set -eu

NAME=$1
FULL_NAME=sisyphus-$NAME
PROJECT=allen-discovery-center-mcovert

gcloud compute \
       --project=$PROJECT \
       instances create $FULL_NAME \
       --zone=us-west1-b \
       --machine-type=n1-standard-2 \
       --subnet=default \
       --network-tier=PREMIUM \
       --maintenance-policy=MIGRATE \
       --service-account=441871726775-compute@developer.gserviceaccount.com \
       --scopes=https://www.googleapis.com/auth/devstorage.read_only,https://www.googleapis.com/auth/logging.write,https://www.googleapis.com/auth/monitoring.write,https://www.googleapis.com/auth/servicecontrol,https://www.googleapis.com/auth/service.management.readonly,https://www.googleapis.com/auth/trace.append \
       --image=sisyphus-base \
       --image-project=$PROJECT \
       --boot-disk-size=200GB \
       --boot-disk-type=pd-standard \
       --boot-disk-device-name=$FULL_NAME \
       --description='sisyphus worker'

# Container image #3: wcm-parca.
# This Dockerfile builds a container image with the wcEcoli Parca (parameter
# calculator) data (ready to run simulations), layered on the wcm-code container
# image.
#
# To build locally in the wcEcoli/ project root directory:
#
#     > docker build -f cloud/docker/parca/Dockerfile -t parca .
#
# To build using the Cloud Build service and store in the Container Registry, in
# the project root directory:
#
#     > gcloud builds submit --timeout=90m --config cloudbuild3-wcm-parca.json
#
# After building locally you can start up a new container from the image:
#
#     > docker run --name parcatest -it -d parca
#
# Once the container exists you can `exec` a shell inside it:
#
#     > docker exec -it parcatest /bin/sh
#
# You will be dropped into an `sh` shell where you can execute commands:
#
#     # nosetests
#
# If this succeeds you should be good to go, e.g.:
#
#     # python runscripts/manual/runSim.py


FROM gcr.io/allen-discovery-center-mcovert/wcm-code:latest

RUN python runscripts/manual/runParca.py

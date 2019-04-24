# Container image #3: wcm-parca.
# This Dockerfile builds a container image with the wcEcoli Parca output (ready
# to run simulations), layered on the wcm-code image.
#
# To build locally from the wcEcoli/ project root directory:
#
#     > docker build -f cloud/docker/parca/Dockerfile -t wcm-parca .
#
# To build using the Cloud Build service and store in the Container Registry,
# from this directory:
#
#     > gcloud builds submit --timeout=90m --tag gcr.io/allen-discovery-center-mcovert/wcm-parca .
#
# After building locally you can start up a new container from the image:
#
#     > docker run --name parcatest -it parca
#
# It will start a shell where you can execute commands:
#
#     # nosetests
#
# If this succeeds you should be good to go, e.g.:
#
#     # python runscripts/manual/runSim.py


FROM gcr.io/allen-discovery-center-mcovert/wcm-code:latest

RUN python runscripts/manual/runParca.py

CMD ["/bin/bash"]

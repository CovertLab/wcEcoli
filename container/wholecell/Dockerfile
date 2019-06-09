# Container image #2: wcm-code.
# This Dockerfile builds a container image with the wcEcoli whole cell model
# code, layered on the wcm-runtime container image.
#
# To build this image from the wcEcoli/ project root directory:
#
#     > docker build -f container/wholecell/Dockerfile -t wcm-code .
#
# After building you can start up a new container from the image:
#
#     > docker run --name wholecelltest -it wcm-code
#
# It will start a shell where you can execute commands:
#
#     # nosetests
#
# If this succeeds you can start running WCM code in the container, e.g.:
#
#     # python runscripts/manual/runFitter.py

FROM wcm-runtime:latest

COPY . /wcEcoli
WORKDIR /wcEcoli

RUN make clean compile
ENV PYTHONPATH=/wcEcoli

CMD ["/bin/bash"]

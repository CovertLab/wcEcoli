# Container image #1: wcm-runtime.
# This Dockerfile builds the runtime environment for the whole cell model.
#
# To build this image locally from the wcEcoli/ project root directory:
#
#     > docker build -f cloud/docker/runtime/Dockerfile -t ${USER}-wcm-runtime .
#
# (To build using the Cloud Build service and store in the Container Registry,
# run `cloud/build.sh`.)
#
# Add option `--build-arg from=ABC` to build from a different base image "ABC"
# but DO NOT USE an alpine base since the simulation math comes out different!
# See https://pythonspeed.com/articles/alpine-docker-python/ for more reasons.
ARG from=python:3.8.3
FROM ${from}

# Add option `--build-arg NO_AVX2=1` to make OpenBLAS 0.3.6+ not use the AVX2
# hardware so OpenBLAS can pass its build self-test and run properly in macOS
# Docker.
#
# NO_AVX2=1 is a fine workaround BUT it costs (reportedly) 20-30% of BLAS
# runtime performance. Also it needs to be set when building an Image that
# might ever run on macOS.
#
# Turn this off once Docker fixes the bug.
# https://github.com/xianyi/OpenBLAS/issues/2244
# https://github.com/docker/for-mac/issues/4576
# https://github.com/machyve/xhyve/issues/171
ARG NO_AVX2=0

RUN apt-get update \
    && apt-get install -y swig gfortran llvm cmake nano

RUN echo "alias ls='ls --color=auto'" >> ~/.bashrc \
    && echo "alias ll='ls -l'" >> ~/.bashrc

# Install openblas (at least 0.3.5, and 0.3.6+ need "NO_AVX2=1" for macOS Docker
# to get correct results) for Theano and maybe numpy/scipy.
# Someday `apt-get install libopenblas-dev=...` will work.
RUN (mkdir -p openblas && cd openblas \
    && curl -SL https://github.com/xianyi/OpenBLAS/archive/v0.3.10.tar.gz | tar -xz \
    && cd OpenBLAS* \
    && make "NO_AVX2=${NO_AVX2}" FC=gfortran \
    && make "NO_AVX2=${NO_AVX2}" PREFIX=/usr install) \
    && rm -r openblas

# Conditionally build numpy and scipy from source to link them to this OpenBLAS
# instead of their embedded copies of OpenBLAS. (See `--no-binary`, below.)
# That can add 50 minutes to the build time!
RUN if [ $NO_AVX2 != 0 ] ; then (echo "[openblas]" \
    && echo "libraries = openblas" \
    && echo "library_dirs = /usr/lib" \
    && echo "include_dirs = /usr/include") > ~/.numpy-site.cfg ; \
    fi

ENV OPENBLAS_NUM_THREADS=1

# Install all the pips within one Docker layer, don't cache the downloads, and
# enable Theano to write into the data dir it gets when running as a user with
# no home dir. (That happens because DockerTask runs tasks inside the Docker
# container with the uid:gid of the user outside the container to get
# read/write/delete access to the tasks's output files.)
COPY py3_requirements.txt /
RUN (b1="" \
    && b2="" \
    && if [ -f ~/.numpy-site.cfg ] ; then b1="--no-binary"; b2="numpy,scipy"; fi \
    && pip install --no-cache-dir --upgrade pip setuptools wheel \
    && pip install --no-cache-dir numpy==1.19.0 $b1 $b2 \
    && pip install --no-cache-dir -r py3_requirements.txt $b1 $b2 \
    && umask 000 && mkdir -p /.theano)

CMD ["/bin/bash"]

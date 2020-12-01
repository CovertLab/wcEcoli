#! /usr/bin/env sh

set -eu

WC_MONGO_DB=$1

mkdir -p /scratch/PI/mcovert/jenkins/fireworks/logs/launchpad
echo "authsource: admin" > my_launchpad.yaml
echo "host: mongodb+srv://${WC_MONGO_USER}:${WC_MONGO_PW}@${WC_MONGO_CLUSTER}/${WC_MONGO_DB}?retryWrites=true&w=majority" >> my_launchpad.yaml
echo "logdir: /scratch/PI/mcovert/jenkins/fireworks/logs/launchpad" >> my_launchpad.yaml
echo "mongoclient_kwargs: {}" >> my_launchpad.yaml
echo "name: null" >> my_launchpad.yaml
echo "password: null" >> my_launchpad.yaml
echo "port: null" >> my_launchpad.yaml
echo "ssl: false" >> my_launchpad.yaml
echo "ssl_ca_certs: null" >> my_launchpad.yaml
echo "ssl_certfile: null" >> my_launchpad.yaml
echo "ssl_keyfile: null" >> my_launchpad.yaml
echo "ssl_pem_passphrase: null" >> my_launchpad.yaml
echo "strm_lvl: INFO" >> my_launchpad.yaml
echo "uri_mode: true" >> my_launchpad.yaml
echo "user_indices: []" >> my_launchpad.yaml
echo "username: null" >> my_launchpad.yaml
echo "wf_user_indices: []" >> my_launchpad.yaml

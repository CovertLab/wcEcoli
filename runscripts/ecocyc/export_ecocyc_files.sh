#! /usr/bin/env bash
set -e

mkdir -p out/temp
if [ $(find "$1" -name 'wcm*' | wc -l) -gt 0 ]; then
  cp -p $(find "$1" -name 'wcm*') out/temp
fi

tar -zcf out/wcm-files.tar.gz -C out/temp .
rm -r out/temp

rclone copy out/wcm-files.tar.gz remote:

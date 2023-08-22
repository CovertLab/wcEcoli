#!/usr/bin/env bash
set -e

module load system rclone

mkdir -p out/temp
if [ $(find "$1" -name 'wcm*' | wc -l) -gt 0 ]; then
  cp -p $(find "$1" -name 'wcm*') out/temp
fi

for file in out/temp/*; do
  filename=$(basename "$file" | sed 's/\(.*\)\..*/\1/')
  MEDIA="${filename##*_}"
  file_basename=$(basename "${file//_$MEDIA/}")
  mkdir -p out/temp/"$MEDIA"/ && \
  mv "$file" "$_$file_basename"
done

tar -zcf out/wcm-files.tar.gz -C out/temp .
rm -r out/temp

rclone copy out/wcm-files.tar.gz remote:

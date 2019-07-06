#!/bin/sh

bowtie_root="${PWD%%bowtie*}/bowtie"

docker run -t -i --rm \
  -v $bowtie_root:/io \
  phusion/holy-build-box-64:latest \
  bash /io/scripts/bowtie-hbb.sh

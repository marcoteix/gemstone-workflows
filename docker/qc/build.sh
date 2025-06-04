#!/bin/bash

set -e

VERSION="1.0.0"

cd docker/qc
echo Building the Docker image

# Build for amd64; base image does not support arm builds
docker buildx build \
    --platform linux/amd64 \
    --tag marcoteix/gemstone-qc:$VERSION .

echo Pushing to Docker Hub
docker push marcoteix/gemstone-qc:$VERSION

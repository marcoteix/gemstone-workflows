#!/bin/bash

VERSION="1.0.0"

cd docker/qc
echo Building the Docker image
docker build -t marcoteix/gemstone-qc:$VERSION .
echo Pushing to Docker Hub
docker push marcoteix/gemstone-qc:$VERSION

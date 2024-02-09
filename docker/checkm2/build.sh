#!/bin/bash

VERSION="1.0.0"

cd docker/checkm2
echo Building the Docker image
docker build -t marcoteix/checkm2:$VERSION .
echo Pushing to Docker Hub
docker push marcoteix/checkm2:$VERSION

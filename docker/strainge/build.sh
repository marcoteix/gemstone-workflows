#!/bin/bash

VERSION="1.0.1"

cd docker/strainge
echo Building the Docker image
docker build -t marcoteix/strainge:$VERSION .
echo Pushing to Docker Hub
docker push marcoteix/strainge:$VERSION
echo Done!

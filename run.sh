#!/bin/bash

# define names
build_name="react"
container_name="react-container"
folder_name="react"  # see Dockerfile

# build
docker build . -t $build_name --no-cache

# run
docker run --name $container_name $build_name julia runEverything.jl

# restart (if the container was stopped)
# docker start $container_name
# docker exec -it $container_name julia runEverything.jl

# copy results
docker cp $container_name:/$folder_name/plots/ .
docker cp $container_name:/$folder_name/results/ .

# clean up
# rm -Rf $folder_name
# docker rm --force $container_name
# docker image rm --force $build_name

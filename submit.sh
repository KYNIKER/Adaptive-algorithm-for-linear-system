
build_name="dmg-project"
container_name="dmg-project-container"

docker build . -t $build_name --no-cache

docker run --name $container_name $build_name julia runEverything.jl


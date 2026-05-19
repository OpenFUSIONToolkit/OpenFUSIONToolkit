# How to use the OFT container
Examples below use `podman`, but you should be able to use `docker` with exactly the same flags and arguments.

## To build
```shell
./docker/build_container.sh docker/Dockerfile.OFT_run_container oft_run
```

## To run a Jupyter lab server in the container
Replace `HOST_PORT` below with a suitable open port on your host system and then replace `8888` with this number in the URL
given to connect to the Jupyter server
```shell
podman run --userns=keep-id -v $(pwd):/home/oft_user -p HOST_PORT:8888 oft_run
```
or with docker
```shell
docker run -u $(id -u):$(id -g) -v $(pwd):/home/oft_user -p HOST_PORT:8888 oft_run
```

## To run a Python script
**Note:** Script must be in the current working directory or a child directory
```shell
podman run --userns=keep-id -v $(pwd):/home/oft_user oft_run script_path.py
```
or with docker
```shell
docker run -u $(id -u):$(id -g) -v $(pwd):/home/oft_user oft_run script_path.py
```

## To run a Jupyter notebook
**Note:** Notebook must be in the current working directory or a child directory
```shell
podman run --userns=keep-id -v $(pwd):/home/oft_user oft_run -j notebook_path.ipynb
```
or with docker
```shell
docker run -u $(id -u):$(id -g) -v $(pwd):/home/oft_user oft_run -j notebook_path.ipynb
```

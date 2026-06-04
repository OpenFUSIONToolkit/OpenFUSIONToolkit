# How to use the OFT container
Examples below use `podman` and `docker` below, be careful about flag differences depending on the tool you use.

## To build
Replace `IMAGE_NAME` below with the suitable dockerfile suffix.
```shell
./docker/build_container.sh -f docker/Dockerfile.IMAGE_NAME -t oft_run
```

## To run a Python script
**Note:** Script must be in the current working directory or a child directory
```shell
podman run --userns=keep-id -v $(pwd):/work oft_run -s script_path.py
```
or with docker
```shell
docker run -u $(id -u):$(id -g) -v $(pwd):/work oft_run -s script_path.py
```

## To run a Jupyter lab server in the container
Replace `HOST_PORT` below with a suitable open port on your host system and then replace `8888` with this number in the URL
given to connect to the Jupyter server.
**Note:** Only supported with `jupyter` or `extra` image.
```shell
podman run --userns=keep-id -v $(pwd):/work -p HOST_PORT:8888 oft_run -j
```
or with docker
```shell
docker run -u $(id -u):$(id -g) -v $(pwd):/work -p HOST_PORT:8888 oft_run -j
```

## To run a Jupyter notebook
**Note:** Only supported with `jupyter` or `extra` image.
**Note:** Notebook must be in the current working directory or a child directory
```shell
podman run --userns=keep-id -v $(pwd):/work oft_run -n notebook_path.ipynb
```
or with docker
```shell
docker run -u $(id -u):$(id -g) -v $(pwd):/work oft_run -n notebook_path.ipynb
```

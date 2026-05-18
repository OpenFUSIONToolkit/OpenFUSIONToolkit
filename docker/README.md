
## To build
```shell
podman build -f docker/Dockerfile.OFT_run_container -t oft_run .
```

## To run
```shell
podman run -p 8889:8888 -u $(id -u):$(id -g) -v $(pwd):/home/oft_user oft_run
```
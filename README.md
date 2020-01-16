# PRYMETIME

## Docker

### Build Image

```shell
docker build --tag prymetime:0.1 .
```

### Run container interactively

```shell
docker run -it --rm prymetime:0.1
```

### Run container with some data

Mount a directory with the `-v` flag. The directory before the `:`
must be an absolute path to a file or directory, and the directory
after the `:` is where it will be mounted inside the container.

```shell
docker run -it --rm \
    -v /path/to/input_dir:/input \
    -v /path/to/output_dir:/output \
    prymetime:0.1
```

# Run the docker container interactively

The entrypoint script can be overridden for debugging using the
`--entrypoint` argument to docker run. Using `/bin/bash` as the
entrypoint starts an interactive shell when the docker image is
run. Here is an example:

```shell
docker run -it --rm \
    -v $(realpath ../wpi-data):/input \
    -v $(realpath output):/output \
    --entrypoint /bin/bash \
    prymetime:0.1
```

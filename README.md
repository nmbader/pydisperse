# Build software for dispersion curve tracing
by [Milad Bader](mailto:nmbader@sep.stanford.edu)

## Docker

Build the docker image (it should take about 20 seconds)
```
docker build -f Dockerfile -t dispersion .
```

Then run a container
```
docker run -it -p 8080:8080 dispersion
```

or with a local volume mounted (for development purposes)
```
docker run -it -v local_host_directory:target_container_directory -p 8080:8080 dispersion
```

By default a bash shell will be opened at /home inside the container.
Run jupyter notebook from within the container
```
jupyter notebook --ip 0.0.0.0 --port 8080 --no-browser --allow-root &
```

Open the browser at *localhost:8080/​* and use the printed token above to authenticate.
Go to /home/notebooks to run some examples.
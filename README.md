# Software for dispersion curve tracing
by [Milad Bader](mailto:nmbader@sep.stanford.edu)

## Description

Compute dispersion functions and trace dispersion curves in elastic isotropic layered media. Currently P-SV modal solutions are covered. This includes the following configurations:
- Surface waves 
- Lamb waves (known also as plate waves)
- Guided waves in a perfect waveguide (infinitly stiff boundaries)
- Guided waves in an embedded waveguide (waveguide sandwiched between two half spaces)

In the last configuration, both normal and leaky modes are considered. For leaky modes, attenuation due to leakage into the half-spaces is also computed for each mode.

Check the pdf documents for more details about the method and algebraic derivation.

## Installation with Docker

Build the docker image (it should take about 2 minutes)
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
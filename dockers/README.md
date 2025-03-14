## Running in a docker

Docker container allows to run opm inside an encapsulated environment. To use it either pull the public [docker image](https://hub.docker.com/r/jafranc/opm-u22-spe11spe) 

```
	docker pull jafranc/opm-u22-spe11csp
```

or build it from the source _omp-spe11csp-ubuntu.Dockerfile_ located at _./dockers_ with command:

```
	cd dockers
	docker build --build-arg IMG=ubuntu --build-arg VERSION=jammy -t <tag-docker> -f omp-spe11csp-ubuntu.Dockerfile .
```

Then using `docker images`, the created or pulled image should be visible. Running the data generator and simulation for a case is then done by

```
	docker run [-e SPE_CASE=<a,b or c>  -e INPUT_PATH=<docker-path-to-input>] -v <host-path-to-result>:/opt/spe11csp/output_csp11<a,b or c> -v <host-path-to-input>:<docker-path-to-input> <tag-docker> 

```

Then you can inspect what has been produced at \<host-path-to-result\>. \<tag-docker\> being the tag you choose when building the image or if pulled from _docker://_ , _jafranc/opm-u22-spe11csp:latest_. Eventually \<a,b or c\> will set a variable in the docker container generating and running the simulation either for csp SPE 11th case a (2D surface conditions), case b (2D reservoir conditions) or c (3D extruded reservoir conditions). 

The _SPE_CASE_ env variable is defaulted to _a_, the _INPUT_PATH_ to _/opt/spe11csp/examples/hello_world/_ so that it runs the simplest example. Eventually it will run the following command:

```bash

sh -c /usr/local/bin/pyopmspe11 -i ${INPUT_PATH}/spe11${SPE_CASE}.txt -o output_csp11${SPE_CASE}

```

Hence the input file name has to be formatted as spe11a.txt (or spe11b.txt or spe11c.txt).

```bash
	docker run -e SPE_CASE=b -v /work/docker/opm/spe11/b/:/opt/spe11csp/output_csp11b  jafranc/opm-u22-spe11csp:latest

```


is a concrete example running _spe11b.txt_ from _hello_world_ and dumping the results files in the bound directory at _/work/docker/opm/spe11/b/_.


Results can be visualized using [ResInsight](https://resinsight.org/) or other viewers.

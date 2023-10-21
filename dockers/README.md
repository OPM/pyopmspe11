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
	docker run -e SPE_CASE=<a,b or c> -v <host-path-to-result>:/opt/spe11csp/output_csp11<a,b or c> <tag-docker> 

```

Then you can inspect what has been produced at \<host-path-to-result\>. \<tag-docker\> being the tag you choose when building the image or if pulled from _docker://_ , _jafranc/opm-u22-spe11csp:latest_. Eventually \<a,b or c\> will set a variable in the docker container generating and running the simulation either for csp SPE 11th case a (2D surface conditions), case b (2D reservoir conditions) or c (3D extruded reservoir conditions). 

Results can be visualized using [ResInsight](https://resinsight.org/) or other viewers.
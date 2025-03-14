ARG IMG
ARG VERSION


FROM ${IMG}:${VERSION} AS system_stage

RUN apt update
RUN apt install -y apt-utils software-properties-common
RUN apt-add-repository ppa:opm/ppa
RUN apt update
RUN apt install -y openssh-server ca-certificates vim
RUN apt install -y mpi-default-bin libopm-simulators-bin

FROM system_stage AS spe11csp_stage

ENV SPE_CASE=a
ENV INPUT_PATH=/opt/spe11csp/examples/hello_world/
ENV CASEFILE=spe11${SPE_CASE}.txt

RUN apt install -y python3 python3-venv python3-pip git
ARG SPE11_DIR=/opt/spe11csp
RUN git clone https://github.com/OPM/pyopmspe11.git ${SPE11_DIR}

WORKDIR ${SPE11_DIR}
RUN python3 -m venv venv
RUN bash -c "source venv/bin/activate"
# Upgrade pip, setuptools, and wheel
RUN pip install --upgrade pip setuptools wheel
# Install the pyopmcsp11 package (in editable mode for contributions/modifications, i.e., pip install -e .)
RUN pip install .
# For contributions/testing/linting, install the dev-requirements
RUN pip install -r dev-requirements.txt

# to make python cmd execute python3.10
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3 2
RUN update-alternatives --config python

CMD ["sh", "-c", "/usr/local/bin/pyopmspe11 -i ${INPUT_PATH}/${CASEFILE} -o output_csp11${SPE_CASE}"]

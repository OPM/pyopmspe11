[![Build Status](https://github.com/OPM/pyopmspe11/actions/workflows/CI.yml/badge.svg)](https://github.com/OPM/pyopmspe11/actions/workflows/CI.yml)
<a href="https://www.python.org/"><img src="https://img.shields.io/badge/python-3.8%20|%203.9%20|%203.10-blue.svg"></a>
[![Code style](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/ambv/black)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue)](https://opensource.org/license/mit/)

# pyopmspe11: A Python framework using OPM Flow for the SPE11 benchmark project

<img src="docs/text/figs/animation.gif" width="830" height="400">

This repository contains scripts to set up a workflow in python for the three cases in the SPE11 project.
Here we use the [_OPM-Flow_](https://opm-project.org/?page_id=19) simulator.

## Installation
You will first need to install
* Flow (https://opm-project.org)

To build dune and the corresponding OPM master branches from source (e.g., you are a macOS user), you can run the script
`./build_dune_and_opm-flow.bash`, which in turn should build flow in the folder 
./build/opm-simulators/bin/flow (for macOS users the dependecies such as boost can be installed using brew or macports).

You can install the requirements in a virtual environment with the following commands:

```bash
# Clone the repo
git clone https://github.com/OPM/pyopmspe11.git
# Get inside the folder
cd pyopmspe11
# Create virtual environment
python3 -m venv vpyopmspe11
# Activate virtual environment
source vpyopmspe11/bin/activate
# Upgrade pip, setuptools, and wheel
pip install --upgrade pip setuptools wheel
# Install the pyopmspe11 package (in editable mode for contributions/modifications; otherwise, pip install .)
pip install -e .
# For contributions/testing/linting, install the dev-requirements
pip install -r dev-requirements.txt
``` 

See the [_CI.yml_](https://github.com/OPM/pyopmspe11/blob/main/.github/workflows/CI.yml) script
for installation of OPM Flow and the pyopmspe11 package in Linux.

## Running pyopmspe11
You can run _pyopmspe11_ as a single command line:
```
pyopmspe11 -i some_input.txt -o some_output_folder
```
Run `pyopmspe11 --help` to see all possible command line 
argument options. Inside the `some_input.txt` file you provide the path to the
flow executable and simulation parameters. See the .txt files in the examples
folders.

## Getting started
See the [_documentation_](https://OPM.github.io/pyopmspe11/introduction.html).

## About pyopmspe11
The pyopmspe11 package is being funded by the [_HPC Simulation Software for the Gigatonne Storage Challenge project_](https://www.norceresearch.no/en/projects/hpc-simulation-software-for-the-gigatonne-storage-challenge) [project number 622059] and [_Center for Sustainable Subsurface Resources (CSSR)_](https://cssr.no) 
[project no. 331841].
This is work in progress. [_Here_](https://www.spe.org/en/csp/) is the link to the spe11 details.
Contributions are more than welcome using the fork and pull request approach.

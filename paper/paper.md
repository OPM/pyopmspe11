---
title: 'pyopmspe11: A Python framework using OPM Flow for the SPE11 benchmark project'
tags:
  - Python
  - CO$_2$
  - OPM Flow
  - CSP11
  - SPE11
authors:
  - name: David Landa-Marbán
    orcid: 0000-0002-3343-1005
    affiliation: 1
    corresponding: true
  - name: Tor H. Sandve
    orcid: 0000-0002-3267-8276
    affiliation: 1
affiliations:
 - name: NORCE Norwegian Research Centre AS, Bergen, Norway
   index: 1
date: 28 June 2024
bibliography: paper.bib
---

# Summary

The imperative to achieve climate change goals and the increasing worldwide demand for energy 
have made carbon capture and storage (CCS) technology more relevant today.
Since utilizing computational models is essential for planning large-scale CCS projects, then it is crucial
to benchmark simulation tools to enhance confidence in the results. Inspired by
a recent validation study for laboratory-scale CO$_2$ storage [@Flemisch:2024], a new comparative
solution project (CSP) was launched to simulate both lab- and field-scale CO$_2$ storage [@Nordbotten:2024].
As part of the Early Access team, the authors have developed and made open the
`pyopmspe11` tool which facilitates reproducible solutions to the benchmark. This tool serves as a 
common starting point for developing and testing new simulation technology. It is expected that its impact 
will extend far beyond the initial benchmark study.

![Generated model by the configuration file `spe11c_cp_ca20e6cells.txt` in the examples folder.](paper.png){ width=100% }

# Statement of need

Geological carbon storage (GCS) applications benefit from both commercial and open-source simulators. 
However, using open-source simulators often involves data preprocessing and postprocessing, which can be 
challenging for those without advanced programming knowledge. Additionally, setting up and running simulations 
requires engineering expertise. To bridge this gap, developers can simplify the setup of numerical studies by 
using user-friendly approaches, such as configuration files. This not only ensures result reproducibility but 
also facilitates flexible testing of different simulator parameters and allows for easy extension to further studies.

`pyopmspe11` is a simplified and flexible Python tool for the three cases in the SPE Comparative Solution
Project via configuration files. `pyopmspe11` relies on the OPM Flow numerical simulator [@Rassmussen:2021], where the
implementation of the CO$_2$ model can be found in @Sandve:2021. The primary contribution of `pyopmspe11` lies in its 
data preprocessing and postprocessing capabilities. It offers flexibility in generating various types of grids, 
including Cartesian, tensor, and corner-point grids. These grids adhere to standard industry formats, making them 
compatible not only with OPM Flow but also with other simulators. Additionally, `pyopmspe11` supports different resolutions, 
having been tested to generate approximately 160 million cells. In the context of data postprocessing, `pyopmspe11` not only 
generates the necessary reporting data as specified by the benchmark, but it also produces .png figures for rapid inspection 
of individual simulations and for making comparisons between different runs (e.g., to assess sensitivities). The postprocessing 
methods efficiently map non-overlapping cell values (both intensive and extensive quantities) between the simulation grid and 
the reporting grid.

# Outlook
`pyopmspe11` is designed for use by researchers, engineers, and students. During the preliminary intercomparison workshops for 
the benchmark study, the authors received positive feedback about the framework. Some participant groups have utilized `pyopmspe11`, particularly 
for grid generation. Additionally, the authors have been contacted to provide specific support on how to use the tool for setting up 
further studies, such as optimizing injection strategies. Looking ahead, the plan for `pyopmspe11`’s future development includes extending 
its functionality to incorporate additional physical models (e.g., hydrogen storage, salt precipitation, biofilm effects).

# Acknowledgements

The authors acknowledge funding from the [Center for Sustainable Subsurface Resources (CSSR)](https://cssr.no), grant nr. 331841, supported by the Research Council of Norway, research partners NORCE Norwegian Research Centre and the University of Bergen, and user partners Equinor ASA, Wintershall Dea Norge AS, Sumitomo Corporation, Earth Science Analytics, GCE Ocean Technology, and SLB Scandinavia. The authors also acknowledge funding from the [HPC Simulation Software for the Gigatonne Storage Challenge project](https://www.norceresearch.no/en/projects/hpc-simulation-software-for-the-gigatonne-storage-challenge), grant nr. 622059, supported by Equinor ASA and CLIMIT DEMO/Gassnova.

# References
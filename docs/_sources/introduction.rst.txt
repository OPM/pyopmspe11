============
Introduction
============

.. image:: ./figs/animation.gif

This documentation describes the content of the **pyopmcsp11** package.
The numerical simulations for the CO2 are performed using the 
`Flow <https://opm-project.org/?page_id=19>`_ simulator. 

Concept
-------
Simplified and flexible testing framework for the three cases in the SPE Comparative Solution Project
using the open-source simulator OPM Flow. This is work in progress and current efforts focus on the
implementation of dispersion and the thermal effects by the caprock. 

Overview
--------

The current implementation supports the following executable with the argument options:

.. code-block:: bash

    pyopmcsp11 -i input.txt -o output

where 

- \-i, \-input: The base name of the :doc:`configuration file <./configuration_file>` ('input.txt' by default).
- \-o, \-output: The base name of the :doc:`output folder <./output_folder>` ('output' by default).

Installation
------------

See the `Github page <https://github.com/daavid00/pyopmcsp11>`_.

.. tip::
    Check the `CI.yml <https://github.com/daavid00/pyopmcsp11/blob/main/.github/workflows/CI.yml>`_ file.
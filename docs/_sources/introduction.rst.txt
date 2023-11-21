============
Introduction
============

.. image:: ./figs/animation.gif

This documentation describes the content of the **pyopmspe11** package.
The numerical simulations for the CO2 are performed using the 
`Flow <https://opm-project.org/?page_id=19>`_ simulator. 

Concept
-------
Simplified and flexible testing framework for the three cases in the SPE Comparative Solution Project
using the open-source simulator OPM Flow. 

Overview
--------

The current implementation supports the following executable with the argument options:

.. code-block:: bash

    pyopmspe11 -i input.txt -o output

where 

- \-i, \-\-input: The base name of the :doc:`configuration file <./configuration_file>` ('input.txt' by default).
- \-o, \-\-output: The base name of the :doc:`output folder <./output_folder>` ('output' by default).
- \-m, \-\-mode: Run the whole framework ('all'), only create decks ('deck'),  only run flow ('flow'), only write benchmark data ('data'), only create plots ('plot'), deck and run ('deck_flow'), deck, run, and plot (deck_flow_plot), or deck, run, and data (deck_flow_data) ('deck_flow' by default).
- \-g, \-\-generate: Write only the 'dense', 'sparse', 'performance', 'dense_performance', 'performance_sparse', 'dense_sparse', or 'all' ('performance_sparse' by default).
- \-r, \-\-resolution: Number of x, y, and z elements to write the data ('100,10,50' by default).
- \-t, \-\-time: Time interval for the spatial maps (spe11a [h]; spe11b/c [y]) ('5' by default).
- \-c, \-\-compare: Generate a common plot for the current folders for 'spe11a', 'spe11b', or 'spe11c' ('' by default).
    
Installation
------------

See the `Github page <https://github.com/OPM/pyopmspe11>`_.

.. tip::
    Check the `CI.yml <https://github.com/OPM/pyopmspe11/blob/main/.github/workflows/CI.yml>`_ file.
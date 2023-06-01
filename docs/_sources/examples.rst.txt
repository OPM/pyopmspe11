********
Examples
********

======
CSP11C 
======

In this example we consider the configuration file in the
:doc:`configuration file<./configuration_file>` section (also available
in the `examples <https://github.com/daavid00/pyopmcsp11/blob/main/examples>`_ folder).

If the configuration file is saved as 'csp11c.txt' and the generated files are to be
saved in a folder called 'csp11c', then this is achieved by the following command:

.. code-block:: bash

    pyopmcsp11 -i csp11c.txt -o csp11c

The following are screenshots of the simulation results:

.. figure:: figs/facies.png
.. figure:: figs/saturation.png
.. figure:: figs/dissolved.png

    Simulation results of the (top) facies, (middle) saturation, and (bottom) dissolved CO2.

This example uses a very coarser grid to run fast (it takes ca. 5 minutes). See the following
section for finer grids.

============
Presentation 
============

Using the same grid size for the reporting of the results in the CSP11 description, the following computational times
were reported for the 2023 SPE Reservoir Simulation Conference (see the presentation `here <https://github.com/Simulation-Benchmarks/11thSPE-CSP/blob/main/description/SPE11%20CSP.pdf>`_, 
where you can also see some of the preliminary simulation results using OPM Flow):

.. code-block:: yaml

    Case      Dimensions [m]      Max. grid size [m]   No. grid cells  Total no. cells  No. active cells   Solver time step [d]¨  Total simulation time [s]
    csp11a^   [2.8,  0.01,  1.2]  [0.01, 0.01, 0.01]    [280, 1, 120]            33600             31034                   1e-5                    2118.30
    csp11b^*  [8400,    1, 1200]  [10,      1,   10]    [842, 1, 120]           101040             93318                     50                    1420.15
    csp11c^*  [8400, 5000, 1350]  [50,     50,   10]  [170, 100, 120]          2040000           1885200                     50                   25450.68


    ^ All three cases were run with 70 MPI processes and 2 threads per MPI process. i.e., 140 cpu cores.
    * csp11b and csp11c have an extra layer [1 m] of grid cells on the left and right boundaries to include the buffer volume  
    ¨ The solver time step is the maximum value allowed by the simulator

Then, the configuration files in the `examples <https://github.com/daavid00/pyopmcsp11/blob/main/examples>`_ folder can be modified to use the same grid sizes.

.. tip::
    By executing flow --help you get an overview of the available flags in the flow simulator to improve/fix convergence issues 
    (i.e., by setting the flag --solver-max-time-step-in-days to the reported values in the above table).


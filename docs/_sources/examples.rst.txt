********
Examples
********

===========
Hello world 
===========

The `examples/hello_world <https://github.com/OPM/pyopmspe11/blob/main/examples/hello_world>`_ folder contains configuration files
with low grid resolution and shorter injetion times (for initial testing of the framework). For example, by executing:

.. code-block:: bash

    pyopmspe11 -i spe11b.txt -o spe11b -m all -g all -t 5 -r 50,1,20

the following is the figure generated related to the temperature in the domain over time:

.. figure:: figs/spe11b_temp_2Dmaps.png

Let us now tight the mass balance tolerance by adding the flag \-\-tolerance-mb=1e-7 in line 2 of the configuration file (1e-6 is the default one).
Then we run the simulations and we save the results in a different output folder:

.. code-block:: bash

    pyopmspe11 -i spe11b.txt -o stricter_mb -m all -g all -t 5 -r 50,1,20

Then, to visualize the comparison between both runs, this can be achived by executing:

.. code-block:: bash

    pyopmspe11 -c spe11b

The following are the comparison figures generated in the compare folder:

.. figure:: figs/spe11b_sparse_data.png
.. figure:: figs/spe11b_performance.png

This example uses a very coarser grid to run fast. See the following section for finer grids.

============
Presentation 
============

Using the same grid size for the reporting of the results in the SPE11 description, the following computational times
were reported for the 2023 SPE Reservoir Simulation Conference (see the presentation `here <https://github.com/Simulation-Benchmarks/11thSPE-CSP/blob/main/description/SPE11%20CSP.pdf>`_, 
where you can also see some of the preliminary simulation results using OPM Flow):

.. code-block:: yaml

    Case      Dimensions [m]      Max. grid size [m]   No. grid cells  Total no. cells  No. active cells   Solver time step [d]¨  Total simulation time [s]
    spe11a^   [2.8,  0.01,  1.2]  [0.01, 0.01, 0.01]    [280, 1, 120]            33600             31034                   1e-5                    2118.30
    spe11b^*  [8400,    1, 1200]  [10,      1,   10]    [842, 1, 120]           101040             93318                     50                    1420.15
    spe11c^*  [8400, 5000, 1350]  [50,     50,   10]  [170, 100, 120]          2040000           1885200                     50                   25450.68


    ^ All three cases were run with 70 MPI processes and 2 threads per MPI process. i.e., 140 cpu cores.
    * spe11b and spe11c have an extra layer [1 m] of grid cells on the left and right boundaries to include the buffer volume  
    ¨ The solver time step is the maximum value allowed by the simulator

Then `examples/finner_grids <https://github.com/OPM/pyopmspe11/blob/main/examples>`_ folder contains configuration files
with grid size of the same order as requested for reporting the spatial data in the benchmark, as well as the required injection schedules. For example, by executing:

.. code-block:: bash

    pyopmspe11 -i spe11b.txt -o spe11b_finner -m all -g all -r 840,1,120 -t 5

the following are some of the generated figures:

.. figure:: figs/spe11b_sparse_data_finner.png
.. figure:: figs/spe11b_performance_finner.png
.. figure:: figs/spe11b_xco2_2Dmaps_finner.png

.. tip::
    By executing flow --help you get an overview of the available flags in the flow simulator to improve/fix convergence issues 
    (i.e., by setting the flag --linear-solver=cprw to changue the linear solver).
********
Examples
********

===========
Hello world 
===========

The `examples <https://github.com/OPM/pyopmspe11/blob/main/examples>`_ folder contains configuration files
with low grid resolution and shorter injection times (for initial testing of the framework). For example, by executing:

.. code-block:: bash

    pyopmspe11 -i spe11b.txt -o spe11b -m all -g all -t 5 -r 50,1,15 -w 1

The following is the figure `spe11b_tco2_2Dmaps`, which shows the CO2 mass in the domain over time (i.e., the simulations results from
the corner-point grid mapped to the equidistance reporting grid of 50 x 15 as defined by the -r flag). You can
compare your example results to this figure to evaluate if your example ran correctly:

.. figure:: figs/spe11b_tco2_2Dmaps.png

Using the :ref:`toml` format, the previous run is equivalent to (this requires a Python version of at least 3.11 [due to `tomllib <https://toml.io/en/>`_]):

.. code-block:: bash

    pyopmspe11 -i spe11b.toml -o spe11b -m all -g all -t 5 -r 50,1,15 -w 1

Let us now change the grid type from corner-point to tensor in line 7 of the configuration file.
Then, we run the simulations and we save the results in a different output folder:

.. code-block:: bash

    pyopmspe11 -i spe11b.txt -o spe11b_tensor_grid -m deck_flow_data -g performance_sparse -t 5 -r 50,1,15 -w 1

Here we have just set the framework to generate the deck, run the simulations, and generate the performance and sparse data.
Then, to visualize the comparison between both runs, this can be achived by executing:

.. code-block:: bash

    pyopmspe11 -c spe11b

The following are some of the figures generated in the compare folder:

.. figure:: figs/spe11b_performance.png
.. figure:: figs/spe11b_sparse_data.png

This example uses a very coarser grid to run fast. See the :doc:`benchmark <./benchmark>` section for finer grids. 

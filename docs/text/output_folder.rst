=============
Output folder
=============

The following screenshot shows an example of generated files in the selected output folder after 
executing **pyopmspe11**.

.. figure:: figs/output.png

    Generated files in the :doc:`hello world <./examples>` example.

The simulation results are saved in the **flow** folder, and
`ResInsight <https://resinsight.org>`_ and `plopm <https://github.com/cssr-tools/plopm>`_ can be used for the visualization.
In addition, some figures are plotted in png format in the **figures** folder.
Then after running **pyopmspe11**, one could modify the generated OPM related files and 
run directly the simulations calling the Flow solvers, e.g., to add tracers, salinity, rock compressibility, etc. 
(see the OPM Flow documentation `here <https://opm-project.org/?page_id=955>`_).

.. tip::
    If you install the dev-requirements.txt, executing **pytest \-\-cov=pyopmspe11 \-\-cov-report term-missing tests/** runs the
    configuration files in the tests with different argument options for **pyopmspe11** in ca. 30 minutes, then one can see the different outputs.

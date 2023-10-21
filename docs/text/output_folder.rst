=============
Output folder
=============

The following screenshot shows the generated files in the selected output folder after 
executing **pyopmspe11**.

.. figure:: figs/output.png

    Generated files after executing **pyopmspe11**.

The simulation results are saved in the **flow** folder, and
`ResInsight <https://resinsight.org>`_ can be used for the visualization.
In addition, some figures are plotted in png format in the **figures** folder.
Then after running **pyopmspe11**, one could modify the generated OPM related files and 
run directly the simulations calling the Flow solvers, e.g., to add tracers 
(see the OPM Flow documentation `here <https://opm-project.org/?page_id=955>`_). 
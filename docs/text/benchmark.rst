*********
Benchmark
*********

.. note::
    These are preliminary results and will be updated up to the final submission deadline (September 20th 2024, see `here <https://connect.spe.org/discussion/update-on-spe11-submission-and-workshops#bm94ca322f-a2f3-421b-9f90-0191376f7b25>`_).
    For example, new simulations for the spe11a r4 and spe11c r4 are still running to handle the mass issue and to include the 1000 years of initialization time, respectively. 

All configuration files are located in the `benchmark <https://github.com/OPM/pyopmspe11/blob/main/benchmark>`_ folder.

======
SPE11A
======

* r1_Cart_1cm                
    Uniform grid of 1 cm size
* r2_Cart_1cm_capmax2500Pa   
    Uniform grid of 1 cm size (maximum capillary pressure of 2500 Pa instead of 95000 Pa following the remarks in the benchmark description)
* r3_cp_1cmish_capmax2500Pa  
    Corner-point grid of ca. 1 cm size (maximum capillary pressure of 2500 Pa instead of 95000 Pa following the remarks in the benchmark description)
* r4_Cart_1mm_capmax2500Pa   
    Uniform grid of 1 mm size (maximum capillary pressure of 2500 Pa instead of 95000 Pa following the remarks in the benchmark description)
 
To run the cases in the terminal:

.. code-block:: bash

    pyopmspe11 -i r1_Cart_1cm.txt -o r1_Cart_1cm -m all -g all -t 1 -r 280,1,120 -w 0.16666666666666666
    pyopmspe11 -i r2_Cart_1cm_capmax2500Pa.txt -o r2_Cart_1cm_capmax2500Pa -m all -g all -t 1 -r 280,1,120 -w 0.16666666666666666
    pyopmspe11 -i r3_cp_1cmish_capmax2500Pa.txt -o r3_cp_1cmish_capmax2500Pa -m all -g all -t 1 -r 280,1,120 -w 0.16666666666666666
    pyopmspe11 -i r4_Cart_1mm_capmax2500Pa.txt -o r4_Cart_1mm_capmax2500Pa -m all -g all -t 1 -r 280,1,120 -w 0.16666666666666666

As mentioned in the `CSP description <https://onepetro.org/SJ/article/29/05/2507/540636/The-11th-Society-of-Petroleum-Engineers>`_, using the maximum value of 2500 Pa instead of
95000 Pa does not significantly impact the results (comparing r1 and r2), and this choice also reduces the simulation time. In addition, the corner-point grid results (r3) compare 
very well to the fine-scale simulations (r4).

----------------
Performance data
----------------

.. figure:: figs/benchmark_spe11a_performance.png

-----------
Sparse data
-----------

.. figure:: figs/benchmark_spe11a_sparse_data.png

------------
Spatial maps
------------

.. figure:: figs/massfracta.png

======
SPE11B
======

* r1_Cart_10m                
    Uniform grid of 10 m size (1 m dx size on left and right boundaries)
* r2_cp_10mish   
    Corner-point grid of ca. 10 m size (1 m dx size on left and right boundaries)
* r3_cp_10mish_convective 
    Corner-point grid of ca. 10 m size (1 m dx size on left and right boundaries) using a subgrid model for convective mixing for facies 2 and 5.
* r4_Cart_1m    
    Uniform grid of 1 m size
 
To run the cases in the terminal:

.. code-block:: bash

    pyopmspe11 -i r1_Cart_10m.txt -o r1_Cart_10m -m all -g all -r 840,1,120 -t 5 -w 0.1
    pyopmspe11 -i r2_cp_10mish.txt -o r2_cp_10mish -m all -g all -r 840,1,120 -t 5 -w 0.1
    pyopmspe11 -i r3_cp_10mish_convective.txt -o r3_cp_10mish_convective -m all -g all -r 840,1,120 -t 5 -w 0.1
    pyopmspe11 -i r4_Cart_1m.txt -o r4_Cart_1m -m all -g all -r 840,1,120 -t 5 -w 0.1

For the box A quantities, the convective model results (r3) compare very well to the fine-scale simulations (r4), which runs ca. 500 times faster. Details on the convective model
will be available in Mykkeltvedt et al., under review. For the implementation in OPM Flow of this model, it is work in progress to handle better the zones in the geological model where 
dissolve CO2 accumulates.

----------------
Performance data
----------------

.. figure:: figs/benchmark_spe11b_performance.png

-----------
Sparse data
-----------

.. figure:: figs/benchmark_spe11b_sparse_data.png

------------
Spatial maps
------------

.. figure:: figs/massfractb.png

======
SPE11C
======

* r1_Cart_50m-50m-10m                
    Grid of [50, 50, 10] m size (1 m dx size on left and right boundaries and 1 m dy size on back and front boundaries)
* r2_cp_50m-50m-8mish   
    Corner-point grid of [50, 50, mean ca. 8] m size (1 m dx size on left and right boundaries and 1 m dy size on back and front boundaries)
* r3_cp_50m-50m-8mish_convective 
    Corner-point grid of [50, 50, mean ca. 8] m size (1 m dx size on left and right boundaries and 1 m dy size on back and front boundaries) using the convective model for facies 2 and 5.
* r4_cp_8m-8mish-8mish    
    Corner-point grid of [8, mean ca. 8, mean ca. 8] m size (1 m dx size on left and right boundaries and 1 m dy size on back and front boundaries)
 
To run the cases in the terminal:

.. code-block:: bash

    pyopmspe11 -i r1_Cart_50m-50m-10m.txt -o r1_Cart_50m-50m-10m -m all -g all -r 168,100,120 -t 0,5,10,15,20,25,30,35,40,45,50,75,100,150,200,250,300,350,400,450,500,600,700,800,900,1000 -w 0.1
    pyopmspe11 -i r2_cp_50m-50m-8mish.txt -o r2_cp_50m-50m-8mish -m all -g all -r 168,100,120 -t 0,5,10,15,20,25,30,35,40,45,50,75,100,150,200,250,300,350,400,450,500,600,700,800,900,1000 -w 0.1
    pyopmspe11 -i r3_cp_50m-50m-8mish_convective.txt -o r3_cp_50m-50m-8mish_convective -m all -g all -r 168,100,120 -t 0,5,10,15,20,25,30,35,40,45,50,75,100,150,200,250,300,350,400,450,500,600,700,800,900,1000 -w 0.1
    pyopmspe11 -i r4_cp_8m-8mish-8mish.txt -o r4_cp_8m-8mish-8mish -m all -g all -r 168,100,120 -t 0,5,10,15,20,25,30,35,40,45,50,75,100,150,200,250,300,350,400,450,500,600,700,800,900,1000 -w 0.1 -u opm

To run the case with more than 100 million cells (r4), it required improvements in the OPM Flow simulator (e.g., `METIS support <https://github.com/OPM/opm-grid/pull/738>`_), as well as in the **pyopmspe11** pre- and postprocessing functionality
(and of course a big computer, the `NORCE <https://www.norceresearch.no>`_ HPC cluster). See `this gif <https://github.com/OPM/pyopmspe11/blob/main/docs/text/figs/pyopmspe11c100Mcells.gif>`_ for visualization in `ParaView <https://www.paraview.org>`_ of the CO2 molar fraction (liquid phase) over time.

----------------
Performance data
----------------

.. figure:: figs/benchmark_spe11c_performance.png

-----------
Sparse data
-----------

.. figure:: figs/benchmark_spe11c_sparse_data.png

------------
Spatial maps
------------

.. figure:: figs/massfractc.png
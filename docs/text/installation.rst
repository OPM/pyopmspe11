============
Installation
============

The following steps work installing the dependencies in Linux via apt-get or in macOS using brew or macports.
While using packages managers such as Anaconda, Miniforge, or Mamba might work, these are not tested.
The supported Python versions are 3.11 to 3.14.

.. _vpyopmspe11:

Python package
--------------

To install the **pyopmspe11** executable from the development version: 

.. code-block:: bash

    pip install git+https://github.com/opm/pyopmspe11.git

If you are interested in a specific version (e.g., v2025.10) or in modifying the source code, then you can clone the repository and 
install the Python requirements in a virtual environment with the following commands:

.. code-block:: console

    # Clone the repo
    git clone https://github.com/opm/pyopmspe11.git
    # Get inside the folder
    cd pyopmspe11
    # For a specific version (e.g., v2025.10), or skip this step (i.e., edge version)
    git checkout v2025.10
    # Create virtual environment (to specific Python, python3.13 -m venv vpyopmspe11)
    python3 -m venv vpyopmspe11
    # Activate virtual environment
    source vpyopmspe11/bin/activate
    # Upgrade pip, setuptools, and wheel
    pip install --upgrade pip setuptools wheel
    # Install the pyopmspe11 package
    pip install -e .
    # For contributions/testing/linting, install the dev-requirements
    pip install -r dev-requirements.txt

.. tip::

    Typing **git tag -l** writes all available specific versions.

.. note::
  
    For not macOS users, to install the (optional but recommended) dependencies used for the figure's LaTeX formatting, execute 
    
    **sudo apt-get install texlive-fonts-recommended texlive-fonts-extra dvipng cm-super**

    For macOS users, see :ref:`macOS`.

OPM Flow
--------
You also need to install:

* OPM Flow (https://opm-project.org, Release 2025.10 or current master branches)

.. tip::

    See the `CI.yml <https://github.com/opm/pyopmspe11/blob/main/.github/workflows/CI.yml>`_ script 
    for installation of OPM Flow (binary packages) and the pyopmspe11 package in Ubuntu. 

Source build in Linux/Windows
+++++++++++++++++++++++++++++
If you are a Linux user (including the Windows subsystem for Linux, see `this link <https://learn.microsoft.com/en-us/windows/python/web-frameworks>`_ 
for a nice tutorial for setting Python environments in WSL), then you could try to build flow (after installing the `prerequisites <https://opm-project.org/?page_id=239>`_) from the master branches with mpi support by running
in the terminal the following lines (which in turn should build flow in the folder ./build/opm-simulators/bin/flow): 

.. code-block:: console

    CURRENT_DIRECTORY="$PWD"

    mkdir build

    for repo in common grid
    do  git clone https://github.com/OPM/opm-$repo.git
        mkdir build/opm-$repo
        cd build/opm-$repo
        cmake -DUSE_MPI=1 -DWITH_NDEBUG=1 -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="$CURRENT_DIRECTORY/build/opm-common;$CURRENT_DIRECTORY/build/opm-grid" $CURRENT_DIRECTORY/opm-$repo
        if [[ $repo == simulators ]]; then
            make -j5 flow
        else
            make -j5 opm$repo
        fi
        cd ../..
    done


.. tip::

    You can create a .sh file (e.g., build_opm_mpi.sh), copy the previous lines, and run in the terminal **source build_opm_mpi.sh**

.. _macOS:

Source build in macOS
+++++++++++++++++++++
For macOS, there are no available binary packages, so OPM Flow needs to be built from source, in addition to the dune libraries and the OPM Python
package (see the `prerequisites <https://opm-project.org/?page_id=239>`_, which can be installed using macports or brew). 
For example, with brew the prerequisites can be installed by:

.. code-block:: console

    brew install boost openblas suite-sparse python@3.14 cmake

In addition, it is recommended to uprade and update your macOS to the latest available versions (the following steps have 
worked for macOS Tahoe 26.0.1 with Apple clang version 17.0.0).
After the prerequisites are installed and the vpyopmspe11 Python environment is created (see :ref:`vpyopmspe11`), 
then building OPM Flow can be achieved with the following bash lines:

.. code-block:: console

    CURRENT_DIRECTORY="$PWD"

    deactivate
    source vpyopmspe11/bin/activate

    for module in common geometry grid istl
    do   git clone https://gitlab.dune-project.org/core/dune-$module.git --branch v2.9.1
        ./dune-common/bin/dunecontrol --only=dune-$module cmake -DCMAKE_DISABLE_FIND_PACKAGE_MPI=1
        ./dune-common/bin/dunecontrol --only=dune-$module make -j5
    done

    mkdir build

    for repo in common grid simulators
    do  git clone https://github.com/OPM/opm-$repo.git
        mkdir build/opm-$repo
        cd build/opm-$repo
        cmake -DUSE_MPI=0 -DWITH_NDEBUG=1 -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="$CURRENT_DIRECTORY/dune-common/build-cmake;$CURRENT_DIRECTORY/dune-grid/build-cmake;$CURRENT_DIRECTORY/dune-geometry/build-cmake;$CURRENT_DIRECTORY/dune-istl/build-cmake;$CURRENT_DIRECTORY/build/opm-common;$CURRENT_DIRECTORY/build/opm-grid" $CURRENT_DIRECTORY/opm-$repo
        if [[ $repo == simulators ]]; then
            make -j5 flow
        else
            make -j5 opm$repo
        fi
        cd ../..
    done

    echo "export PATH=\$PATH:$CURRENT_DIRECTORY/build/opm-simulators/bin" >> $CURRENT_DIRECTORY/vpyopmspe11/bin/activate

    deactivate
    source vpyopmspe11/bin/activate

This builds OPM Flow, and it exports the path to the flow executable (i.e., executing in the terminal **which flow** should print the path).

For macOS, the LaTeX dependency can be installed from https://www.tug.org/mactex/.

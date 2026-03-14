============
Installation
============

The following steps work installing the dependencies in Linux via apt-get or in macOS using brew or macports.
While using packages managers such as Anaconda, Miniforge, or Mamba might work, these are not tested.
The supported Python versions are 3.11 to 3.14.

`ResInsight <https://resinsight.org>`_ and `plopm <https://github.com/cssr-tools/plopm>`_ are used for the visualization of the results.

.. note::

    There are binary packages for Linux and Windows to install Resinsight, see the `ResInsight Documentation <https://resinsight.org/releases/>`_. For macOS users, you could try to install it using `brew <https://brew.sh>`_ by executing:

    .. code-block:: bash

        brew install cssr-tools/opm/resinsight
    
    Then, you should be able to open resinsight by typing in the terminal **resinsight**.

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

    For macOS users, the LaTeX dependency can be installed from https://www.tug.org/mactex/.

OPM Flow
--------
You also need to install:

* OPM Flow (https://opm-project.org, Release 2025.10 or current master branches)

Binary packages
+++++++++++++++

See the `downloading and installing <https://opm-project.org/?page_id=36>`_ OPM Flow online documentation for 
instructions to install the binary packages in Ubuntu and Red Hat Enterprise Linux, and for other platforms which are
supported either via source builds or through running a virtual machine.

.. tip::

    See the `CI.yml <https://github.com/opm/pyopmspe11/blob/main/.github/workflows/CI.yml>`_ script 
    for installation of OPM Flow (binary packages), LaTeX (optional) libraries, and the pyopmspe11 package in Ubuntu.

Source build in Linux/Windows
+++++++++++++++++++++++++++++
If you are a Linux user (including the Windows subsystem for Linux, see `this link <https://learn.microsoft.com/en-us/windows/python/web-frameworks>`_ 
for a nice tutorial for setting Python environments in WSL), then you could try to build flow (after installing the `prerequisites <https://opm-project.org/?page_id=239>`_) from the master branches with mpi support by running
in the terminal the following lines (which in turn should build flow in the folder ./build/opm-simulators/bin/flow): 

.. code-block:: console

    CURRENT_DIRECTORY="$PWD"

    mkdir build

    for repo in common grid simulators
    do  git clone https://github.com/OPM/opm-$repo.git
        mkdir build/opm-$repo
        cd build/opm-$repo
        cmake -DWITH_NDEBUG=1 -DCMAKE_BUILD_TYPE=Release $CURRENT_DIRECTORY/opm-$repo
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

Brew formula for macOS
++++++++++++++++++++++
For macOS, there are no available binary packages, so OPM Flow needs to be built from source. Recently, a formula to build flow using brew has
been added in `https://github.com/cssr-tools/homebrew-opm <https://github.com/cssr-tools/homebrew-opm>`_. 
Then, you can try to install flow (v2026.02) by simply typing:

.. code-block:: console

    brew install cssr-tools/opm/opm-simulators

You can check if the installation of OPM Flow succeded by typing in the terminal **flow \-\-help**.

.. tip::
    See the actions in the `cssr-tools/homebrew-opm <https://github.com/cssr-tools/homebrew-opm/actions>`_ repository. 

Source build in macOS
+++++++++++++++++++++
See `this repository <https://github.com/daavid00/OPM-Flow_macOS>`_ dedicated to build OPM Flow from source in the latest macOS (GitHub actions), and tested with **pycopm**, a repository in cssr-tools.

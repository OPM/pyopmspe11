************
Contributing
************

Contributions are more than welcome using the fork and pull request approach 🙂 (if you are not familiar with this approach, 
please visit `GitHub Docs PRs <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests>`_ for an extended 
documentation about collaborating with pull request; also, looking at previous merged pull requests helps to get familiar with this).

============
Ground Rules 
============

- We use Black code formatting
- We use Pylint
- We document our code

==========================
Contribute to the software
==========================

This requires a Python version of at least 3.11 (due to `tomllib <https://toml.io/en/>`_).

#. Work on your own fork of the main repo
#. In the main repo execute:

    #. **pip install -r dev-requirements.txt** (this installs the `dev-requirements.txt <https://github.com/cssr-tools/pyopmspe11/blob/main/dev-requirements.txt>`_; in addition, both opm Python and LaTeX are required, then for not macOs users run **pip install opm** and **sudo apt-get install texlive-fonts-recommended texlive-fonts-extra dvipng cm-super**, or else follow the instructions in `macOS installation <https://cssr-tools.github.io/pyopmspe11/installation.html#source-build-in-macos>`_)
    #. **black src/ tests/** (this formats the code)
    #. **pylint src/ tests/** (this analyses the code, and might rise issues that need to be fixed before the pull request)
    #. **mypy --ignore-missing-imports src/ tests/** (this is a static checker, and might rise issues that need to be fixed before the pull request)
    #. **pytest --cov=pyopmspe11 --cov-report term-missing tests/** (this runs locally the tests, and might rise issues that need to be fixed before the pull request)
    #. **pushd docs & make html** (this generates the documentation, and might rise issues that need to be fixed before the pull request; if the build succeeds and if the contribution changes the documentation, then delete all content from the `docs <https://github.com/cssr-tools/pyopmspe11/tree/main/docs>`_ folder except `Makefile <https://github.com/OPM/pyopmspe11/blob/main/docs/Makefile>`_, `text <https://github.com/OPM/pyopmspe11/blob/main/docs/text>`_, and `.nojekyll <https://github.com/OPM/pyopmspe11/blob/main/docs/.nojekyll>`_, after copy all contents from the docs/_build/html/ folder, and finally paste them in the `docs <https://github.com/cssr-tools/pyopmspe11/tree/main/docs>`_ folder)
    
    .. tip::
        See the `CI.yml <https://github.com/cssr-tools/pyopmspe11/blob/main/.github/workflows/CI.yml>`_ script and the `Actions <https://github.com/cssr-tools/pyopmspe11/actions>`_ for installation of pyopmspe11, OPM Flow (binary packages), and dependencies, as well as the execution of the six previous steps in Ubuntu 24.10 using Python3.11.

#. Squash your commits into a single commit (see this `nice tutorial <https://gist.github.com/lpranam/4ae996b0a4bc37448dc80356efbca7fa>`_ if you are not familiar with this)
#. Push your commit and make a pull request
#. The maintainers will review the pull request, and if the contribution is accepted, then it will be merge to the main repo

============================
Reporting issues or problems
============================

#.  Issues or problems can be raised by creating a `new issue <https://github.com/cssr-tools/pyopmspe11/issues>`_ in the repository GitHub page (if you are not familiar with this approach, please visit `GitHub Docs Issues <https://docs.github.com/en/issues/tracking-your-work-with-issues>`_).
#.  We will try to answer as soon as possible, but also any user is more than welcome to answer.

============
Seek support
============

#.  The preferred approach to seek support is to raise an Issue as described in the previous lines.
#.  We will try to answer as soon as possible, but also any user is more than welcome to answer.

- An alternative approach is to send an email to any of the `mantainers <https://github.com/cssr-tools/pyopmspe11/blob/main/pyproject.toml>`_.
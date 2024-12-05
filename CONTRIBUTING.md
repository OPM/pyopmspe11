# Contributing

Contributions are more than welcome using the fork and pull request approach 🙂 (if you are not familiar with this approach, please visit [_GitHub Docs PRs_](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests) for an extended documentation about collaborating with pull request; also, looking at previous merged pull requests helps to get familiar with this).

## Ground Rules

- We use Black code formatting
- We use Pylint
- We document our code

## Contribute to the software
This requires a Python version of at least 3.11 (due to [_tomllib_](https://toml.io/en/)).

1. Work on your own fork of the main repo
1. In the main repo execute:
    1. **pip install -r dev-requirements.txt** (this installs the [_dev-requirements.txt_](https://github.com/cssr-tools/pyopmspe11/blob/main/dev-requirements.txt); in addition, both opm Python and LaTeX are required, then for not macOs users run **pip install opm** and **sudo apt-get install texlive-fonts-recommended texlive-fonts-extra dvipng cm-super**, or else follow the instructions in [_macOS installation_](https://cssr-tools.github.io/pyopmspe11/installation.html#source-build-in-macos))
    1. **black src/ tests/** (this formats the code)
    1. **pylint src/ tests/** (this analyses the code, and might rise issues that need to be fixed before the pull request)
    1. **mypy --ignore-missing-imports src/ tests/** (this is a static checker, and might rise issues that need to be fixed before the pull request)
    1. **pytest --cov=pyopmspe11 --cov-report term-missing tests/** (this runs locally the tests, and might rise issues that need to be fixed before the pull request)
    1. **pushd docs & make html** (this generates the documentation, and might rise issues that need to be fixed before the pull request; if the build succeeds and if the contribution changes the documentation, then delete all content from the [_docs_](https://github.com/cssr-tools/pyopmspe11/tree/main/docs) folder except [_Makefile_](https://github.com/OPM/pyopmspe11/blob/main/docs/Makefile), [_text_](https://github.com/OPM/pyopmspe11/blob/main/docs/text), and [_.nojekyll_](https://github.com/OPM/pyopmspe11/blob/main/docs/.nojekyll), after copy all contents from the docs/_build/html/ folder, and finally paste them in the [_docs_](https://github.com/cssr-tools/pyopmspe11/tree/main/docs) folder)
    * Tip: See the [_CI.yml_](https://github.com/cssr-tools/pyopmspe11/blob/main/.github/workflows/CI.yml) script and the [_Actions_](https://github.com/cssr-tools/pyopmspe11/actions) for installation of pyopmspe11, OPM Flow (binary packages), and dependencies, as well as the execution of the six previous steps in Ubuntu 24.10 using Python3.11.
1. Squash your commits into a single commit (see this [_nice tutorial_](https://gist.github.com/lpranam/4ae996b0a4bc37448dc80356efbca7fa) if you are not familiar with this)
1. Push your commit and make a pull request
1. The maintainers will review the pull request, and if the contribution is accepted, then it will be merge to the main repo 

## Reporting issues or problems

1. Issues or problems can be raised by creating a [_new issue_](https://github.com/cssr-tools/pyopmspe11/issues) in the repository GitHub page (if you are not familiar with this approach, please visit [_GitHub Docs Issues_](https://docs.github.com/en/issues/tracking-your-work-with-issues)).
1. We will try to answer as soon as possible, but also any user is more than welcome to answer. 

## Seek support

1. The preferred approach to seek support is to raise an issue as described in the previous lines.
1. We will try to answer as soon as possible, but also any user is more than welcome to answer.
- An alternative approach is to send an email to any of the [_mantainers_](https://github.com/cssr-tools/pyopmspe11/blob/main/pyproject.toml). 
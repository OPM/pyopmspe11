# SPDX-FileCopyrightText: 2023 NORCE
# SPDX-License-Identifier: MIT

"""pyopmspe11: A Python framework using OPM Flow for the SPE11 benchmark project"""

from setuptools import find_packages, setup

with open("README.md", "r", encoding="utf8") as file:
    long_description = file.read()

with open("requirements.txt", "r", encoding="utf8") as file:
    install_requires = file.read().splitlines()

with open("dev-requirements.txt", "r", encoding="utf8") as file:
    dev_requires = file.read().splitlines()

setup(
    name="pyopmspe11",
    version="2023.10-pre",
    install_requires=install_requires,
    extras_require={"dev": dev_requires},
    setup_requires=["setuptools_scm"],
    description="pyopmspe11: A Python framework using OPM Flow for the SPE11 benchmark project",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/OPM/pyopmspe11",
    author="David Landa-Marbán",
    mantainer="David Landa-Marbán",
    mantainer_email="dmar@norceresearch.no",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Topic :: Scientific/Engineering",
        "License :: OSI Approved :: MIT",
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    keywords="csp11 spe11 co2 flow python",
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    license="MIT",
    python_requires=">=3.8, <4",
    entry_points={
        "console_scripts": [
            "pyopmspe11=pyopmspe11.core.pyopmspe11:main",
        ]
    },
    include_package_data=True,
)

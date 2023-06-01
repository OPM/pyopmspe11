# SPDX-FileCopyrightText: 2023 NORCE
# SPDX-License-Identifier: MIT

"""pyopmcsp11: A Python framework using OPM Flow for the CSP11 benchmark project"""

from setuptools import find_packages, setup

with open("README.md", "r", encoding="utf8") as file:
    long_description = file.read()

with open("requirements.txt", "r", encoding="utf8") as file:
    install_requires = file.read().splitlines()

with open("dev-requirements.txt", "r", encoding="utf8") as file:
    dev_requires = file.read().splitlines()

setup(
    name="pyopmcsp11",
    version="0.0.1",
    install_requires=install_requires,
    extras_require={"dev": dev_requires},
    setup_requires=["setuptools_scm"],
    description="pyopmcsp11: A Python framework using OPM Flow for the CSP11 benchmark project",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/dmar/pyopmcsp11",
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
    keywords="spe11 co2 flow python",
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    license="MIT",
    python_requires=">=3.8, <4",
    entry_points={
        "console_scripts": [
            "pyopmcsp11=pyopmcsp11.core.pyopmcsp11:main",
        ]
    },
    include_package_data=True,
)

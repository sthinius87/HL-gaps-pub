#!/usr/bin/env python

"""The setup script."""
import re
from pathlib import Path

from setuptools import find_packages, setup

VERSION = "0.1.0"

module_dir = Path(__file__).resolve().parent
with open(module_dir / "hl_gaps_pub" / "__init__.py", encoding="utf-8") as version_file:
    for line in version_file:
        match = re.match(r'__version__ = "(.*)"', line)
        if match is not None:
            _version = match.group(1)
            break
    else:
        raise RuntimeError(
            f"Could not determine package version from {version_file.name} !"
        )
    if _version != VERSION:
        raise RuntimeError(
            f"Version mismatch, check {version_file.name} and {module_dir/'setup.py'} !"
        )

with open("README.rst", encoding="utf-8") as readme_file:
    readme = readme_file.read()

with open("HISTORY.rst", encoding="utf-8") as history_file:
    history = history_file.read()

 # https://github.com/pypa/setuptools/issues/2769
install_requirements = [
    "setuptools", "click", "environs", ]
setup_requirements = [ "pytest-runner", ]
test_requirements = [ "pytest>=3", ]

setup(
    author="Sascha Thinius",
    author_email="sascha.thinius@ifam.fraunhofer.de",
    python_requires=">=3.9",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",  
    ],
    description="High throughput tight binding calculation of electronic HOMO-LUMO gaps and its prediction for natural compounds",
    entry_points={
        "console_scripts": [
            "HLgap=hl_gaps_pub.cli:get_hl_gap",
        ],
    },   
    install_requires=install_requirements,
    license="MIT license",
    long_description=readme + "\n\n" + history,
    include_package_data=True,
    package_data={
        "hl_gaps_pub": [
            "py.typed",
        ]
    },
    keywords="HL-gaps-pub",
    name="HL-gaps-pub",
    packages=find_packages(include=["hl_gaps_pub", "hl_gaps_pub.*"]),
    setup_requires=setup_requirements,
    test_suite="tests",
    tests_require=test_requirements,
    url="https://HL_gaps_pub/ifam418/HL-gaps-pub",
    version=VERSION,
    zip_safe=False,
)


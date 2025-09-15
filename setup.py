#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages
from pathlib import Path

RELEASE = "1.0.0"

readme = Path("README.rst").read_text(encoding="utf-8") if Path("README.rst").exists() else ""
history = Path("HISTORY.rst").read_text(encoding="utf-8") if Path("HISTORY.rst").exists() else ""

# Runtime dependencies
requirements = [
    "MDAnalysis>=2.8.0",
    "dask[dataframe]>=2024.1.0",
    "numpy>=1.23",
    "pandas>=1.5",
]

# Optional extras (install with: pip install .[tests])
extras_require = {
    "tests": [
        "pytest>=7",
        "pytest-xdist>=3",
        "numpy",  # keep for assert helpers
        # If you publish ENPMDATests to PyPI and want to pull a pinned version, add it here.
        # Otherwise tests in this repo use the local test package under tests/ENPMDATests.
        # "ENPMDATests=={}".format(RELEASE),
    ],
}

setup(
    name="ENPMDA",
    version=RELEASE,
    description="Parallel analysis for ensemble simulations",
    long_description=readme + "\n\n" + history,
    long_description_content_type="text/x-rst",
    author="Yuxuan Zhuang",
    author_email="wsygzyx@gmail.com",
    url="https://github.com/yuxuanzhuang/ENPMDA",
    license="GNU General Public License v3",
    packages=find_packages(include=["ENPMDA", "ENPMDA.*"]),
    include_package_data=True,
    zip_safe=False,
    python_requires=">=3.10",
    install_requires=requirements,
    extras_require=extras_require,
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Topic :: Scientific/Engineering",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
    keywords="ENPMDA MDAnalysis Dask molecular-dynamics",
)
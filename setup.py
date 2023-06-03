# -*- encoding: utf-8 -*-
"""
SageMath package for computing various moduli space invariants.
"""

from setuptools import setup

with open("README.rst") as f:
    long_description = f.read()

version = {}
with open("msinvar/version.py") as f:
    exec(f.read(), version)

setup(
    name="msinvar",
    author="Sergey Mozgovoy",
    author_email="mozhov@gmail.com",
    description="SageMath package for computing various moduli space invariants",
    long_description=long_description,
    version=version["version"],
    url="https://github.com/smzg/msinvar",
    license="GNU General Public License, version 2",
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Mathematics",
        "License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)",
        "Programming Language :: Python",
    ],
    keywords="SageMath, moduli spaces, invariants",
    packages=["msinvar"],
    install_requires=[],
)

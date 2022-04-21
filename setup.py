#!/usr/bin/env python

from setuptools import find_packages, setup

exec(open("src/theiform/_version.py").read())

setup(
    name="theiform",
    version=__version__,  # type: ignore
    author="Soichiro Hattori",
    author_email="soichiro.hattori@gmail.com",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    license="MIT",
    description="A lightweight Python package to download and interact with ZTF light curve data",
    long_description=open("README.rst").read(),
    install_requires=["numpy>=1.20.2", "astropy>=4.2.1", "requests>=2.25.1"],
)

import os

from setuptools import find_packages, setup

setup(
    name="april",
    version="0.1.0",
    author="Soichiro Hattori",
    author_email="soichiro.hattori@gmail.com",
    packages=find_packages(where="src"),
    license="MIT",
    description="A lightweight Python package to download and interact with ZTF light curve data",
    long_description=open("README.rst").read(),
    install_requires=[
        "numpy>=1.20.2",
        "astropy>=4.2.1",
        "requests>=2.25.1"
    ],
)
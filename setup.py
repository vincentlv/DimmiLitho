import io
import os
import re

from setuptools import find_packages
from setuptools import setup


def get_install_requires():
    with open("requirements.txt", "r") as f:
        return [line.strip() for line in f.readlines() if not line.startswith("-")]


setup(
    name="litho",
    version="0.0.2",
    url="https://github.com/joamatab/dimmilitho",
    license="MIT",
    author="vincentlv",
    description="litho simulation",
    packages=find_packages(exclude=("tests",)),
    install_requires=get_install_requires(),
    python_requires=">=3.6",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.7",
    ],
)

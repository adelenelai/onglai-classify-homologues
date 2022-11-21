#!/usr/bin/env python

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="onglai-classify-homologues",
    version="1.0.0",
    author="Adelene Lai",
    author_email="adelene.lai@uni.lu",
    maintainer="Adelene Lai",
    maintainer_email="adelene.lai@uni.lu",
    description="A cheminformatics algorithm to classify homologous series.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/adelenelai/onglai-classify-homologues",
    packages=setuptools.find_packages(),
    license="Apache",
    install_requires=[
        "pandas",
        "numpy",
        "matplotlib",
        "pytest",
        "rdkit",
        "datamol",
    ],
    package_data={"onglai-classify-homologues":["nextgen_classify_homologues*.*", "utils*.*"]},
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3 :: Only",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.5",
)

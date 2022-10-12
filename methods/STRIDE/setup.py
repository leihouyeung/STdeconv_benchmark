# -*- coding: utf-8 -*-
# @Author: Dongqing Sun
# @E-mail: Dongqingsun96@gmail.com
# @Date:   2021-06-10 15:25:08
# @Last Modified by:   Dongqing Sun
# @Last Modified time: 2021-08-03 21:36:20


import sys,os

try:
    from setuptools import setup, find_packages
except ImportError:
    print("Could not load setuptools. Please install the setuptools package.")

exec(open('src/STRIDE/version.py').read())

with open("requirements.txt") as f:
    requirements = f.read().splitlines()

def main():
    setup(
        name = "stridespatial",
        package_dir = {'':'src'},
        version = __version__,
        packages = find_packages(where="src"),
        scripts = ['bin/STRIDE'],
        package_data={
            "":["*.txt"]
        },
        install_requires = requirements,
        setup_requires = requirements,
        include_package_data = True,
        author = "Dongqing Sun",
        author_email = "Dongqingsun96@gmail.com",
        description = "STRIDE (Spatial TRanscrIptomics DEconvolution by topic modelling) is a cell-type deconvolution tool for spatial transcriptomics. ",
        license = "GPL-3.0",
        url = "https://github.com/dongqingsun96/STRIDE",
        
        classifiers = [
            "Development Status :: 4 - Beta",
            "Environment :: Console",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
            "Natural Language :: English",
            "Programming Language :: Python :: 3",
            "Topic :: Scientific/Engineering :: Bio-Informatics"
        ],
        python_requires=">=3.7",
    )

if __name__ == "__main__":
    main()

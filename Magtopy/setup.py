#!/usr/bin/env python
"""
The script for building/installing packages
"""

import setuptools

def readreadme(filename):
    with open(filename, "r") as fid:
        long_description = fid.read()
    return long_description

long_description_magtopy = "**Magtopy** - Magnetic Tomography Reconstruction Algorithm in PyCUDA"
packages_list = setuptools.find_packages()

if __name__ == "__main__":
    packages_list = setuptools.find_packages()
    long_description_magtopy = readreadme("./README.md")
    if packages_list is None:
        sys.exit("Failed to fetch packages")

setuptools.setup(
    name="magtopy",
    version="1.0.0",
    author="Marisel Di Pietro Martinez",
    author_email="dipietro AT cpfs DOT mpg DOT de",
    url="https://gitlab.com/magtopy/magtopy",
    package_dir={"magtopy": "magtopy"},
    packages=packages_list,
    description="Magnetic Tomography Reconstruction Algorithm",
    long_description=readreadme("./README.md"),
    license="LICENCE",
    install_requires=[
        "numpy>=1.19.5",
        "scipy>=1.5.4",
        "matplotlib>=3.0.3",
        "pycuda>=2020.1",
        "ipython>=7.5.0",
        "pytest>=4.5.0",
        "pytest-runner>=5.3.0",
        ],
)

print("      ___           ___           ___                       ___           ___")
print("     /__/\         /  /\         /  /\          ___        /  /\         /  /\        ___")
print("    |  |::\       /  /::\       /  /:/_        /  /\      /  /::\       /  /::\      /__/|")
print("    |  |:|:\     /  /:/\:\     /  /:/ /\      /  /:/     /  /:/\:\     /  /:/\:\    |  |:|")
print("  __|__|:|\:\   /  /:/~/::\   /  /:/_/::\    /  /:/     /  /:/  \:\   /  /:/~/:/    |  |:|")
print(" /__/::::| \:\ /__/:/ /:/\:\ /__/:/__\/\:\  /  /::\    /__/:/ \__\:\ /__/:/ /:/   __|__|:|")
print(" \  \:\~~\__\/ \  \:\/:/__\/ \  \:\ /~~/:/ /__/:/\:\   \  \:\ /  /:/ \  \:\/:/   /__/::::\ ")
print("  \  \:\        \  \::/       \  \:\  /:/  \__\/  \:\   \  \:\  /:/   \  \::/       ~\~~\:\ ")
print("   \  \:\        \  \:\        \  \:\/:/        \  \:\   \  \:\/:/     \  \:\         \  \:\ ")
print("    \  \:\        \  \:\        \  \::/          \__\/    \  \::/       \  \:\         \__\/")
print("     \__\/         \__\/         \__\/                     \__\/         \__\/")
print("Magtopy.")

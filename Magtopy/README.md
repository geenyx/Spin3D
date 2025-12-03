# Magnetic Reconstruction from Tomography

PyCUDA library to reconstruct the three components of the magnetization using as input a set of projections from tomographic measurements.

Author: Marisel Di Pietro Martínez

If you use this code, please cite: Di Pietro Martínez et al., Phys. Rev. B 107, 094425 (2023) [https://doi.org/10.1103/PhysRevB.107.094425]

Installation
------------

Installation should be as simple as:

    sudo python3 setup.py install

or, for local installation, using the flag --user:

    python3 setup.py install --user

Test
----

Run

    python setup.py pytest

Dependencies
------------

Magtopy depends on the following python packages:

* python >= 3.6
* numpy
* scipy
* pycuda
* matplotlib

Attention
--------

The library is based on CUDA C, so data MUST be given in [row-major order](https://en.wikipedia.org/wiki/Row-_and_column-major_order).

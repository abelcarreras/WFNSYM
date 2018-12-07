[![Build Status](https://travis-ci.org/abelcarreras/WFNSYM.svg?branch=development)](https://travis-ci.org/abelcarreras/WFNSYM)
[![Coverage Status](https://coveralls.io/repos/github/abelcarreras/WFNSYM/badge.svg?branch=development)](https://coveralls.io/github/abelcarreras/WFNSYM?branch=development)

WFNSYM
=========
Software to calculate continuous symmetry measures of
the electronic wave function of molecules


Installation instructions
---------------------------------------------------------

1. Requirements
  - Lapack & Blas libraries
  - Fortran77 compiler
  - cmake 2.6
  - (optional) MKL


2a. Install as standalone binary
   ```
   ./configure (see --help for available options)
   cd build
   make install
   ```
2b. Compile as a python module
   ```
   cd python
   python setup.py install --user
   ```

Authors
--------------------------------------------------------

This software has been developed by David Casanova
<br>Python module by Abel Carreras

The theoretical background implemented in this software is described in:
<br>Casanova D, Alemany P. Phys Chem Chem Phys. 2010;12(47):15523–9. 
<br>Casanova D, Alemany P, Falceto A, Carreras A, Alvarez S. J Comput Chem 2013;34(15):1321–31.

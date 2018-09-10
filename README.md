
WFNSYM
=========
Software to calculate continuous symmetry measures of
the electronic wave function of molecules


Installation instructions
---------------------------------------------------------

1. Requirements
  - Lapack & Blas libraries
  - cmake 2.6
  - (optional) MKL


2a. Install as standalone binary:
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
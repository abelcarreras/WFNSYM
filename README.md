[![Build Status](https://travis-ci.org/abelcarreras/WFNSYM.svg?branch=master)](https://travis-ci.org/abelcarreras/WFNSYM)
[![Coverage Status](https://coveralls.io/repos/github/abelcarreras/WFNSYM/badge.svg?branch=master)](https://coveralls.io/github/abelcarreras/WFNSYM?branch=master)
[![PyPI version](https://badge.fury.io/py/wfnsympy.svg)](https://badge.fury.io/py/wfnsympy)

WFNSYM
======
This software calculates continuous symmetry measures (CSM) of 
the electronic wave function of molecules


Installation instructions
-------------------------

1a. Requirements
  - LAPACK & Blas libraries
  - Fortran77 compiler (g77/ifort/gfortran)
  - cmake 2.6

1b. Additional requirements for python module
  - Python 2.7.x/3.5+
  - numpy
  - scipy
  - C compiler

2a. Install as standalone binary
   ```shell
   ./configure (see --help for available options)
   cd build
   make install
   ```
2b. Install as a python module
   ```shell
   cd python
   python setup.py install --user
   ```

Simple python API 
-----------------

```python
from wfnsympy import WfnSympy

basis = {'name': 'STO-3G', 
         'primitive_type': 'gaussian', 
         'atoms': [{'symbol': 'O',
                    'shells': [{'shell_type': 's',
                                'p_exponents': [130.70932, 23.808861, 6.4436083],
                                'con_coefficients': [0.154328969, 0.535328136, 0.444634536],
                                'p_con_coefficients': [0.0, 0.0, 0.0]},
                               {'shell_type': 'sp',
                                'p_exponents': [5.0331513, 1.1695961, 0.380389], 
                                'con_coefficients': [-0.0999672287, 0.399512825, 0.700115461], 
                                'p_con_coefficients': [0.155916268, 0.607683714, 0.391957386]}]},  
                   {'symbol': 'H',
                    'shells': [{'shell_type': 's', 
                                'p_exponents': [3.42525091, 0.62391373, 0.1688554], 
                                'con_coefficients': [0.154328971, 0.535328142, 0.444634542], 
                                'p_con_coefficients': [0.0, 0.0, 0.0]}]},     
                   {'symbol': 'H',
                    'shells': [{'shell_type': 's', 
                                'p_exponents': [3.42525091, 0.62391373, 0.1688554],
                                'con_coefficients': [0.154328971, 0.535328142, 0.444634542],
                                'p_con_coefficients': [0.0, 0.0, 0.0]}]}]}

mo_coefficients = [[ 0.9942164, 0.0258468, 0.0000000, 0.0000000,-0.0041640,-0.0055837,-0.0055837],
                   [ 0.2337666,-0.8444565, 0.0000000, 0.0000000, 0.1228297,-0.1555932,-0.1555932],
                   [ 0.0000000, 0.0000000, 0.6126923, 0.0000000, 0.0000000,-0.4492216, 0.4492216],
                   [-0.1040333, 0.5381536, 0.0000000, 0.0000000, 0.7558802,-0.2951071,-0.2951071],
                   [ 0.0000000, 0.0000000, 0.0000000,-1.0000000, 0.0000000, 0.0000000, 0.0000000],
                   [-0.1258185, 0.8201209, 0.0000000, 0.0000000,-0.7635388,-0.7691551,-0.7691551],
                   [ 0.0000000, 0.0000000, 0.9598001, 0.0000000, 0.0000000, 0.8146297,-0.8146297]]

wf_results = WfnSympy(coordinates=[[ 0.00000000, 0.00000000, -0.04280085],
                                   [-0.75810741, 0.00000000, -0.67859957], 
                                   [ 0.75810741, 0.00000000, -0.67859957]],
                      symbols=['O', 'H', 'H'],
                      basis=basis,
                      alpha_mo_coeff=mo_coefficients,
                      charge=0,
                      multiplicity=1,
                      group='C2v')

wf_results.print_CSM()
wf_results.print_ideal_group_table()
wf_results.print_overlap_mo_alpha()
wf_results.print_overlap_wf()
```

Authors
-------

This software has been developed by David Casanova
<br>Python module by Abel Carreras & Efrem Bernuz

The theoretical background implemented in this software is described in:
<br>Casanova D, Alemany P. Phys Chem Chem Phys. 2010;12(47):15523–9. 
<br>Casanova D, Alemany P, Falceto A, Carreras A, Alvarez S. J Comput Chem 2013;34(15):1321–31.

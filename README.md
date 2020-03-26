[![Build Status](https://travis-ci.org/abelcarreras/WFNSYM.svg?branch=development)](https://travis-ci.org/abelcarreras/WFNSYM)
[![Coverage Status](https://coveralls.io/repos/github/abelcarreras/WFNSYM/badge.svg?branch=development)](https://coveralls.io/github/abelcarreras/WFNSYM?branch=development)

WFNSYM
=========
This software allows to calculate continuous symmetry measures of 
the electronic wave function of molecules


Installation instructions
---------------------------------------------------------

1a. Requirements
  - Lapack & Blas libraries
  - Fortran77 compiler (g77/ifort/gfortran)
  - cmake 2.6

1b. Additional requirements for python module
  - Python 2.7.x/3.5+
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

Examples 
--------
**Simple pythonic API to define your input**

```python
from wfnsympy import WfnSympy
from wfnsympy.file_io import get_data_from_file_fchk


basis = {'name': 'STO-3G', 
         'primitive_type': 'gaussian', 
         'atoms': [{'shells': 
                      [{'shell_type': 's', 
                        'functions': 1, 
                        'p_exponents': [130.70932, 23.808861, 6.4436083],
                        'con_coefficients': [0.154328969, 0.535328136, 0.444634536],
                        'p_con_coefficients': [0.0, 0.0, 0.0]},
                       {'shell_type': 'sp', 
                        'functions': 4, 
                        'p_exponents': [5.0331513, 1.1695961, 0.380389], 
                        'con_coefficients': [-0.0999672287, 0.399512825, 0.700115461], 
                        'p_con_coefficients': [0.155916268, 0.607683714, 0.391957386]}], 
                        'symbol': 'O', 
                        'atomic_number': 8}, 
                        
                        {'shells': [{'shell_type': 's', 
                        'functions': 1, 'p_exponents': [3.42525091, 0.62391373, 0.1688554], 
                        'con_coefficients': [0.154328971, 0.535328142, 0.444634542], 
                        'p_con_coefficients': [0.0, 0.0, 0.0]}], 
                        'symbol': 'H', 'atomic_number': 1}, 
                        
                        {'shells': [{'shell_type': 's', 
                        'functions': 1,
                        'p_exponents': [3.42525091, 0.62391373, 0.1688554],
                        'con_coefficients': [0.154328971, 0.535328142, 0.444634542],
                        'p_con_coefficients': [0.0, 0.0, 0.0]}],
                        'symbol': 'H', 'atomic_number': 1}]}


molecule = get_data_from_file_fchk('molecule.fchk')
wf_results = WfnSympy(coordinates=molecule['coordinates'],
                      symbols=molecule['symbols'],
                      basis=molecule['basis'],
                      alpha_mo_coeff=molecule['mo_coefficients']['alpha'],
                      charge=0,
                      multiplicity=1,
                      group='C6v')

wf_results.print_CSM()
wf_results.print_ideal_group_table()
wf_results.print_overlap_mo_alpha()
wf_results.print_overlap_wf()
```

Authors
--------------------------------------------------------

This software has been developed by David Casanova
<br>Python module by Abel Carreras & Efrem Bernuz

The theoretical background implemented in this software is described in:
<br>Casanova D, Alemany P. Phys Chem Chem Phys. 2010;12(47):15523–9. 
<br>Casanova D, Alemany P, Falceto A, Carreras A, Alvarez S. J Comput Chem 2013;34(15):1321–31.

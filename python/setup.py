from numpy.distutils.core import setup, Extension
from distutils.dir_util import copy_tree
from distutils.errors import DistutilsFileError
import os

# Make python package


travis = bool('TRAVIS' in os.environ)
if travis:
    print('Testing in travis')
    copy_tree('../src', './src', update=True)
    copy_tree('../include', './include', update=True)

s_dir = 'src/'

wfnsymlib = Extension('wfnsympy.WFNSYMLIB',
                      # extra_compile_args=['-std=c99'],
                      #include_dirs=include_dirs_numpy,
                      include_dirs=['include'],
                      libraries=['lapack', 'blas'],
                      sources=['WFNSYMLIB.pyf',
                               s_dir + 'VRoutines.F',
                               s_dir + 'aos_product.F',
                               s_dir + 'get_basis_lib.F',
                               s_dir + 'get_dim.F',
                               s_dir + 'get_mos.F',
                               s_dir + 'group_dim.F',
                               s_dir + 'group_table.F',
                               s_dir + 'make_molden.F',
                               s_dir + 'make_s.F',
                               s_dir + 'make_uhfno.F',
                               s_dir + 'print_routines.F',
                               s_dir + 'lib_main.F',
                               s_dir + 'overlap.F',
                               s_dir + 'read_routines.F',
                               s_dir + 'sym_routines.F'])

setup(name='wfnsympy',
      version='0.1',
      description='wfnsympy',
      author='Abel Carreras',
      author_email='abelcarreras83@gmail.com',
      packages=['wfnsympy'],
      ext_modules=[wfnsymlib])

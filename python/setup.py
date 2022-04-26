from numpy.distutils.core import setup, Extension
from distutils.dir_util import copy_tree
from distutils.errors import DistutilsFileError
import os, sys
from shutil import copyfile


def get_version_number():
    main_ns = {}
    for line in open('wfnsympy/__init__.py', 'r').readlines():
        if not(line.find('__version__')):
            exec(line, main_ns)
            return main_ns['__version__']


if bool('TRAVIS_WFNSYM' in os.environ):
    print('Testing in travis')
    copy_tree('../src', './src', update=True)
    copy_tree('../include', './include', update=True)
    copyfile(os.path.abspath('../README.md'), os.path.abspath('./README.md'))

if os.path.isfile('../README.md') and not os.path.isfile('./README.md'):
    copyfile(os.path.abspath('../README.md'), os.path.abspath('./README.md'))

if os.path.isdir('./src') and os.path.isdir('./include'):
    print('exists')
    s_dir = './src/'
    i_dir = './include/'
else:
    print('not exists')
    s_dir = '../src/'
    i_dir = '../include/'

wfnsymlib = Extension('wfnsympy.WFNSYMLIB',
                      # extra_compile_args=['-ffixed-form', '-ffixed-line-length-none'],
                      # include_dirs=include_dirs_numpy,
                      extra_f77_compile_args=['-ffixed-line-length-0'],
                      include_dirs=[i_dir],
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

qsymlib = Extension('wfnsympy.QSYMLIB',
                    extra_f77_compile_args=['-ffixed-line-length-0'],
                    include_dirs=[i_dir],
                    libraries=['lapack', 'blas'],
                    sources=['QSYMLIB.pyf',
                             s_dir + 'center_dens.F',
                             s_dir + 'dens_basis_lib.F',
                             s_dir + 'make_dens.F',
                             s_dir + 'group_dim.F',
                             s_dir + 'norma.F',
                             s_dir + 'overlap.F',
                             s_dir + 'denslib.F',
                             s_dir + 'self_similarity.F',
                             s_dir + 'sym_transform.F',
                             s_dir + 'sym_overlap.F',
                             s_dir + 'sym_routines.F',
                             s_dir + 'VRoutines.F'])

setup(name='wfnsympy',
      version=get_version_number(),
      description='wfnsympy',
      # long_description=open('README.md').read(),
      # long_description_content_type='text/markdown',
      author='Abel Carreras',
      author_email='abelcarreras83@gmail.com',
      install_requires=['numpy', 'scipy'],
      packages=['wfnsympy'],
      include_package_data=True,
      ext_modules=[wfnsymlib, qsymlib])

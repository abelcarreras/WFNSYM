from numpy.distutils.core import setup, Extension


wfnsymlib = Extension('wfnsympy.WFNSYMLIB',
                      # extra_compile_args=['-std=c99'],
                      #include_dirs=include_dirs_numpy,
                      include_dirs=['../include'],
                      libraries=['lapack', 'blas'],
                      sources=['WFNSYMLIB.pyf',
                               '../src/VRoutines.F',
                               '../src/aos_product.F',
                               '../src/get_basis_lib.F',
                               '../src/get_dim.F',
                               '../src/get_mos.F',
                               '../src/group_dim.F',
                               '../src/group_table.F',
                               '../src/make_molden.F',
                               '../src/make_s.F',
                               '../src/make_uhfno.F',
                               '../src/print_routines.F',
                               '../src/lib_main.F',
                               '../src/overlap.F',
                               '../src/read_routines.F',
                               '../src/sym_routines.F'])

setup(name='wfnsympy',
      version='0.1',
      description='wfnsympy',
      author='Abel Carreras',
      author_email='abelcarreras83@gmail.com',
      packages=['wfnsympy'],
      ext_modules=[wfnsymlib])

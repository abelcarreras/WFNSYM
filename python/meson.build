project('wfnsymlib', 'c', 'fortran', version: '0.1.0')

py = import('python').find_installation(pure: false)
py_dep = py.dependency()

cc = meson.get_compiler('c')
# blas_dep = cc.find_library('openblas', required : true)
blas_dep = dependency('blas', required : true)
lapack_dep = dependency('lapack', required : true)

#locate libraries quadmath (needed for windows)
quadmath_dep = meson.get_compiler('fortran').find_library('quadmath')


# fortran source files
fortran_sources = [
    'WFNSYMLIBmodule.c',
    'src/VRoutines.F',
    'src/aos_product.F',
    'src/get_basis_lib.F',
    'src/get_dim.F',
    'src/get_mos.F',
    'src/group_dim.F',
    'src/group_table.F',
    'src/make_molden.F',
    'src/make_s.F',
    'src/make_uhfno.F',
    'src/print_routines.F',
    'src/lib_main.F',
    'src/overlap.F',
    'src/read_routines.F',
    'src/sym_routines.F'
]

# fortran source files
fortran_sources_2 = [
    'QSYMLIBmodule.c',
    'src/center_dens.F',
    'src/dens_basis_lib.F',
    'src/make_dens.F',
    'src/group_dim.F',
    'src/norma.F',
    'src/overlap.F',
    'src/denslib.F',
    'src/self_similarity.F',
    'src/sym_transform.F',
    'src/sym_overlap.F',
    'src/sym_routines.F',
]



# Fortran compiler option
fortran_flags = ['-ffixed-line-length-0', '-O2']

incdir_numpy = run_command(py,
  ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()

incdir_f2py = run_command(py,
    ['-c', 'import os; os.chdir(".."); import numpy.f2py; print(numpy.f2py.get_include())'],
    check : true
).stdout().strip()

inc_np = include_directories(incdir_numpy, incdir_f2py, '../include')

wfnsymlib = py.extension_module('WFNSYMLIB',
                                fortran_sources,
                                # incdir_f2py / 'fortranobject.c',
                                join_paths(incdir_f2py, 'fortranobject.c'),
                                include_directories: inc_np,
                                dependencies : [py_dep, quadmath_dep, blas_dep, lapack_dep],
                                fortran_args: fortran_flags,
                                install : true,
                                # link_args: ['-static-libgfortran', '-Wl,-rpath,@loader_path/../wfnsymlib.dylibs'],
                                install_dir : get_option('prefix')
)


qsymlib = py.extension_module('QSYMLIB',
                               fortran_sources_2,
                               # incdir_f2py / 'fortranobject.c',
                               join_paths(incdir_f2py, 'fortranobject.c'),
                               include_directories: inc_np,
                               dependencies : [py_dep, quadmath_dep, blas_dep, lapack_dep],
                               fortran_args: fortran_flags,
                               install : true,
                               # link_args: ['-static-libgfortran', '-Wl,-rpath,@loader_path/../wfnsymlib.dylibs'],
                               install_dir : get_option('prefix')
)


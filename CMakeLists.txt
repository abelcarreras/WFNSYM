# CMake project file for WFNSYM

cmake_minimum_required (VERSION 2.6)
project (WFNSYM)
enable_language (Fortran)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# default installation
get_filename_component (default_prefix "." ABSOLUTE)
set (CMAKE_INSTALL_PREFIX ${default_prefix} CACHE STRING
      "Choose the installation directory; by default it installs in the WFNSYM directory."
      FORCE)

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  # gfortran
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -fbacktrace")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  # ifort
  set (CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O0 -g -traceback -debug extended")
elseif (Fortran_COMPILER_NAME MATCHES "g77")
  # g77
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3 -m32")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -m32")
else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")

include_directories("${CMAKE_SOURCE_DIR}/include")

# MKL library
if(WITH_MKL)
  message(STATUS "Use MKL library")
  link_directories("/opt/intel/mkl/lib/")
  include_directories("/opt/intel/compilers_and_libraries/mac/include/intel64/")
#  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -lmkl_def -lmkl_core -lmkl_lapack -lmkl_intel_thread  -lmkl_sequential -lmkl_intel_lp64 -lguide  -lpthread")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -lmkl_core -lmkl_intel_thread -lmkl_sequential -lmkl_intel_lp64 -lpthread")

else(WITH_MKL)
  link_directories("/usr/lib")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -llapack -lblas")   
endif(WITH_MKL)


# build executables
set( EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/EXE)
set( LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/LIB)
set (EXECUTABLES "WFNSYM")

add_executable ("WFNSYM"
                "src/VRoutines.F"
                "src/aos_product.F"
                "src/get_basis.F"
                "src/get_dim.F"
                "src/get_mos.F"
                "src/group_dim.F"
                "src/group_table.F"
                "src/make_molden.F"
                "src/make_s.F"
                "src/make_uhfno.F"
                "src/overlap.F"
                "src/print_routines.F"
                "src/prog_main.F"
                "src/read_routines.F"
                "src/sym_routines.F")

ADD_LIBRARY ("WFNSYMLIB" "STATIC"
                "src/VRoutines.F"
                "src/aos_product.F"
                "src/get_basis_lib.F"
                "src/get_dim.F"
                "src/get_mos.F"
                "src/group_dim.F"
                "src/group_table.F"
                "src/make_molden.F"
                "src/make_s.F"
                "src/make_uhfno.F"
                "src/overlap.F"
                "src/print_routines.F"
                "src/lib_main.F"
                "src/read_routines.F"
                "src/sym_routines.F")


# install executables and scripts
install (TARGETS ${EXECUTABLES}
         RUNTIME DESTINATION "bin")

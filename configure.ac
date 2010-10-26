AC_PREREQ([2.63])
AC_INIT(libplump, 0.1, j.gasthaus@ucl.ac.uk)
AM_INIT_AUTOMAKE([1.10 no-define])
LT_INIT
LT_LANG([C++])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_FILES([Makefile
                 src/Makefile
                 bindings/Makefile 
                 bindings/python/Makefile
                 bindings/python/setup.py])

AC_ARG_WITH([python],
              [AS_HELP_STRING([--with-python],
                [use SWIG to generate a Python binding])],
                [],
                [pythonbinding=no])

AC_PATH_PROG([SWIG], [swig])
AM_CONDITIONAL([HAVE_SWIG], [test "$SWIG"])
if test "$SWIG" -a "x$pythonbinding" != xno
then
    AX_PKG_SWIG(1.3, [], [])
    SWIG_ENABLE_CXX
    SWIG_PYTHON
fi

AM_CONDITIONAL([PYTHONBINDING], test "$SWIG" -a "x$pythonbinding" != xno)

if test  "x$pythonbinding" != xno -a ! "$SWIG"
then
    AC_MSG_WARN([Cannot build Python binding -- SWIG not found.])
fi

# Checks for programs.
AC_PROG_CXX

# Checks for libraries.
AX_BOOST_BASE([1.3.5],[],[])
AX_BOOST_PROGRAM_OPTIONS
AX_BOOST_FILESYSTEM
AX_BOOST_SERIALIZATION
AX_BOOST_IOSTREAMS

# Checks for header files.
AC_CHECK_HEADERS([stdint.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_INT32_T
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_LIB([m],[cos])
AC_CHECK_LIB([gslcblas],[cblas_dgemm])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm])

AC_OUTPUT

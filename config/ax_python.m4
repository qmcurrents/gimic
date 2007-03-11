dnl @synopsis AX_PYTHON
dnl
dnl This macro does a complete Python development environment check.
dnl
dnl It recurses through several python versions (from 2.1 to 2.4 in
dnl this version), looking for an executable. When it finds an
dnl executable, it looks to find the header files and library.
dnl
dnl It sets PYTHON_BIN to the name of the python executable,
dnl PYTHON_INCLUDE_DIR to the directory holding the header files, and
dnl PYTHON_LIB to the name of the Python library.
dnl
dnl This macro calls AC_SUBST on PYTHON_BIN (via AC_CHECK_PROG),
dnl PYTHON_INCLUDE_DIR and PYTHON_LIB.
dnl
dnl @category InstalledPackages
dnl @author Michael Tindal <mtindal@paradoxpoint.com>
dnl @version 2004-09-20
dnl @license GPLWithACException

AC_DEFUN([AX_PYTHON],
[AC_MSG_CHECKING(for python build information)
AC_MSG_RESULT([])
for python in python2.4 python2.3 python2.2 python2.1 python; do
AC_CHECK_PROGS(PYTHON_BIN, [$python])
ax_python_bin=$PYTHON_BIN
if test x$ax_python_bin != x; then
   AC_CHECK_LIB($ax_python_bin, main, ax_python_lib=$ax_python_bin, ax_python_lib=no)
   AC_CHECK_HEADER([$ax_python_bin/Python.h],
   [[ax_python_header=`locate $ax_python_bin/Python.h | sed -e s,/Python.h,,`]],
   ax_python_header=no)
   if test $ax_python_lib != no; then
     if test $ax_python_header != no; then
       break;
     fi
   fi
fi
done
if test x$ax_python_bin = x; then
   ax_python_bin=no
fi
if test x$ax_python_header = x; then
   ax_python_header=no
fi
if test x$ax_python_lib = x; then
   ax_python_lib=no
fi

AC_MSG_RESULT([  results of the Python check:])
AC_MSG_RESULT([    Binary:      $ax_python_bin])
AC_MSG_RESULT([    Library:     $ax_python_lib])
AC_MSG_RESULT([    Include Dir: $ax_python_header])

if test x$ax_python_header != xno; then
  PYTHON_INCLUDE_DIR=$ax_python_header
  AC_SUBST(PYTHON_INCLUDE_DIR)
fi
if test x$ax_python_lib != xno; then
  PYTHON_LIB=$ax_python_lib
  AC_SUBST(PYTHON_LIB)
fi
])dnl

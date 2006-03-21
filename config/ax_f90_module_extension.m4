dnl @synopsis AX_F90_MODULE_EXTENSION
dnl
dnl Find Fortran 90 modules file extension. The module extension is
dnl stored in the cached variable ax_f90_modext, or "unknown" if the
dnl extension cannot be found.
dnl
dnl @category Fortran
dnl @author Luc Maisonobe <luc@spaceroots.org>
dnl @version 2005-06-17
dnl @license AllPermissive

AC_DEFUN([AX_F90_MODULE_EXTENSION],[
AC_CACHE_CHECK([fortran 90 modules extension],
ax_f90_modext,
[AC_LANG_PUSH(Fortran)
i=0
while test \( -f tmpdir_$i \) -o \( -d tmpdir_$i \) ; do
  i=`expr $i + 1`
done
mkdir tmpdir_$i
cd tmpdir_$i
AC_COMPILE_IFELSE([module conftest_module
   contains
   subroutine conftest_routine
   write(*,'(a)') 'gotcha!'
   end subroutine conftest_routine
   end module conftest_module
  ],
  [ax_f90_modext=`ls | sed -n 's,conftest_module\.,,p'`
   ax_f90_modcase=lowercase
   if test x$ax_f90_modext = x ; then
     ax_f90_modcase=uppercase
dnl Some F90 compilers put module filename in uppercase letters
     ax_f90_modext=`ls | sed -n 's,CONFTEST_MODULE\.,,p'`
     if test x$ax_f90_modext = x ; then
       ax_f90_modext=unknown
	   ax_f90_modcase=unknown
     fi
   fi
  ],
  [ax_f90_modext=unknown; ax_f90_modcase=unknown])
cd ..
rm -fr tmpdir_$i
FC_MODEXT=$ax_f90_modext
FC_MODCASE=$ax_f90_modcase
AC_SUBST([FC_MODEXT])
AC_SUBST([FC_MODCASE])
AC_LANG_POP(Fortran)
])])

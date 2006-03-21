m4_include([config/acx_blas.m4])
m4_include([config/acx_lapack.m4])
m4_include([config/acx_blas95.m4])
m4_include([config/acx_lapack95.m4])
m4_include([config/ax_fc_search_path.m4])
m4_include([config/acx_fc_library_setup.m4])
m4_include([config/ax_f90_header.m4])
m4_include([config/ax_f90_internal_headmod.m4])
m4_include([config/ax_f90_library.m4])
m4_include([config/ax_f90_library_setup.m4])
m4_include([config/ax_f90_module_extension.m4])
m4_include([config/ax_f90_module_flag.m4])
m4_include([config/ax_f90_module.m4])
m4_include([config/check_gnu_make.m4])
m4_include([config/acx_mpi.m4])
m4_include([config/acx_getkw.m4])
m4_include([config/acx_f90_mpi.m4])

AC_DEFUN([ACX_BUILD_FLAGS],[. ./config/$1.conf])

AC_DEFUN([ACX_SUBST_BUILD_FLAGS],[
AC_SUBST(fcflags, $fcflags)
AC_SUBST(fdebug, $fdebug)
AC_SUBST(fprof, $fprof)
AC_SUBST(frange, $frange)
AC_SUBST(fldflags, $fldflags)

AC_SUBST(cflags, $cflags)
AC_SUBST(cdebug, $cdebug)
AC_SUBST(cprof, $cprof)
AC_SUBST(crange, $crange)
AC_SUBST(cldflags, $cldflags)

]) #ACX_SUBST_BUILD_FLAGS

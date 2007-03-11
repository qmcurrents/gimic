#
# Process this file with autoconf to produce a configure script.
#
# Written by Jonas Juselius <jonas@iki.fi>
#

AC_INIT([GIMIC],[1.0],[jonas@iki.fi])
AC_CONFIG_SRCDIR(gimic.f90)
AC_CONFIG_HEADER(config.h)
AC_CONFIG_AUX_DIR(config)
AC_PREFIX_DEFAULT([/usr/local])
AC_PROG_MAKE_SET
AC_CANONICAL_HOST

AC_ARG_ENABLE([mpi], AS_HELP_STRING(--enable-mpi,Enable the parallel version))

#AC_ARG_WITH([aces2], AS_HELP_STRING(--with-aces2=/usr/local/aces2/lib,path to the ACES2 libraries))

# Check for programs.

CHECK_GNU_MAKE
if test "$ifGNUmake" = "#"; then
	AC_MSG_ERROR([No working version of GNU make found!])
fi
AC_PROG_RANLIB
AC_PROG_INSTALL
AC_CHECK_PROGS(LD,ld,false)
AC_CHECK_PROGS(AR,ar,false)
AC_CHECK_PROGS(STRIP,strip,true)
AC_CHECK_PROGS(ETAGS,etags,true)
AC_CHECK_PROGS(CTAGS,ctags,true)
AC_PROG_LN_S

AC_CACHE_SAVE

# ACX_ macros are defined in aclocal.m4
# Setup compilers for different architechtures
case $host in 
	*i*86*-linux*)
		try_cc="gcc icc cc"
		try_fort="ifort ifc pgf90"
		build_type="i386_linux"
		;;
	*ia*64*-linux*)
		try_cc="gcc icc cc"
		try_fort="ifort ifc gfortran pgf90 f90"
		build_type="ia64_linux"
		;;
	*alpha*-linux*)
		try_cc="ccc cc gcc"
		try_fort="fort f90"
		build_type="linux_alpha"
		;;
	*dec*-osf*)
		try_cc="ccc cc gcc"
		try_fort="fort f90"
		build_type="dec_osf"
		;;
	*ibm*-aix*)
		try_cc="cc gcc"
		try_fort="xlf90 f90"
		build_type="ibm_aix"
		;;
	*sgi*-irix*)
		try_cc="cc gcc"
		try_fort="f90"
		build_type="sgi_irix"
		;;
	*hppa*-hpux*)
		try_cc="cc gcc"
		try_fort="f90"
		build_type="hp_hpux"
		;;
	*sparc*-solaris*)
		try_cc="cc gcc"
		try_fort="f90"
		build_type="sun_solaris"
		;;
	*) 	
		AC_MSG_WARN([Unknown architecture! This might not work...])
		try_cc="cc gcc ccc cl"
		try_fort="f90 xlf90 ifc pgf90 epcf90 f95 xlf95 lf95 fort"
		build_type="unknown"
		;;
esac

AC_PROG_CC([$try_cc])
AC_PROG_FC([$try_fort])

# we build these from scratch, to avoid stupid defaults...
CFLAGS=""
FCFLAGS=""

ACX_BUILD_FLAGS([$build_type])
AC_CACHE_SAVE

AX_FC_SEARCH_PATH

AC_FC_LIBRARY_LDFLAGS
AC_FC_WRAPPERS

AC_LANG_PUSH([Fortran])
AC_FC_SRCEXT([f90])
AC_LANG_POP([Fortran])

AC_CACHE_SAVE
AC_LANG_PUSH([Fortran])
ACX_BLAS(,AC_MSG_WARN([No working blas library found!]))
ACX_LAPACK(,AC_MSG_WARN([No working lapack library found!]))
AC_LANG_POP([Fortran])
ACX_BLAS95(,AC_MSG_WARN([No working lapack95 library found!]))
ACX_LAPACK95(,AC_MSG_WARN([No working lapack95 library found!]))

ax_save_LIBS="$LIBS"
AC_MSG_NOTICE([Checking if ACES2 libraries are available])
AX_F90_LIBRARY_SETUP([aces2],[],[],[liblibr.a],[call crapsi], [$FC_SEARCH_PATH])
LIBS="$ax_save_LIBS"

if test "x$enable_mpi" = "xyes"; then
	AC_CACHE_SAVE
	AC_MSG_NOTICE([Checking if MPI is available])
	AC_LANG_PUSH([Fortran])
	ACX_MPI([FC=$MPIFC])
	AC_LANG_POP
	fcflags="$fcflags -DHAVE_MPI"
fi
AC_CACHE_SAVE

AX_PYTHON

AC_CACHE_SAVE

ACX_SUBST_BUILD_FLAGS 
#AC_SUBST(ACES2_LIBS, $ACES2_LIBS)
AC_SUBST(top_bindir,['$(topdir)'])
AC_SUBST(top_includedir,['$(topdir)/include'])
AC_SUBST(top_libdir,['$(topdir)/lib'])

AC_CONFIG_FILES([Makefile make.config make.rules gimic])
AC_OUTPUT

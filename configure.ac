AC_INIT([seqwise], [m4_esyscmd_s([git describe --always])])
AC_CONFIG_AUX_DIR([.])
AM_INIT_AUTOMAKE
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_RANLIB
AM_SILENT_RULES([no])

m4_include([m4/ax_lib_hdf5.m4])
m4_include([m4/ax_openmp.m4])


tldevel_check_f=$ac_abs_confdir/tldevel/README

AC_CHECK_FILE([$tldevel_check_f], [],[have_tldevel=no])


if test "x${have_tldevel}" = xno; then
echo "------------------------------------------"
echo " Unable to find pthreads on this system.  "
echo " Building a single-threaded version.      "
echo "------------------------------------------"
AC_MSG_ERROR([ASFSAF])
fi



AC_ARG_ENABLE(openmp, [AC_HELP_STRING([--enable-openmp],[use OpenMP directives for parallelism])], enable_openmp=$enableval, enable_openmp=no)

if test "$enable_openmp" = "yes"; then
AC_DEFINE(HAVE_OPENMP,1,[Define to enable OpenMP])
AX_OPENMP([], [AC_MSG_ERROR([don't know how to enable OpenMP])])
fi


AC_ARG_WITH([openmp],
AS_HELP_STRING([--without-openmp], [Do not use OpenMP]))
AS_IF([test "x$with_openmp" = "xno"],
[],
[AX_OPENMP(
AC_DEFINE([HAVE_OPENMP],[1],[Defined if OpenMP should and can be used])
)]
)
AC_SUBST(OPENMP_CFLAGS)
# according to http://beagle-lib.googlecode.com/svn-history/r824/trunk/configure.ac:
# "OpenMP checker only defines for C when compiling both C and C++"
OPENMP_CXXFLAGS=$OPENMP_CFLAGS
AC_SUBST(OPENMP_CXXFLAGS)



PKG_PROG_PKG_CONFIG([0.25]) # check and set $PKG_CONFIG

PKG_CHECK_MODULES([CAIRO],    [cairo],    [enable_cairo=yes],    [AS_IF([test "$enable_cairo" = "yes"],
      [AC_MSG_ERROR([cairo required for plotting, but not found.])],
      [enable_cairo=no])])
AM_CONDITIONAL([ENABLE_CAIRO], [test "$enable_cairo" = "yes"])



# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [log])
# FIXME: Replace `main' with a function in `-lpthread':
AC_CHECK_LIB([pthread], [pthread_create])

AC_CHECK_HEADERS([stdlib.h string.h float.h limits.h stddef.h stdlib.h sys/time.h sys/timeb.h unistd.h])



LIB_TLDEVELDIR="tldevel"
AC_SUBST(LIB_TLDEVELDIR)

AC_FUNC_ERROR_AT_LINE

AC_TYPE_SIZE_T

AC_TYPE_INT8_T
AC_TYPE_INT32_T
AC_TYPE_INT64_T

AC_TYPE_UINT8_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T



AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([arc4random arc4random_uniform])
AC_CHECK_FUNCS([floor gettimeofday pow sqrt getcwd mkdir])
AC_FUNC_LSTAT_FOLLOWS_SLASHED_SYMLINK


AC_CANONICAL_HOST
case $host_os in
darwin* )
DEAD_STRIP="-Wl,-dead_strip"
;;
*)
DEAD_STRIP="-Wl,--gc-sections -Wl,--as-needed"
;;
esac
AC_SUBST(DEAD_STRIP)




AC_ARG_ENABLE(debugging,[AS_HELP_STRING([--enable-debugging],[include debugging code])
AS_HELP_STRING([--enable-debugging=x],[also set diagnostics verbosity level to <x> (1-3)])],
enable_debugging=$enableval, enable_debugging="no")

case $enable_debugging in
yes)  AC_DEFINE(DEBUGLEVEL, 0,[No debugging. ]);;
1)  AC_DEFINE(DEBUGLEVEL, 1,[Defines debugging level 1.]);;
2)  AC_DEFINE(DEBUGLEVEL, 2,[Defines debugging level 2.]);;
3)  AC_DEFINE(DEBUGLEVEL, 3,[Defines debugging level 3.]);;
no)  AC_DEFINE(DEBUGLEVEL, 0,[No debugging.]);;
*)  AC_MSG_ERROR([Unknown argument to --enable-debugging: $enable_debugging]);;
esac

if test "$enable_debugging" != "no"; then
AC_DEFINE(DEBUG,1,[Defines debugging .])
CFLAGS="-ggdb -pedantic -Wall -std=gnu11 -fdata-sections -ffunction-sections"
else
CFLAGS="-O3 -pedantic -funroll-loops -Wall -std=gnu11 -fdata-sections -ffunction-sections"
DEBUG=0
fi

ac_configure_args="$ac_configure_args PACKAGE_VERSION=$PACKAGE_VERSION   PACKAGE_NAME=$PACKAGE_NAME"

AC_ARG_ENABLE(valgrind-tests,
[--enable-valgrind-tests run self tests under valgrind],
opt_valgrind_tests="$enableval", opt_valgrind_tests="no")


OPTS="-q --trace-children=no --error-exitcode=1 --leak-check=yes"
# Run self-tests under valgrind?
if test "$opt_valgrind_tests" = "yes"; then
AC_DEFINE(DEBUG,1,[Defines debugging .])
CFLAGS="-ggdb -pedantic -Wall -m64 -std=gnu11"
AC_CHECK_PROGS(VALGRIND, valgrind)
VALGRIND="$VALGRIND $OPTS"
else
VALGRIND=
fi
AC_MSG_CHECKING([whether self tests are run under valgrind])
AC_MSG_RESULT($opt_valgrind_tests)


AX_LIB_HDF5([serial])
if test "$with_hdf5" = "yes"; then
CC="$HDF5_CC"
else
AC_MSG_ERROR([

Unable to find HDF5, we need parallel HDF5.
Install on fedora:
sudo dnf install hdf5 hdf5-devel

])
fi


ac_configure_args="$ac_configure_args PACKAGE_VERSION=$PACKAGE_VERSION   PACKAGE_NAME=$PACKAGE_NAME"

AC_CONFIG_FILES([Makefile
src/Makefile])
AC_CONFIG_SUBDIRS([tldevel])

AC_OUTPUT



if test "$enable_cairo" = no; then
echo "------------------------------------------"
echo " Unable to find cairo  on this system.    "
echo " Motif plotting code will not be built.   "
echo " To install cairo:                        "
echo " Fedora: sudo dnf install cairo-devel     "
echo " Ubuntu: sudo apt install libcairo2-dev   "
echo "------------------------------------------"
fi


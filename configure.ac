AC_INIT([seqer], [0.8.0])

AC_CONFIG_AUX_DIR([.])
AM_INIT_AUTOMAKE


AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_LIBTOOL
dnl AC_PROG_RANLIB
AM_SILENT_RULES([no])

# Checks for programs.
AC_PROG_CXX
AC_PROG_AWK
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET

LT_INIT

AC_CONFIG_MACRO_DIR([m4])

m4_include([m4/m4_ax_valgrind_check.m4])
m4_include([m4/ax_check_compile_flag.m4])
m4_include([m4/ax_with_tldevel.m4])
m4_include([m4/ax_set_compile_flags.m4])
m4_include([m4/ax_lib_hdf5.m4])
m4_include([m4/ax_openmp.m4])

dnl AC_CONFIG_HEADERS([src/config_seqer.h])

AX_ENABLE_DEBUG

AX_LIB_HDF5([serial])
if test "$with_hdf5" = "yes"; then

dnl CC="$HDF5_CC"
CFLAGS="${CFLAGS} ${HDF5_CFLAGS}"
CPPFLAGS="${CPPFLAGS} ${HDF5_CPPFLAGS}"
LDFLAGS="${LDFLAGS} ${HDF5_LDFLAGS}"
LIBS="${LIBS} ${HDF5_LIBS}"
else
AC_MSG_WARN([

Unable to find HDF5, we need parallel HDF5.
Install on fedora:
sudo dnf install hdf5 hdf5-devel

])
fi

AM_CONDITIONAL([ENABLE_HDF5], [test "$with_hdf5" = "yes"])

# test for working C11 _Generic support
AN_IDENTIFIER([_Generic], [AC_C__GENERIC])
AC_DEFUN([AC_C__GENERIC],
[AC_CACHE_CHECK([for _Generic], ac_cv_c__Generic,
[AC_COMPILE_IFELSE(
   [AC_LANG_SOURCE(
      [[int
        main (int argc, char **argv)
        {
          int a = _Generic (argc, int: argc = 1);
          int *b = &_Generic (argc, default: argc);
          char ***c = _Generic (argv, int: argc, default: argv ? &argv : 0);
          _Generic (1 ? 0 : b, int: a, default: b) = &argc;
          _Generic (a = 1, default: a) = 3;
          return a + !b + !c;
        }
      ]])],
   [ac_cv_c__Generic=yes],
   [ac_cv_c__Generic=no])])
if test $ac_cv_c__Generic = no; then
echo "-------------------------------------------"
echo " Unable to compile C11 code on this system!"
echo "-------------------------------------------"
echo ""
echo "Install a new(er) version of GCC (e.g. GCC 4.9+)"
echo "or clang!"
echo ""
AC_MSG_ERROR([ERROR: C compiler too old!])

AC_MSG_ERROR([We could not compile a simple C11 example. Install a new(er) compiler version (e.g.  GCC 4.9).])
fi
])

# Check for C11 _Generic support
AC_C__GENERIC

AC_C_INLINE

AX_WITH_TLDEVEL

AC_MSG_NOTICE([$TLDEVEL_CPPFLAGS])
AC_MSG_NOTICE([$TLDEVEL_LDFLAGS])
AC_MSG_NOTICE([$TLDEVEL_LIB])



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



dnl PKG_PROG_PKG_CONFIG([0.25]) # check and set $PKG_CONFIG

dnl PKG_CHECK_MODULES([CAIRO],    [cairo],    [enable_cairo=yes],    [AS_IF([test "$enable_cairo" = "yes"],
dnl       [AC_MSG_ERROR([cairo required for plotting, but not found.])],
dnl       [enable_cairo=no])])
dnl AM_CONDITIONAL([ENABLE_CAIRO], [test "$enable_cairo" = "yes"])



# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [log])
# FIXME: Replace `main' with a function in `-lpthread':
AC_CHECK_LIB([pthread], [pthread_create])

AC_CHECK_HEADERS([stdlib.h string.h float.h limits.h stddef.h stdlib.h sys/time.h sys/timeb.h unistd.h])

AC_FUNC_ERROR_AT_LINE

AC_TYPE_SIZE_T

AC_TYPE_INT8_T
AC_TYPE_INT16_T
AC_TYPE_INT32_T
AC_TYPE_INT64_T

AC_TYPE_UINT8_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T


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


AX_VALGRIND_DFLT([memcheck], [on])
AX_VALGRIND_CHECK()

AC_SUBST([DATE]) DATE="`date`"

AC_SUBST([USER]) USER="`whoami`@`hostname`"

AC_SUBST([UNAME_INFO]) UNAME_INFO=`uname -n`

AC_SUBST([COMPILER]) COMPILER=`$CC --version | head -n 1`


AC_CONFIG_FILES([Makefile
src/Makefile])

AC_OUTPUT

AC_MSG_NOTICE([---------------------------------------------------------])
AC_MSG_NOTICE([Configuration complete - $PACKAGE_NAME-$PACKAGE_VERSION])
AC_MSG_NOTICE([])
AC_MSG_NOTICE([General Information:])
AC_MSG_NOTICE([               Version: $VERSION])
AC_MSG_NOTICE([         Configured on: $DATE])
AC_MSG_NOTICE([         Configured by: $USER])
AC_MSG_NOTICE([           Host system: $host_cpu-$host_vendor-$host_os])
AC_MSG_NOTICE([     Uname information: $UNAME_INFO])
AC_MSG_NOTICE([                CFLAGS: $CFLAGS])
AC_MSG_NOTICE([               CPFLAGS: $CPPFLAGS   $TLDEVEL_CPPFLAGS])
AC_MSG_NOTICE([                  LIBS: $LIBS $TLDEVEL_LIB])
AC_MSG_NOTICE([               LDFLAGS: $LDFLAGS $TLDEVEL_LDFLAGS])
AC_MSG_NOTICE([    Installation point: $prefix])
AC_MSG_NOTICE([---------------------------------------------------------])



if test "$enable_cairo" = no; then
echo "------------------------------------------"
pecho " Unable to find cairo  on this system.    "
echo " Motif plotting code will not be built.   "
echo " To install cairo:                        "
echo " Fedora: sudo dnf install cairo-devel     "
echo " Ubuntu: sudo apt install libcairo2-dev   "
echo "------------------------------------------"
fi


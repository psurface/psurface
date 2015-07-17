AC_INIT([psurface],[2.0.0],[sander@igpm.rwth-aachen.de])
AC_PREREQ([2.60]) # FIXME: what version do we really need?
AC_CONFIG_SRCDIR([.])
AC_CONFIG_AUX_DIR([build-aux])
AM_INIT_AUTOMAKE([-Wall -Werror foreign silent-rules])
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])   # needed for some non-standard archivers
LT_INIT
AC_PROG_CXX
AX_CXX_COMPILE_STDCXX_11(noext, mandatory)

AC_CONFIG_MACRO_DIR([m4])

CHECK_FOR_AMIRAMESH

CHECK_FOR_HDF5


# {{{ Handle --enable-assertions
AC_ARG_ENABLE([assertions],
    AS_HELP_STRING([--enable-assertions], [Enable run-time assertions]))

AS_IF([test "x$enable_assertions" = "xyes"],
      ,
      [CPPFLAGS="$CPPFLAGS -DNDEBUG"])
# }}}

# Write config.h file containing settings as preprocessor macros
AC_CONFIG_HEADERS([config.h])

AC_CONFIG_FILES([
  Makefile
  test/Makefile
  psurface.pc
])
AC_OUTPUT

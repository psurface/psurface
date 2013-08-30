AC_INIT([psurface],[1.3.1],[sander@igpm.rwth-aachen.de])
AC_PREREQ([2.60]) # FIXME: what version do we really need?
AC_CONFIG_SRCDIR([.])
AC_CONFIG_AUX_DIR([build-aux])
AM_INIT_AUTOMAKE([-Wall -Werror foreign silent-rules])
AM_PROG_AR   # needed for some non-standard archivers
LT_INIT

AC_PROG_CXX

# {{{ Check for amiramesh
AC_ARG_WITH([amiramesh],
    AS_HELP_STRING([--without-amiramesh], [Do not use libamiramesh]))

# $1: $with_amiramesh
# $2: action if found
# $3: action if not found
AC_DEFUN([CHECK_FOR_AMIRAMESH],[
        AC_LANG_PUSH([C++])
        ac_save_CPPFLAGS="$CPPFLAGS"
        ac_save_LDFLAGS="$LDFLAGS"
        ac_save_LIBS="$LIBS"
        LIBS="$LIBS -lamiramesh"
        AS_IF([test "x$1" != x],[CPPFLAGS="$CPPFLAGS -DHX_HAS_STDIOSTREAM -I$1/include"; LDFLAGS="$LDFLAGS -L$1/lib"])
        AC_CHECK_HEADER([amiramesh/AmiraMesh.h],
                        [AC_LINK_IFELSE([
                           AC_LANG_PROGRAM([[#include "amiramesh/AmiraMesh.h"]],
                                           [[AmiraMesh* am = AmiraMesh::read("test");]]
                                          )],
                           [$2; AMIRA_DIR="$1"],
                           $3
                        )],
                        $3)
        CPPFLAGS="$ac_save_CPPFLAGS"
        LDFLAGS="$ac_save_LDFLAGS"
        LIBS="$ac_save_LIBS"
        AC_LANG_POP([C++])
])

AS_IF([test "x$with_amiramesh" != "xno"],
      [CHECK_FOR_AMIRAMESH(["$with_amiramesh"], [have_amiramesh=yes], [have_amiramesh=no])],
      [have_amiramesh=no])

AS_IF([test "x$have_amiramesh" = "xyes"],
      [CPPFLAGS="$CPPFLAGS -DHAVE_AMIRAMESH -DHX_HAS_STDIOSTREAM -I${AMIRA_DIR}/include";
       LIBS="$LIBS -lamiramesh";
       LDFLAGS="$LDFLAGS -L${AMIRA_DIR}/lib"],
      [AS_IF([test "x$with_amiramesh" = "xyes"],
             [AC_MSG_ERROR([libamiramesh requested but not found])
      ])
])
# }}}


# {{{ Check for hdf5 
AC_ARG_WITH([hdf5],
    AS_HELP_STRING([--without-hdf5], [Do not use libhdf5]))

# $4: $with_hdf5
# $5: action if found
# $6: action if not found
AC_DEFUN([CHECK_FOR_HDF5],[
        AC_LANG_PUSH([C++])
        ac_save_CPPFLAGS="$CPPFLAGS"
        ac_save_LDFLAGS="$LDFLAGS"
        ac_save_LIBS="$LIBS"
        LIBS="$LIBS -lhdf5"
        AS_IF([test "x$with_hdf5" != x],[CPPFLAGS="$CPPFLAGS -I$with_hdf5/include"; LDFLAGS="$LDFLAGS -L$with_hdf5/lib"])
        AC_CHECK_HEADER([hdf5.h],
                        [AC_LINK_IFELSE([
                           AC_LANG_PROGRAM([[#include "hdf5.h"]],
                                           [[" $with_hdf5 $3 $4"]] 
                                          )],
                           [[have_hdf5=yes]; HDF5_DIR=""$with_hdf5""],
                           [have_hdf5=no] 
                        )],
                        [have_hdf5=no])
        CPPFLAGS="$ac_save_CPPFLAGS"
        LDFLAGS="$ac_save_LDFLAGS"
        LIBS="$ac_save_LIBS"
        AC_LANG_POP([C++])
])

AS_IF([test "x$with_hdf5" != "xno"],
      [CHECK_FOR_HDF5(["$with_hdf5"], [have_hdf5=yes], [have_hdf5=no])],
      [have_hdf5=no])

AS_IF([test "x$have_hdf5" = "xyes"],
      [CPPFLAGS="$CPPFLAGS -I${HDF5_DIR}/include";
       LIBS="$LIBS -lhdf5";
       LDFLAGS="$LDFLAGS -L${HDF5_DIR}/lib"],
      [AS_IF([test "x$with_hdf5" = "xyes"],
             [AC_MSG_ERROR([libhdf5 requested but not found])
      ])
])
# }}}

# {{{ Handle --enable-assertions
AC_ARG_ENABLE([assertions],
    AS_HELP_STRING([--enable-assertions], [Enable run-time assertions]))

AS_IF([test "x$enable_assertions" = "xyes"],
      ,
      [CPPFLAGS="$CPPFLAGS -DNDEBUG"])
# }}}

# FIXME: Settings are currently passed via CPPFLAGS
#AC_CONFIG_HEADERS([config.h])

AC_CONFIG_FILES([
  Makefile
])
AC_OUTPUT

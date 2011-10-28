AC_INIT([psurface],[1.3.0],[sander@mi.fu-berlin.de])
AC_PREREQ([2.60]) # FIXME: what version do we really need?
AC_CONFIG_SRCDIR([src])
AC_CONFIG_AUX_DIR([build-aux])
AM_INIT_AUTOMAKE([-Wall -Werror foreign silent-rules])
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
  src/Makefile
])
AC_OUTPUT

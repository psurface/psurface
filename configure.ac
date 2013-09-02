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
  AC_ARG_WITH(hdf5,
    [AC_HELP_STRING([--with-hdf5=PATH],[directory with HDF5 inside])],
    [# expand tilde / other stuff
    eval with_hdf5=$with_hdf5
    ]
  )
  AC_ARG_WITH(hdf5_libs,
    [AC_HELP_STRING([--with-hdf5-libs=LIBS],[additional libraries needed to link hdf5 programs. Those might be needed if your hdf5 library is static. Possible values are: -lz or -lz -lsz.])],[])
 # store values
 ac_save_CPPFLAGS="$CPPFLAGS"
 ac_save_LDFLAGS="$LDFLAGS"
 ac_save_LIBS="$LIBS"

 # start building variables

 # use special HDF5-lib-path if it's set
 if test x$with_hdf5 != x ; then
   #  extract absolute path
   if test -d $with_hdf5; then
     eval with_hdf5=`cd $with_hdf5 ; pwd`
   else
     AC_MSG_ERROR([HDF5-directory $with_hdf5 does not exist])
   fi
   _dune_hdf5_libpath="-L$with_hdf5/lib"
   _dune_hdf5_incpath="-I$with_hdf5/include"
 else
   _dune_hdf5_libpath=""
   _dune_hdf5_incpath=""
 fi

 CPPFLAGS="$CPPFLAGS $_dune_hdf5_incpath"

 direct_HDF5_CPPFLAGS="$_dune_hdf5_incpath"
 nodep_HDF5_CPPFLAGS="$_dune_hdf5_incpath"
 HDF5_CPPFLAGS="$_dune_hdf5_incpath"
 direct_HDF5_LDFLAGS=""
 nodep_HDF5_LDFLAGS=""
 HDF5_LDFLAGS=""
 direct_HDF5_LIBS=""
 nodep_HDF5_LIBS=""
 HDF5_LIBS=""
 HDF5_PARALLEL=0
 # test if we are parallel
 AC_CHECK_DECL(H5_HAVE_PARALLEL, [dnl
	CPPFLAGS="$CPPFLAGS $DUNEMPICPPFLAGS"
	direct_HDF5_CPPFLAGS="$HDF5_CPPFLAGS $DUNEMPICPPFLAGS"
	HDF5_CPPFLAGS="$HDF5_CPPFLAGS \${DUNEMPICPPFLAGS}"
	LDFLAGS="$LDFLAGS $DUNEMPILDFLAGS"
	direct_HDF5_LDFLAGS="$HDF5_LDFLAGS $DUNEMPILDFLAGS"
	HDF5_LDFLAGS="$HDF5_LDFLAGS \${DUNEMPILDFLAGS}"
	LIBS="$DUNEMPILIBS $LIBS"
	direct_HDF5_LIBS="$DUNEMPILIBS $HDF5_LIBS"
	HDF5_LIBS="\${DUNEMPILIBS} $HDF5_LIBS"
	HDF5_PARALLEL=1],[],[#include"H5pubconf.h"])

 # test for an arbitrary header
 AC_CHECK_HEADER([hdf5.h], 
   [HAVE_HDF5=1],
   [HAVE_HDF5=0])
	
 # Just for the configure check.  In the end, -L has to go into LIBS.
 LDFLAGS="$LDFLAGS $_dune_hdf5_libpath"
 # test for lib
 if test x$HAVE_HDF5 = x1 ; then
   AC_CHECK_LIB([hdf5], [H5open],
     [
       direct_HDF5_LIBS="$_dune_hdf5_libpath -lhdf5 $with_hdf5_libs $direct_HDF5_LIBS"
       nodep_HDF5_LIBS="$_dune_hdf5_libpath -lhdf5 $with_hdf5_libs $nodep_HDF5_LIBS"
       HDF5_LIBS="$_dune_hdf5_libpath -lhdf5 $with_hdf5_libs $HDF5_LIBS"
     ],
     [HAVE_HDF5=0], ["$with_hdf5_libs"])
 fi

 # pre-set variable for summary
 with_hdf5="no"

 # did we succeed?
 if test x$HAVE_HDF5 = x1 ; then
   AC_DEFINE(HAVE_HDF5, 1, [Define to 1 if hdf5 was found])

   # proudly show in summary
   with_hdf5="yes"
   if test x"$HDF5_PARALLEL" = x1; then
     with_hdf5="$with_hdf5 (parallel)"
   fi
 else
   # clear variables
   direct_HDF5_CPPFLAGS=
   nodep_HDF5_CPPFLAGS=
   HDF5_CPPFLAGS=
   direct_HDF5_LDFLAGS=
   nodep_HDF5_LDFLAGS=
   HDF5_LDFLAGS=
   direct_HDF5_LIBS=
   nodep_HDF5_LIBS=
   HDF5_LIBS=
   HDF5_PARALLEL=0
 fi

 AC_SUBST([HDF5_CPPFLAGS])
 AC_SUBST([HDF5_LDFLAGS])
 AC_SUBST([HDF5_LIBS])

 # reset values					    
 LIBS="$ac_save_LIBS"
 LDFLAGS="$ac_save_LDFLAGS"
 CPPFLAGS="$ac_save_CPPFLAGS"

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

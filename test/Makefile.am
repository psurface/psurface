# -*- tab-width: 4; indent-tabs-mode: nil -*-
# $Id$

# Magic variable: all programs in TESTS are run when 'make check' is called.
TESTS = gmshiotest \
        simplifytest

# programs just to build when "make check" is used
check_PROGRAMS = $(TESTS)

# define the programs (in alphabetical order)
AM_CPPFLAGS= -I$(top_srcdir)/include/psurface -DPSURFACE_STANDALONE
gmshiotest_SOURCES = gmshiotest.cpp
gmshiotest_CPPFLAGS = $(AM_CPPFLAGS)
gmshiotest_LDADD = $(top_builddir)/libpsurface.la
gmshiotest_LDFLAGS = $(AM_LDFLAGS)

simplifytest_SOURCES = simplifytest.cpp
simplifytest_CPPFLAGS = $(AM_CPPFLAGS)
simplifytest_LDADD = $(top_builddir)/libpsurface.la
simplifytest_LDFLAGS = $(AM_LDFLAGS)

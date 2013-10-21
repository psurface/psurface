SUBDIRS = test

ACLOCAL_AMFLAGS = -I m4

bin_PROGRAMS =  psurface-convert

include_psurface_HEADERS = \
	$(top_srcdir)/AmiraMeshIO.h \
	$(top_srcdir)/Box.h \
	$(top_srcdir)/CircularPatch.h \
	$(top_srcdir)/ContactMapping.h \
	$(top_srcdir)/DirectionFunction.h \
	$(top_srcdir)/DomainPolygon.h \
	$(top_srcdir)/Domains.h \
	$(top_srcdir)/EdgeIntersectionFunctor.h \
	$(top_srcdir)/GlobalNodeIdx.h \
	$(top_srcdir)/HxParamToolBox.h \
	$(top_srcdir)/IntersectionPrimitiveCollector.h \
	$(top_srcdir)/IntersectionPrimitive.h \
	$(top_srcdir)/MultiDimOctree.h \
	$(top_srcdir)/NodeBundle.h \
	$(top_srcdir)/Node.h \
	$(top_srcdir)/NormalProjector.h \
	$(top_srcdir)/PathVertex.h \
	$(top_srcdir)/PlaneParam.h \
	$(top_srcdir)/PointIntersectionFunctor.h \
	$(top_srcdir)/psurfaceAPI.h \
	$(top_srcdir)/PSurfaceFactory.h \
	$(top_srcdir)/PSurface.h \
	$(top_srcdir)/PSurfaceSmoother.h \
	$(top_srcdir)/QualityRequest.h \
	$(top_srcdir)/SparseMatrix.h \
	$(top_srcdir)/StaticMatrix.h \
	$(top_srcdir)/StaticVector.h \
	$(top_srcdir)/SurfaceBase.h \
	$(top_srcdir)/SurfaceParts.h \
        $(top_srcdir)/TargetSurface.h \
	$(top_srcdir)/Triangulator.h \
	$(top_srcdir)/VertexHeap.h \
	$(top_srcdir)/Hdf5IO.h \
	$(top_srcdir)/VtkIO.h \
	$(top_srcdir)/GmshIO.h \
	$(top_srcdir)/streams.hh \
	$(top_srcdir)/indent.hh \
	$(top_srcdir)/common.hh \
	$(top_srcdir)/vtuwriter.hh \
	$(top_srcdir)/dataarraywriter.hh \
	$(top_srcdir)/b64enc.hh
include_psurfacedir = $(includedir)/psurface

lib_LTLIBRARIES= libpsurface.la

libpsurface_la_SOURCES= \
	AmiraMeshIO.cpp \
	CircularPatch.cpp \
	ContactMapping.cpp \
	DomainPolygon.cpp \
	Domains.cpp \
	HxParamToolBox.cpp \
	IntersectionPrimitiveCollector.cpp \
	Iterators.cpp \
	NormalProjector.cpp \
	PlaneParam.cpp \
	PSurface.cpp \
	PSurfaceFactory.cpp \
	PSurfaceSmoother.cpp \
	SurfaceBase.cpp \
	TargetSurface.cpp \
	Triangulator.cpp \
	VtkIO.cpp \
	GmshIO.cpp

if HAVE_HDF5
libpsurface_la_SOURCES += Hdf5IO.cpp
endif

AM_CPPFLAGS= -I$(top_srcdir)/include/psurface -DPSURFACE_STANDALONE
libpsurface_la_CPPFLAGS = $(AM_CPPFLAGS) $(HDF5_CPPFLAGS) $(AMIRAMESH_CPPFLAGS)
libpsurface_la_LIBADD = $(HDF5_LIBS) $(HDF5_LDFLAGS) $(AMIRAMESH_LIBS)
libpsurface_la_LDFLAGS = $(AM_LDFLAGS) $(HDF5_LIBS) $(HDF5_LDFLAGS) $(AMIRAMESH_LDFLAGS)

psurface_convert_SOURCES =  psurface-convert.cpp
psurface_convert_CPPFLAGS = $(AM_CPPFLAGS) $(HDF5_CPPFLAGS) $(AMIRAMESH_CPPFLAGS)
psurface_convert_LDADD = $(HDF5_LIBS) $(HDF5_LDFLAGS)  $(AMIRAMESH_LIBS) libpsurface.la
psurface_convert_LDFLAGS = $(AM_LDFLAGS) $(HDF5_LIBS) $(HDF5_LDFLAGS)  $(AMIRAMESH_LDFLAGS)

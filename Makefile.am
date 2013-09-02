noinst_PROGRAMS =  psurface_convert
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
	$(top_srcdir)/SurfacePath.h \
	$(top_srcdir)/SurfacePathSet.h \
        $(top_srcdir)/TargetSurface.h \
	$(top_srcdir)/Triangulator.h \
	$(top_srcdir)/VertexHeap.h \
	$(top_srcdir)/hdf5IO.h \
	$(top_srcdir)/hdf5IO_simple.h
include_psurfacedir = $(includedir)/psurface

lib_LTLIBRARIES= libpsurface.la

libpsurface_la_SOURCES= \
	AmiraMeshIO.cpp \
	CircularPatch.cpp \
	ContactMapping.cpp \
	DomainPolygon.cpp \
	Domains.cpp \libpsurface.lalibpsurface.la
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
	Triangulator.cpp
AM_CPPFLAGS= -I$(top_srcdir)/include/psurface -DPSURFACE_STANDALONE
libpsurface_la_CPPFLAGS = $(AM_CPPFLAGS) $(HDF5_CPPFLAGS) 
libpsurface_la_LIBADD = $(HDF5_LIBS) $(HDF5_LDFLAGS)
libpsurface_la_LDFLAGS = $(AM_LDFLAGS) $(HDF5_LIBS) $(HDF5_LDFLAGS)
psurface_convert_SOURCES =  AmiraMeshIO.h \
	Iterators.cpp \
	PlaneParam.cpp \
	PlaneParam.h \
	streams.hh \
	indent.hh \
	common.hh \
	vtuwriter.hh \
	dataarraywriter.hh \
	PSurface.h \
	SurfacePathSet.h \
	Domains.h \
	SurfaceBase.h \
	psurfaceAPI.h \
	TargetSurface.h \
	AmiraMeshIO.h \
	SurfaceParts.h \
	psurface_convert.h \
	PSurfaceFactory.h \
	StaticVector.h \
	Domains.cpp \
	AmiraMeshIO.cpp \
	SurfaceBase.cpp \
	TargetSurface.cpp \
	PSurface.cpp \
	PSurfaceFactory.cpp \
	Node.h \
	StaticMatrix.h \
	SparseMatrix.h \
	psurface_convert.cpp\
	psurface_convert.h
psurface_convert_CPPFLAGS = $(AM_CPPFLAGS) $(HDF5_CPPFLAGS) 
psurface_convert_LDADD = $(HDF5_LIBS) $(HDF5_LDFLAGS)
psurface_convert_LDFLAGS = $(AM_LDFLAGS) $(HDF5_LIBS) $(HDF5_LDFLAGS)

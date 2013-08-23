bin_PROGRAMS =  hdf5IO\
		hdf5IOS\
		VTKIO\
		gmshIO\
		psurface_convert
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
hdf5IO_SOURCES = hdf5test.cpp
hdf5IO_CPPFLAGS = $(AM_CPPFLAGS) \
	$(DUNEMPICPPFLAGS) \
	$(AMIRAMESH_CPPFLAGS) \
	-I/usr/local/hdf5/include/\
	-I{top_srcdir}/psurface
hdf5IO_LDADD = \
	$(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) \
	/usr/local/hdf5/lib/libhdf5_cpp.a \
	/usr/local/hdf5/lib/libhdf5_cpp.la \
	/usr/local/hdf5/lib/libhdf5_hl_cpp.a \
	/usr/local/hdf5/lib/libhdf5_hl_cpp.la \
	/usr/local/hdf5/lib/libhdf5_cpp.so.8.0.0 \
	/usr/local/hdf5/lib/libhdf5_hl_cpp.so.8.0.0 \
	/usr/local/lib/libpsurface.la
hdf5IO_LDFLAGS = $(AM_LDFLAGS) \
	$(AMIRAMESH_LDFLAGS) \
	-L/usr/local/hdf5/lib \
	-L/usr/local/lib/ 
hdf5IOS_SOURCES = hdf5test_simple.cpp
hdf5IOS_CPPFLAGS = $(AM_CPPFLAGS) \
	$(DUNEMPICPPFLAGS) \
	$(AMIRAMESH_CPPFLAGS) \
	-I/usr/local/hdf5/include/\
	-I{top_srcdir}/psurface
hdf5IOS_LDADD = \
	$(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) \
	/usr/local/hdf5/lib/libhdf5_cpp.a \
	/usr/local/hdf5/lib/libhdf5_cpp.la \
	/usr/local/hdf5/lib/libhdf5_hl_cpp.a \
	/usr/local/hdf5/lib/libhdf5_hl_cpp.la \
	/usr/local/hdf5/lib/libhdf5_cpp.so.8.0.0 \
	/usr/local/hdf5/lib/libhdf5_hl_cpp.so.8.0.0 \
	/usr/local/lib/libpsurface.la
hdf5IOS_LDFLAGS = $(AM_LDFLAGS) \
	$(AMIRAMESH_LDFLAGS) \
	-L/usr/local/hdf5/lib \
	-L/usr/local/lib/
VTKIO_SOURCES = vtuIO.h	\
		vtuIOtest.cpp \
		b64enc.hh	\
		common.hh	\
		dataarraywriter.hh \
		indent.hh	\
		streams.hh \
		vtuwriter.hh 
VTKIO_CPPFLAGS = $(AM_CPPFLAGS) \
	$(DUNEMPICPPFLAGS) \
	$(AMIRAMESH_CPPFLAGS) \
	-I{top_srcdir}/psurface
VTKIO_LDADD = \
	$(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) \
	/usr/local/lib/libpsurface.la
VTKIO_LDFLAGS = $(AM_LDFLAGS) \
	$(AMIRAMESH_LDFLAGS) \
	-L/usr/local/hdf5/lib \
	-L/usr/local/lib/
gmshIO_SOURCES = gmshIO.h	\
		gm.cpp 
gmshIO_CPPFLAGS = $(AM_CPPFLAGS) \
	$(DUNEMPICPPFLAGS) \
	$(AMIRAMESH_CPPFLAGS) \
	-I/usr/local/hdf5/include/\
	-I{top_srcdir}/psurface
gmshIO_LDADD = \
	$(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) \
	/usr/local/hdf5/lib/libhdf5_cpp.a \
	/usr/local/hdf5/lib/libhdf5_cpp.la \
	/usr/local/hdf5/lib/libhdf5_hl_cpp.a \
	/usr/local/hdf5/lib/libhdf5_hl_cpp.la \
	/usr/local/hdf5/lib/libhdf5_cpp.so.8.0.0 \
	/usr/local/hdf5/lib/libhdf5_hl_cpp.so.8.0.0 \
	/usr/local/lib/libpsurface.la	
gmshIO_LDFLAGS = $(AM_LDFLAGS) \
	$(AMIRAMESH_LDFLAGS) \
	-L/usr/local/hdf5/lib \
	-L/usr/local/lib/
psurface_convert_SOURCES = psurface_convert.cpp\
			   psurface_convert.h
psurface_convert_CPPFLAGS = $(AM_CPPFLAGS) \
	$(DUNEMPICPPFLAGS) \
	$(AMIRAMESH_CPPFLAGS) \
	-I/usr/local/hdf5/include/\
	-I{top_srcdir}/psurface
psurface_convert_LDADD = \
	$(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) \
	/usr/local/hdf5/lib/libhdf5_cpp.a \
	/usr/local/hdf5/lib/libhdf5_cpp.la \
	/usr/local/hdf5/lib/libhdf5_hl_cpp.a \
	/usr/local/hdf5/lib/libhdf5_hl_cpp.la \
	/usr/local/hdf5/lib/libhdf5_cpp.so.8.0.0 \
	/usr/local/hdf5/lib/libhdf5_hl_cpp.so.8.0.0 \
	/usr/local/lib/libpsurface.la	
psurface_convert_LDFLAGS = $(AM_LDFLAGS) \
	$(AMIRAMESH_LDFLAGS) \
	-L/usr/local/hdf5/lib \
	-L/usr/local/lib/


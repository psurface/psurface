#!/bin/bash

# Clean directory structure
rm -rf libpsurface
mkdir libpsurface
mkdir libpsurface/include
mkdir libpsurface/include/psurface
mkdir libpsurface/include/contact
mkdir libpsurface/src
mkdir libpsurface/lib

# Copy the licence file
cp COPYING                  libpsurface

cp TargetSurface.h.standalone                         libpsurface/include/TargetSurface.h

cp psurface.h               libpsurface/include/psurface
cp psurface.h               libpsurface/include/
cp CircularPatch.h          libpsurface/include/psurface
cp DomainPolygon.h          libpsurface/include/psurface
cp Domains.h                libpsurface/include/psurface
cp GlobalNodeIdx.h          libpsurface/include/psurface
cp McPointerSurfaceParts.h  libpsurface/include/psurface
cp McSurfaceBase.h          libpsurface/include/psurface
cp Node.h                   libpsurface/include/psurface
cp NodeBundle.h             libpsurface/include/psurface
cp Parametrization.h        libpsurface/include/psurface
cp PlaneParam.h             libpsurface/include/psurface
cp SurfacePath.h            libpsurface/include/psurface
cp SurfacePathSet.h         libpsurface/include/psurface
cp Box.h                    libpsurface/include/psurface
cp MultiDimOctree.h         libpsurface/include/psurface
cp PointIntersectionFunctor.h  libpsurface/include/psurface
cp StaticVector.h           libpsurface/include/psurface
cp StaticMatrix.h           libpsurface/include/psurface
cp SparseMatrix.h           libpsurface/include/psurface


cp psurface.cpp                     libpsurface/src
cp CircularPatch.cpp                libpsurface/src
cp DomainPolygon.cpp                libpsurface/src
cp Domains.cpp                      libpsurface/src
cp Iterators.cpp                    libpsurface/src
cp Parametrization.cpp              libpsurface/src
cp PlaneParam.cpp                   libpsurface/src
cp TargetSurface.cpp.standalone     libpsurface/src/TargetSurface.cpp
cp debugCode.cpp                    libpsurface/src

cp contact.h                    libpsurface/include/psurface
cp ContactToolBox.h             libpsurface/include/psurface
cp IntersectionPrimitive.h      libpsurface/include/psurface
cp ContactBoundary.h            libpsurface/include/psurface
cp NormalProjector.h            libpsurface/include/psurface
cp ContactMapping.h             libpsurface/include/psurface

cp contact.cpp                  libpsurface/src
cp NormalProjector.cpp          libpsurface/src
cp buildContactSurface.cpp      libpsurface/src
cp extractMergedGrid.cpp        libpsurface/src
cp ContactMapping.cpp           libpsurface/src

cp Makefile.standalone      libpsurface/src/Makefile

tar zvcf libpsurface.tar.gz libpsurface

#!/bin/bash

# Clean directory structure
rm -rf libpsurface
mkdir libpsurface
mkdir libpsurface/include
mkdir libpsurface/include/mclib
mkdir libpsurface/include/psurface
mkdir libpsurface/include/contact
mkdir libpsurface/include/hxfield
mkdir libpsurface/include/hxsurface
mkdir libpsurface/src
mkdir libpsurface/lib

cp local_mclib/include/hxfield/HxFieldWinDLLApi.h     libpsurface/include/hxfield
cp local_mclib/include/hxfield/oint.h                 libpsurface/include/hxfield

cp local_mclib/include/hxsurface/HxSurfaceWinDLLApi.h libpsurface/include/hxsurface
cp TargetSurface.h.standalone                         libpsurface/include/TargetSurface.h

cp local_mclib/include/mclib/McBox2f.h                libpsurface/include/mclib
cp local_mclib/include/mclib/McBox3f.h                libpsurface/include/mclib
cp local_mclib/include/mclib/McDVector.h              libpsurface/include/mclib
cp local_mclib/include/mclib/McFHeap.h                libpsurface/include/mclib
cp local_mclib/include/mclib/McHandable.h             libpsurface/include/mclib
cp local_mclib/include/mclib/McHashTable.h            libpsurface/include/mclib
cp local_mclib/include/mclib/McMat3f.h                libpsurface/include/mclib
cp local_mclib/include/mclib/McMat4f.h                libpsurface/include/mclib
cp local_mclib/include/mclib/McPrimType.h             libpsurface/include/mclib
cp local_mclib/include/mclib/McProgressInterface.h    libpsurface/include/mclib
cp local_mclib/include/mclib/McRotation.h             libpsurface/include/mclib
cp local_mclib/include/mclib/McSArray.h               libpsurface/include/mclib
cp local_mclib/include/mclib/McSparseMatrix.h         libpsurface/include/mclib
cp local_mclib/include/mclib/McStdlib.h               libpsurface/include/mclib
cp local_mclib/include/mclib/McStdio.h                libpsurface/include/mclib
cp local_mclib/include/mclib/McVec2f.h                libpsurface/include/mclib
cp local_mclib/include/mclib/McVec2i.h                libpsurface/include/mclib
cp local_mclib/include/mclib/McVec3f.h                libpsurface/include/mclib
cp local_mclib/include/mclib/McVec3f_impl.h           libpsurface/include/mclib
cp local_mclib/include/mclib/McVec3d.h                libpsurface/include/mclib
cp local_mclib/include/mclib/McVec3d_impl.h           libpsurface/include/mclib
cp local_mclib/include/mclib/McVec3i.h                libpsurface/include/mclib
cp local_mclib/include/mclib/McVec4f.h                libpsurface/include/mclib
cp local_mclib/include/mclib/McVec4i.h                libpsurface/include/mclib
cp local_mclib/include/mclib/McWildMatch.h            libpsurface/include/mclib
cp local_mclib/include/mclib/McWinDLLApi.h            libpsurface/include/mclib

cp local_mclib/McDArray.h                 libpsurface/include/mclib
cp local_mclib/McSmallArray.h             libpsurface/include/mclib
cp local_mclib/McOctree.h                 libpsurface/include/mclib
cp local_mclib/McVec2d.h                  libpsurface/include/mclib
cp local_mclib/McMat3d.h                  libpsurface/include/mclib
cp MyMcVec3f.h                            libpsurface/include/mclib

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
cp psurfaceAPI.h            libpsurface/include/psurface

cp local_mclib/McMat3f.cpp                libpsurface/src
cp local_mclib/McHashTable.cpp            libpsurface/src


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

cp contact.cpp                  libpsurface/src
cp local_mclib/McMat3d.cpp      libpsurface/src
cp NormalProjector.cpp          libpsurface/src
cp buildContactSurface.cpp      libpsurface/src
cp extractMergedGrid.cpp        libpsurface/src

cp Makefile.standalone      libpsurface/src/Makefile

tar zvcf libpsurface.tar.gz libpsurface

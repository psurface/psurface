#!/bin/bash

AMIRA_ROOT=/home/haile/sander/amira-ZIB-2007-1

# Clean directory structure
rm -rf libpsurface
mkdir libpsurface
mkdir libpsurface/include
mkdir libpsurface/include/mclib
mkdir libpsurface/include/parametrization
mkdir libpsurface/include/contact
mkdir libpsurface/include/hxfield
mkdir libpsurface/include/hxsurface
mkdir libpsurface/src
mkdir libpsurface/lib

cp $AMIRA_ROOT/include/hxfield/HxFieldWinDLLApi.h     libpsurface/include/hxfield
cp $AMIRA_ROOT/include/hxfield/oint.h                 libpsurface/include/hxfield

cp $AMIRA_ROOT/include/hxsurface/HxSurfaceWinDLLApi.h libpsurface/include/hxsurface
cp TargetSurface.h.standalone                         libpsurface/include/TargetSurface.h

cp $AMIRA_ROOT/include/mclib/McBox2f.h                libpsurface/include/mclib
cp $AMIRA_ROOT/include/mclib/McBox3f.h                libpsurface/include/mclib
cp $AMIRA_ROOT/include/mclib/McDVector.h              libpsurface/include/mclib
cp $AMIRA_ROOT/include/mclib/McFHeap.h                libpsurface/include/mclib
cp $AMIRA_ROOT/include/mclib/McHandable.h             libpsurface/include/mclib
cp $AMIRA_ROOT/include/mclib/McHashTable.h            libpsurface/include/mclib
cp $AMIRA_ROOT/include/mclib/McMat3f.h                libpsurface/include/mclib
cp $AMIRA_ROOT/include/mclib/McMat4f.h                libpsurface/include/mclib
cp $AMIRA_ROOT/include/mclib/McPrimType.h             libpsurface/include/mclib
cp $AMIRA_ROOT/include/mclib/McProgressInterface.h    libpsurface/include/mclib
cp $AMIRA_ROOT/include/mclib/McRotation.h             libpsurface/include/mclib
cp $AMIRA_ROOT/include/mclib/McSArray.h               libpsurface/include/mclib
cp $AMIRA_ROOT/include/mclib/McSparseMatrix.h         libpsurface/include/mclib
cp $AMIRA_ROOT/include/mclib/McStdlib.h               libpsurface/include/mclib
cp $AMIRA_ROOT/include/mclib/McVec2f.h                libpsurface/include/mclib
cp $AMIRA_ROOT/include/mclib/McVec2i.h                libpsurface/include/mclib
cp $AMIRA_ROOT/include/mclib/McVec3f.h                libpsurface/include/mclib
cp $AMIRA_ROOT/include/mclib/McVec3f_impl.h           libpsurface/include/mclib
cp $AMIRA_ROOT/include/mclib/McVec3d.h                libpsurface/include/mclib
cp $AMIRA_ROOT/include/mclib/McVec3d_impl.h           libpsurface/include/mclib
cp $AMIRA_ROOT/include/mclib/McVec3i.h                libpsurface/include/mclib
cp $AMIRA_ROOT/include/mclib/McVec4f.h                libpsurface/include/mclib
cp $AMIRA_ROOT/include/mclib/McVec4i.h                libpsurface/include/mclib
cp $AMIRA_ROOT/include/mclib/McWildMatch.h            libpsurface/include/mclib
cp $AMIRA_ROOT/include/mclib/McWinDLLApi.h            libpsurface/include/mclib

cp local_mclib/McException.h              libpsurface/include/mclib
cp local_mclib/McDArray.h                 libpsurface/include/mclib
cp local_mclib/McSmallArray.h             libpsurface/include/mclib
cp local_mclib/McOctree.h                 libpsurface/include/mclib
cp local_mclib/McVec2d.h                  libpsurface/include/mclib
cp local_mclib/McMat3d.h                  libpsurface/include/mclib
cp local_mclib/MyMcVec3f.h                libpsurface/include/mclib

cp psurface.h               libpsurface/include/parametrization
cp psurface.h               libpsurface/include/
cp CircularPatch.h          libpsurface/include/parametrization
cp DomainPolygon.h          libpsurface/include/parametrization
cp Domains.h                libpsurface/include/parametrization
cp GlobalNodeIdx.h          libpsurface/include/parametrization
cp McPointerSurfaceParts.h  libpsurface/include/parametrization
cp McSurfaceBase.h          libpsurface/include/parametrization
cp Node.h                   libpsurface/include/parametrization
cp NodeBundle.h             libpsurface/include/parametrization
cp Parametrization.h        libpsurface/include/parametrization
cp PlaneParam.h             libpsurface/include/parametrization
cp SurfacePath.h            libpsurface/include/parametrization
cp SurfacePathSet.h         libpsurface/include/parametrization
cp parametrizationAPI.h     libpsurface/include/parametrization

cp local_mclib/McMat3f.cpp                libpsurface/src
cp local_mclib/McHashTable.cpp            libpsurface/src
cp local_mclib/McException.cpp            libpsurface/src


cp psurface.cpp                     libpsurface/src
cp CircularPatch.cpp                libpsurface/src
cp DomainPolygon.cpp                libpsurface/src
cp Domains.cpp                      libpsurface/src
cp Iterators.cpp                    libpsurface/src
cp Parametrization.cpp              libpsurface/src
cp PlaneParam.cpp                   libpsurface/src
cp TargetSurface.cpp.standalone     libpsurface/src/TargetSurface.cpp
cp debugCode.cpp                    libpsurface/src

cp contact.h                    libpsurface/include/contact
cp ContactToolBox.h             libpsurface/include/contact
cp IntersectionPrimitive.h      libpsurface/include/contact
cp ContactBoundary.h            libpsurface/include/contact
cp NormalProjector.h            libpsurface/include/contact

cp contact.cpp                  libpsurface/src
cp local_mclib/McMat3d.cpp      libpsurface/src
cp NormalProjector.cpp          libpsurface/src
cp buildContactSurface.cpp      libpsurface/src
cp extractMergedGrid.cpp        libpsurface/src

cp Makefile.standalone      libpsurface/src/Makefile

tar zvcf libpsurface.tar.gz libpsurface

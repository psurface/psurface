#!/bin/bash

BASENAME=libpsurface-1.1.0

# Clean directory structure
rm -rf $BASENAME
mkdir $BASENAME
mkdir $BASENAME/include
mkdir $BASENAME/include/psurface
mkdir $BASENAME/include/contact
mkdir $BASENAME/src
mkdir $BASENAME/lib

# Copy the licence file
cp COPYING                      $BASENAME

cp TargetSurface.h.standalone   $BASENAME/include/TargetSurface.h

cp psurface.h                   $BASENAME/include/psurface
cp psurface.h                   $BASENAME/include/
cp CircularPatch.h              $BASENAME/include/psurface
cp DomainPolygon.h              $BASENAME/include/psurface
cp Domains.h                    $BASENAME/include/psurface
cp GlobalNodeIdx.h              $BASENAME/include/psurface
cp McPointerSurfaceParts.h      $BASENAME/include/psurface
cp SurfaceBase.h                $BASENAME/include/psurface
cp Node.h                       $BASENAME/include/psurface
cp NodeBundle.h                 $BASENAME/include/psurface
cp PSurface.h                   $BASENAME/include/psurface
cp PlaneParam.h                 $BASENAME/include/psurface
cp SurfacePath.h                $BASENAME/include/psurface
cp SurfacePathSet.h             $BASENAME/include/psurface
cp Box.h                        $BASENAME/include/psurface
cp MultiDimOctree.h             $BASENAME/include/psurface
cp PointIntersectionFunctor.h   $BASENAME/include/psurface
cp StaticVector.h               $BASENAME/include/psurface
cp StaticMatrix.h               $BASENAME/include/psurface
cp SparseMatrix.h               $BASENAME/include/psurface


cp psurface.cpp                 $BASENAME/src
cp CircularPatch.cpp            $BASENAME/src
cp DomainPolygon.cpp            $BASENAME/src
cp Domains.cpp                  $BASENAME/src
cp Iterators.cpp                $BASENAME/src
cp PSurface.cpp                 $BASENAME/src
cp PlaneParam.cpp               $BASENAME/src
cp TargetSurface.cpp.standalone $BASENAME/src/TargetSurface.cpp
cp debugCode.cpp                $BASENAME/src

cp contact.h                    $BASENAME/include/psurface
cp ContactToolBox.h             $BASENAME/include/psurface
cp IntersectionPrimitive.h      $BASENAME/include/psurface
cp ContactBoundary.h            $BASENAME/include/psurface
cp NormalProjector.h            $BASENAME/include/psurface
cp ContactMapping.h             $BASENAME/include/psurface

cp contact.cpp                  $BASENAME/src
cp NormalProjector.cpp          $BASENAME/src
cp buildContactSurface.cpp      $BASENAME/src
cp extractMergedGrid.cpp        $BASENAME/src
cp ContactMapping.cpp           $BASENAME/src

cp Makefile.standalone          $BASENAME/src/Makefile

tar zvcf $BASENAME.tar.gz $BASENAME

/////////////////////////////////////////////////////////////////
/*
 * $Id: TargetSurface.cpp.standalone,v 1.1 2007/10/17 13:16:55 sander Exp $
 *
 * $Log: TargetSurface.cpp.standalone,v $
 *
 * This new class contains the few methods of Surface that I really need.
 * When this package is compiled with -DPSURFACE_STANDALONE, it is used
 * as the target surface of the parametrization.  When compiled within
 * Amira, the standard Surface class is still used (the real one, no
 * local copy).  That way I can still watch target surfaces with the
 * normal SurfaceView module.
 * mailtoauthor: sander@zib.de
 *
 *
 */
/////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <ctype.h>

#include "TargetSurface.h"

using namespace psurface;

void Surface::getBoundingBox(float bbox[6]) const
{
    if (points.size()) {
        bbox[0] = bbox[1] = points[0][0];
        bbox[2] = bbox[3] = points[0][1];
        bbox[4] = bbox[5] = points[0][2];
    } else {
        for (int i=0; i<6; i++)
            bbox[i] = 0;
    }

    for (int i=0 ; i<points.size() ; i++) {
        for (int k=0 ; k<3 ; k++) {
            if (bbox[k*2]>points[i][k])
                bbox[k*2]=points[i][k];
            if (bbox[k*2+1]<points[i][k])
                bbox[k*2+1]=points[i][k];
        }
    }
}


void Surface::computeTrianglesPerPoint()
{
    int nPoints = points.size();
    int nTriangles = triangles.size();

    trianglesPerPoint.resize(nPoints);

    for (int k=0; k<nPoints; k++)
        trianglesPerPoint[k].resize(0);
    
    for (int i=0; i<nTriangles; i++) {
        Surface::Triangle& tri = triangles[i];
        trianglesPerPoint[tri.points[0]].push_back(i);
        trianglesPerPoint[tri.points[1]].push_back(i);
        trianglesPerPoint[tri.points[2]].push_back(i);
    }
}

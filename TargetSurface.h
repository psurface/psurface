/////////////////////////////////////////////////////////////////
/*
 * $Id: TargetSurface.h.standalone,v 1.1 2007/10/17 13:16:55 sander Exp $
 *
 * $Log: TargetSurface.h.standalone,v $
 *
 * This new class contains the few methods of Surface that I really need.
 * When this package is compiled with -DPSURFACE_STANDALONE, it is used
 * as the target surface of the parametrization.  When compiled within
 * Amira, the standard Surface class is still used (the real one, no
 * local copy).  That way I can still watch target surfaces with the
 * normal SurfaceView module.
 * mailtoauthor: sander@zib.de
 *
 */
/////////////////////////////////////////////////////////////////
#ifndef TARGET_SURFACE_H
#define TARGET_SURFACE_H

#include <stdio.h>

#include <vector>
#include "StaticVector.h"

namespace psurface {

/** This class implements unstructured, triangulated surfaces.
    
    This implementation works with indices instead of pointers, this allows
    easy realloc. The surface contains one big array of @c points, and of
    @c triangles. */

class Surface { 

  public:

    // ==================== Subclass declaration ====================

    /**@name Subclass declaration */ //@{
    
    /// This class represents a surface triangle.
    class Triangle {
      public:
        /** Indices of the three points of the triangle. <b> Note:</b> In the
            file format, the first point in the surface's point array is 1.
            In this data structure, the first point has the index 0.
            Conversion is done when writing and reading. */
        std::array<int,3> points;

    };
    
    //@}

    // ======================= Member variables =====================

    /**@name Member variables */ //@{

    /// Array of all surface points (required).
    std::vector<StaticVector<float,3> > points;

    /// Array of all surface triangles (required).
    std::vector<Triangle> triangles;       

    /** Stores for each point all incident triangles. To initialize this
        array the method @c computeTrianglesPerPoint() has to be called
        explicitely. */
    std::vector< std::vector<int> > trianglesPerPoint;

    //@}

    // ======================== Member methods ======================

    /**@name Member methods */ //@{

    /// Returns bounding box of point coordinates.
    void getBoundingBox(float box[6]) const;

    //@}

    // ====================== Internal methods ======================

    /**@name Internal methods */ //@{

    /** Initializes the array @c trianglePerPoint. */
    void computeTrianglesPerPoint();

    //@}

};

} // namespace psurface

#endif


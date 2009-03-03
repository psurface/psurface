#ifndef CONTACT_TOOLBOX_H
#define CONTACT_TOOLBOX_H

#include <mclib/McBox3f.h>
#include <mclib/McVec3f.h>
#include <mclib/McVec3i.h>

#include "Node.h"
#include "IntersectionPrimitive.h"

#include <vector>

class Parametrization;
class DomainTriangle;
class ContactBoundary;
class Surface;

namespace ContactToolBox {

    void buildContactSurface(Parametrization* cPar,
                             const Surface* surf1,  const Surface* surf2,
                             float epsilon, 
                             void (*obsDirections)(const double* pos, double* dir)
                             );

    void contactOracle(const Surface* surf1, const Surface* surf2,
                       std::vector<int>& contactNodes1, std::vector<int>& contactNodes2,
                       float epsilon);

    McVec3f getClosestPointOnTriangle(const McVec3f& p0,
                                                  const McVec3f& p1,
                                                  const McVec3f& p2,          
                                      const McVec3f& candidate);

    void computeContactPatch(const Surface* surf, ContactBoundary& cBound);

    bool isCompletelyCovered(Parametrization* cPar, int tri, const DomainTriangle* cT);

    /**@name Routines for the extraction of the merged grid */
    //@{

    class IntersectionAlt {
    public:
        McVec2f pos;

        McVec2f localTargetCoords;

    };

    void extractMergedGrid(Parametrization* cPar,
                           std::vector<IntersectionPrimitive>& mergedGrid);

    //@}
};
#endif

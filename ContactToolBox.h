#ifndef CONTACT_TOOLBOX_H
#define CONTACT_TOOLBOX_H

#include <psurface/Box.h>
#include <psurface/StaticVector.h>

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

    StaticVector<float,3> getClosestPointOnTriangle(const StaticVector<float,3>& p0,
                                                  const StaticVector<float,3>& p1,
                                                  const StaticVector<float,3>& p2,          
                                      const StaticVector<float,3>& candidate);

    void computeContactPatch(const Surface* surf, ContactBoundary& cBound);

    bool isCompletelyCovered(Parametrization* cPar, int tri, const DomainTriangle* cT);

    /**@name Routines for the extraction of the merged grid */
    //@{

    class IntersectionAlt {
    public:
        StaticVector<float,2> pos;

        StaticVector<float,2> localTargetCoords;

    };

    void extractMergedGrid(Parametrization* cPar,
                           std::vector<IntersectionPrimitive>& mergedGrid);

    //@}
};
#endif

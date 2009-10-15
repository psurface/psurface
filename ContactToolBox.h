#ifndef CONTACT_TOOLBOX_H
#define CONTACT_TOOLBOX_H

#include <psurface/Box.h>
#include <psurface/StaticVector.h>

#include "Node.h"
#include "IntersectionPrimitive.h"

#include <vector>

template <int dim, class ctype>
class PSurface;
class DomainTriangle;
class ContactBoundary;
class Surface;


namespace ContactToolBox {

    void buildContactSurface(PSurface<2,float>* cPar,
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

    bool isCompletelyCovered(PSurface<2,float>* cPar, int tri, const DomainTriangle* cT);

    /** \brief Extraction of the merged grid */
    void extractMergedGrid(PSurface<2,float>* cPar,
                           std::vector<IntersectionPrimitive<float> >& mergedGrid);

    //@}
};
#endif

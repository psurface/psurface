#ifndef CONTACT_TOOLBOX_H
#define CONTACT_TOOLBOX_H

#include <vector>

#include <psurface/StaticVector.h>

template <int dim, class ctype>
class PSurface;
class ContactBoundary;
class Surface;
template <int dimworld, class ctype>
class DirectionFunction;

template <class ctype>
struct ContactToolBox {

    static void buildContactSurface(PSurface<2,ctype>* cPar,
                             const Surface* surf1,  const Surface* surf2,
                             ctype epsilon, 
                                    const DirectionFunction<3,ctype>* domainDirection,
                                    const DirectionFunction<3,ctype>* targetDirection
                             );

    static void contactOracle(const Surface* surf1, const Surface* surf2,
                       std::vector<int>& contactNodes1, std::vector<int>& contactNodes2,
                       ctype epsilon);

    static StaticVector<ctype,3> getClosestPointOnTriangle(const StaticVector<ctype,3>& p0,
                                                  const StaticVector<ctype,3>& p1,
                                                  const StaticVector<ctype,3>& p2,          
                                      const StaticVector<ctype,3>& candidate);

    static void computeContactPatch(const Surface* surf, ContactBoundary& cBound);

};
#endif

#ifndef CONTACT_TOOLBOX_H
#define CONTACT_TOOLBOX_H

#include <vector>

#include <psurface/StaticVector.h>

template <int dim, class ctype>
class PSurface;
class Surface;
template <int dimworld, class ctype>
class DirectionFunction;

template <class ctype>
struct ContactToolBox {

    static void buildContactSurface(PSurface<2,ctype>* cPar,
                                    const Surface* surf1,  const Surface* surf2,
                                    const DirectionFunction<3,ctype>* domainDirection,
                                    const DirectionFunction<3,ctype>* targetDirection
                                    );

};
#endif

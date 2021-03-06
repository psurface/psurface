#ifndef INTERSECTION_PRIMITIVE_COLLECTOR_H
#define INTERSECTION_PRIMITIVE_COLLECTOR_H

#include <vector>
#include "IntersectionPrimitive.h"

#include "psurfaceAPI.h"

namespace psurface {

template <int dim, class ctype>
class PSurface;

template <class ctype>
class PSURFACE_API IntersectionPrimitiveCollector {

public:

    /** \brief Collection of the merged grid */
    static void collect(PSurface<2,ctype>* psurface,
                        std::vector<IntersectionPrimitive<2,ctype> >& mergedGrid);

    /** \brief Collection of the merged grid */
    static void collect(const PSurface<1,ctype>* psurface,
                        std::vector<IntersectionPrimitive<1,ctype> >& mergedGrid);

};

} // namespace psurface

#endif

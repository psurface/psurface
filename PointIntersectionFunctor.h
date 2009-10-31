#ifndef POINT_INTERSECTION_FUNCTOR_H
#define POINT_INTERSECTION_FUNCTOR_H

#include <psurface/StaticVector.h>

/** \brief Functor class needed to insert StaticVector<.,3> objects into a MultiDimOctree
 */
template <class ctype>
struct PointIntersectionFunctor
{

    bool operator()(const std::tr1::array<ctype,3>& lower, 
                    const std::tr1::array<ctype,3>& upper, const StaticVector<ctype,3>& item) const {
        return (lower[0] <= item[0]) && (item[0] <= upper[0])
            && (lower[1] <= item[1]) && (item[1] <= upper[1])
            && (lower[2] <= item[2]) && (item[2] <= upper[2]);
    }

};


#endif

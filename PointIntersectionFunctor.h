#ifndef POINT_INTERSECTION_FUNCTOR_H
#define POINT_INTERSECTION_FUNCTOR_H

#include <mclib/McVec3f.h>

/** \brief Functor class needed to insert McVec3f objects into a MultiDimOctree
 */
struct PointIntersectionFunctor
{
    bool operator()(const Box<float, 3>& box, const McVec3f* item) const {
        return (box.lower()[0] <= item->x) && (item->x <= box.upper()[0])
            && (box.lower()[1] <= item->y) && (item->y <= box.upper()[1])
            && (box.lower()[2] <= item->z) && (item->z <= box.upper()[2]);
    }

    bool operator()(const std::tr1::array<float,3>& lower, 
                    const std::tr1::array<float,3>& upper, const McVec3f& item) const {
        return (lower[0] <= item.x) && (item.x <= upper[0])
            && (lower[1] <= item.y) && (item.y <= upper[1])
            && (lower[2] <= item.z) && (item.z <= upper[2]);
    }

};


#endif

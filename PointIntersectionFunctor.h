#ifndef POINT_INTERSECTION_FUNCTOR_H
#define POINT_INTERSECTION_FUNCTOR_H

#include <psurface/StaticVector.h>

/** \brief Functor class needed to insert StaticVector<float,3> objects into a MultiDimOctree
 */
struct PointIntersectionFunctor
{
    bool operator()(const Box<float, 3>& box, const StaticVector<float,3>* item) const {
        return (box.lower()[0] <= (*item)[0]) && ((*item)[0] <= box.upper()[0])
            && (box.lower()[1] <= (*item)[1]) && ((*item)[1] <= box.upper()[1])
            && (box.lower()[2] <= (*item)[2]) && ((*item)[2] <= box.upper()[2]);
    }

    bool operator()(const std::tr1::array<float,3>& lower, 
                    const std::tr1::array<float,3>& upper, const StaticVector<float,3>& item) const {
        return (lower[0] <= item[0]) && (item[0] <= upper[0])
            && (lower[1] <= item[1]) && (item[1] <= upper[1])
            && (lower[2] <= item[2]) && (item[2] <= upper[2]);
    }

};


#endif

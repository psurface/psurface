#ifndef PSURFACE_DIRECTION_FUNCTION_H
#define PSURFACE_DIRECTION_FUNCTION_H

#include "StaticVector.h"

namespace psurface {

/** \brief Abstract base class for direction fields on simplicial surfaces */
template <int dimworld, class ctype>
struct DirectionFunction
{
    virtual ~DirectionFunction() {};
};

/** \brief Abstract base class for direction fields that are given by closed-form expressions. 
*/
template <int dimworld, class ctype>
struct AnalyticDirectionFunction
    : public DirectionFunction<dimworld,ctype>
{
    /** \brief Return a direction for a given world position */
    virtual StaticVector<ctype,dimworld> operator()(const StaticVector<ctype,dimworld>& position) const = 0;
};

/** \brief Abstract base class for direction fields that are given per vertex.
*/
template <int dimworld, class ctype>
struct DiscreteDirectionFunction
    : public DirectionFunction<dimworld,ctype>
{
    /** \brief Return a direction for a given vertex index */
    virtual StaticVector<ctype,dimworld> operator()(size_t index) const = 0;
};

} // namespace psurface

#endif

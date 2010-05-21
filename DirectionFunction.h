#ifndef PSURFACE_DIRECTION_FUNCTION_H
#define PSURFACE_DIRECTION_FUNCTION_H

#include <psurface/StaticVector.h>

/** \brief Abstract base class for direction fields on simplicial surfaces */
template <int dimworld, class ctype>
struct DirectionFunction
{
    /** \brief Abstract virtual destructor, to prevent people 
        from creating objects of this class directly.
    */
    virtual ~DirectionFunction() = 0;
};

/** \brief Abstract base class for direction fields that are given by closed-form expressions. 
*/
template <int dimworld, class ctype>
struct AnalyticDirectionFunction
    : public DirectionFunction<dimworld,ctype>
{
    virtual StaticVector<ctype,dimworld> operator()(const StaticVector<ctype,dimworld>& position) const = 0;
};

/** \brief Abstract base class for direction fields that are given per vertex.
*/
template <int dimworld, class ctype>
struct DiscreteDirectionFunction
    : public DirectionFunction<dimworld,ctype>
{
    virtual StaticVector<ctype,dimworld> operator()(const StaticVector<ctype,dimworld>& position) const = 0;
};



#endif

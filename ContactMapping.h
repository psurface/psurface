#ifndef CONTACT_MAPPING_HH
#define CONTACT_MAPPING_HH

#include <vector>
#include <iostream>

// Check for VC9 / VS2008 with installed feature pack.
#if defined(_MSC_VER) && (_MSC_VER>=1500)
    #if defined(_CPPLIB_VER) && _CPPLIB_VER>=505
        #include <array>
    #else
        #error Please install the Visual C++ 2008 Feature Pack for TR1 support.
    #endif
#else
    #include <tr1/array>
#endif

#include "StaticVector.h"
#include "PSurface.h"
#include "IntersectionPrimitive.h"
#include "IntersectionPrimitiveCollector.h"

template <int dimworld, class ctype>
struct DirectionFunction;

template <int dim, class ctype>
class ContactMapping {};

template <class ctype>
class ContactMapping<2,ctype>
{
public:
    void build(const std::vector<std::tr1::array<ctype,2> >& coords1,  ///< The vertices of the first surface as \f$x_0 ,y_0 ,z_0, x_1, y_1, z_1 ...\f$
               const std::vector<std::tr1::array<int,2> >& tri1,       ///< The triangles of the first surface
               const std::vector<std::tr1::array<ctype,2> >& coords2,  ///< The vertices of the second surface
               const std::vector<std::tr1::array<int,2> >& tri2,
               const DirectionFunction<2,ctype>* domainDirection = NULL,
               const DirectionFunction<2,ctype>* targetDirection = NULL
               );

    void getOverlaps(std::vector<IntersectionPrimitive<1,ctype> >& overlaps)
    {
        IntersectionPrimitiveCollector<ctype>::collect(&psurface_, overlaps);
    }

    // /////////////////////////////////////////////

private:

    void computeDiscreteDomainDirections(const DirectionFunction<2,ctype>* direction,
                                         std::vector<StaticVector<ctype,2> >& normals);

    void computeDiscreteTargetDirections(const std::vector<std::tr1::array<int,2> >& elements,
                                         const DirectionFunction<2,ctype>* direction,
                                         std::vector<StaticVector<ctype,2> >& normals);

    PSurface<1,ctype> psurface_;

};


template <class ctype>
class ContactMapping<3,ctype>
{
public:	

    void build(const std::vector<std::tr1::array<ctype,3> >& coords1,  ///< The vertices of the first surface as \f$x_0 ,y_0 ,z_0, x_1, y_1, z_1 ...\f$
               const std::vector<std::tr1::array<int,3> >& tri1,       ///< The triangles of the first surface
               const std::vector<std::tr1::array<ctype,3> >& coords2,  ///< The vertices of the second surface
               const std::vector<std::tr1::array<int,3> >& tri2,       ///< The triangles of the second surface
               const DirectionFunction<3,ctype>* domainDirection = NULL,
               const DirectionFunction<3,ctype>* targetDirection = NULL
               );

    void getOverlaps(std::vector<IntersectionPrimitive<2,ctype> >& overlaps) {
        IntersectionPrimitiveCollector<ctype>::collect(&psurface_, overlaps);
    }

private:

    PSurface<2,ctype> psurface_;

};

#endif

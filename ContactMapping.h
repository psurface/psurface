#ifndef CONTACT_MAPPING_HH
#define CONTACT_MAPPING_HH

#include <vector>
#include <iostream>
#include <tr1/array>

#include <psurface/StaticVector.h>
#include <psurface/PSurface.h>
#include <psurface/ContactToolBox.h>
#include <psurface/IntersectionPrimitive.h>
#include <psurface/IntersectionPrimitiveCollector.h>

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
               ctype epsilon,   ///< The estimate maximum deformation for the contact oracle
               void (*obsDirections)(const double* pos, double* dir)
               );

    void getOverlaps(std::vector<IntersectionPrimitive<1,ctype> >& overlaps)
    {
        IntersectionPrimitiveCollector<ctype>::collect(&psurface_, overlaps);
    }

    // /////////////////////////////////////////////

private:

    PSurface<1,ctype> psurface_;

    std::vector<StaticVector<ctype, 2> > domainNormals;

    std::vector<StaticVector<ctype, 2> > targetNormals;

};


template <class ctype>
class ContactMapping<3,ctype>
{
public:	

    ContactMapping()
        : surface1_(NULL), surface2_(NULL)
    {}
    
    ~ContactMapping();

    void build(const std::vector<std::tr1::array<ctype,3> >& coords1,  ///< The vertices of the first surface as \f$x_0 ,y_0 ,z_0, x_1, y_1, z_1 ...\f$
               const std::vector<std::tr1::array<int,3> >& tri1,       ///< The triangles of the first surface
               const std::vector<std::tr1::array<ctype,3> >& coords2,  ///< The vertices of the second surface
               const std::vector<std::tr1::array<int,3> >& tri2,       ///< The triangles of the second surface
               ctype epsilon,   ///< The estimate maximum deformation for the contact oracle
               void (*obsDirections)(const double* pos, double* dir)
               );

    void getOverlaps(std::vector<IntersectionPrimitive<2,ctype> >& overlaps) {
        IntersectionPrimitiveCollector<ctype>::collect(&psurface_, overlaps);
    }

private:

    PSurface<2,ctype> psurface_;

    Surface* surface1_;
    Surface* surface2_;

};

#endif

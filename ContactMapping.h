#ifndef CONTACT_MAPPING_HH
#define CONTACT_MAPPING_HH

#include <vector>
#include <iostream>
#include <tr1/array>

#include <psurface/StaticVector.h>
#include <psurface/PSurface.h>
#include <psurface/contact.h>

template <int dim, class ctype>
class ContactMapping {};

template <class ctype>
class ContactMapping<2,ctype>
{
public:
    void build(const std::vector<std::tr1::array<double,2> >& coords1,  ///< The vertices of the first surface as \f$x_0 ,y_0 ,z_0, x_1, y_1, z_1 ...\f$
               const std::vector<std::tr1::array<int,2> >& tri1,       ///< The triangles of the first surface
               const std::vector<std::tr1::array<double,2> >& coords2,  ///< The vertices of the second surface
               const std::vector<std::tr1::array<int,2> >& tri2,
               float epsilon,   ///< The estimate maximum deformation for the contact oracle
               void (*obsDirections)(const double* pos, double* dir)
               );

    void getOverlaps(std::vector<IntersectionPrimitive<1,float> >& overlaps);

    // /////////////////////////////////////////////

private:

    PSurface<1,double> psurface_;

    std::vector<StaticVector<double, 2> > domainNormals;

    std::vector<StaticVector<double, 2> > targetNormals;

};


template <class ctype>
class ContactMapping<3,ctype>
{
public:	

    ~ContactMapping() {
        // delete the contact surface after use
        deleteContactSurface();
    }

    void build(const std::vector<std::tr1::array<double,3> >& coords1,  ///< The vertices of the first surface as \f$x_0 ,y_0 ,z_0, x_1, y_1, z_1 ...\f$
               const std::vector<std::tr1::array<int,3> >& tri1,       ///< The triangles of the first surface
               const std::vector<std::tr1::array<double,3> >& coords2,  ///< The vertices of the second surface
               const std::vector<std::tr1::array<int,3> >& tri2,       ///< The triangles of the second surface
               float epsilon,   ///< The estimate maximum deformation for the contact oracle
               void (*obsDirections)(const double* pos, double* dir)
               ) {
        buildContactMapping(coords1, tri1, coords2, tri2,
                            epsilon,
                            obsDirections);
    }

    void getOverlaps(std::vector<IntersectionPrimitive<2,float> >& overlaps) {
        getMergedGrid(overlaps);
    }

};

#endif

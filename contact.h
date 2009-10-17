#ifndef CONTACT_INTERFACE_H
#define CONTACT_INTERFACE_H

/*! \file 
 *
 * This is a wrapper for the standalone contact handling library.  It
 * defines the interfaces to access the ContactSurface class and related
 * facilities.  It is actually quite simple since, even though it's 
 * a complex piece of code, there are few different things you want
 * to do with it.  Thus, few functions suffice.
 */

#include <vector>

#include <psurface/IntersectionPrimitive.h>


/** \brief Compute a mapping in normal direction from surface 1
 */
void buildContactMapping(const std::vector<std::tr1::array<double,3> >& coords1,  ///< The vertices of the first surface as \f$x_0 ,y_0 ,z_0, x_1, y_1, z_1 ...\f$
                         const std::vector<std::tr1::array<int,3> >& tri1,       ///< The triangles of the first surface
                         const std::vector<std::tr1::array<double,3> >& coords2,  ///< The vertices of the second surface
                         const std::vector<std::tr1::array<int,3> >& tri2,       ///< The triangles of the second surface
                         float epsilon,   ///< The estimate maximum deformation for the contact oracle
                         void (*obsDirections)(const double* pos, double* dir)
                         );


/** \brief Computes the overlaps of the images of contact surface triangles on the 
 * nonmortar side.
 *
 * Computes the overlaps of the images of contact surface triangles on the 
 * nonmortar side. They are returned as an array of IntPrimitives.
 */
void getMergedGrid(std::vector<IntersectionPrimitive<2,float> >& overlaps);

/** \brief Properly removes the current intermediate surface.
 *
 * Properly removes the current intermediate surface.  Does not remove the
 * memory allocated in getMergedGrid()!
 */
void deleteContactSurface();

#endif 

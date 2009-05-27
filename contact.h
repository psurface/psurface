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


/** This little struct is used to get the resulting basis function
 * support overlaps out of the library.  It is a direct copy of
 * IntersectionPrimitive.  (In fact, it is supposed to be binary
 * compatible.  I put it here because the actual IntersectionPrimitive
 * class uses stuff from the McLib, and I didn't want to expect
 * calling programs to know about it.
 */
struct IntPrimitive 
{
    
    /** This is the geometric position of the basis function support
     * overlap on the intermediate surface.  It is a little triangle.
     * Therefore, its shape can be described as three vectors
     * from \f$ R^3\f$ each.
     */
    float points[3][3];
    
    
    /** An IntPrimitive always represents the overlap of two basis functions
     * restricted to the image of one mortar and one nonmortar triangle.
     * The indices of those two triangles are given in this array.
     */
    int tris[2];
    
    /** This array marks the exact parts of the mortar and the nonmortar
     * triangles whose overlap is represented by a given IntPrimitive.
     * For each of the two triangles it marks the preimage of the
     * overlap in barycentric coordinates.
     *
     * Interpret it like this:  The first index selects nonmortar resp.
     * mortar side.  The second index tells which of the three points
     * of the IntPrimitive is to be considered.  The third one selects
     * the coefficient in the barycentric coordinate system.
     */
    float localCoords[2][3][2];
    
};

/** \brief Compute a direct contact mapping without intermediate surface
 * an intermediate surface.
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
void getMergedGrid(std::vector<IntPrimitive>& overlaps);

/** \brief Properly removes the current intermediate surface.
 *
 * Properly removes the current intermediate surface.  Does not remove the
 * memory allocated in getMergedGrid()!
 */
void deleteContactSurface();

#endif 

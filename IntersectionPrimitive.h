#ifndef INTERSECTIONPRIMITIVE_H
#define INTERSECTIONPRIMITIVE_H

#include <McVec2f.h>
#include <McVec3f.h>
#include <McSArray.h>


/** This class represents a part of the overlap of a basis function on the
 * mortar and the nonmortar side.  The exact part stored is the following:
 * Consider a FE basis function on the nonmortar side and one on the mortar
 * side.  There supports are sets of triangles on \f$\Gamma_m\f$ resp. \f$\Gamma_n\f$.
 * Pick one triangles from each set.  If their images under the projections
 * \f$\phi_m\f$ resp. \f$\phi_n\f$ overlaps, this overlap will generally be
 * distributed over several triangles of the intermediate surface.  Its
 * restriction to an intermediate surface triangle is a plane convex polygon
 * with no more than six vertices.  Triangulate this polygon and each resulting
 * triangle shall be represented by an IntersectionPrimitive.
 *
 * \see IntPrimitive
 */
class IntersectionPrimitive {

public:


    /** This is the geometric position of the basis function support
     * overlap on the intermediate surface.  It is a little triangle.
     * Therefore, its shape can be described as three vectors
     * from \f$ R^3\f$ each.
     */
    McSArray<McVec3f, 3> points;

    /** An IntersectionPrimitive always represents the overlap of two basis functions
     * restricted to the image of one mortar and one nonmortar triangle.
     * The indices of those two triangles are given in this array.
     */
    McSArray<int, 2> tris;

    /** This array marks the exact parts of the mortar and the nonmortar
     * triangles whose overlap is represented by a given IntersectionPrimitive.
     * For each of the two triangles it marks the preimage of the
     * overlap in barycentric coordinates.
     *
     * Interpret it like this:  The first index selects nonmortar resp.
     * mortar side.  The second index tells which of the three points
     * of the IntersectionPrimitive is to be considered.  
     */
    McSArray<McSArray<McVec2f, 3> , 2> localCoords;

};


#endif

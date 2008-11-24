/////////////////////////////////////////////////////////////////
/*
 * $Id: Parametrization.h,v 1.2 2007/10/18 16:05:58 sander Exp $
 *
 * $Log: Parametrization.h,v $
 * Revision 1.2  2007/10/18 16:05:58  sander
 * stuff from 'contact' merged, renamed to psurface
 *
 * Revision 1.1  2007/10/17 13:16:55  sander
 * moved here from the ZIB server
 *
 * Revision 1.21  2006/03/01 10:26:31  bzfsande
 * include TargetSurface.h when compiled as a standalone library and hxsurface/Surface.h within Amira
 * mailtoauthor: sander@zib.de
 *
 * Revision 1.20  2005/04/26 15:34:42  bzfsande
 * use relative include path
 * mailtoauthor: sander@zib.de
 *
 * Revision 1.19  2004/10/01 13:10:53  bzfsande
 * object duplication implemented, assert replaced by exception-throwing
 * mailtoauthor: sander@zib.de
 *
 * Revision 1.18  2003/10/09 12:07:37  bzfsande
 * different fixes
 *
 * Revision 1.17  2003/08/27 12:22:19  bzfsande
 * fixes
 *
 * Revision 1.16  2003/06/30 12:44:44  bzfsande
 * generalizations needed for the contact library
 *
 * Revision 1.15  2003/06/05 13:01:33  bzfsande
 * introduces ghost nodes as a fifth node type
 * also, the access to the node domain positions is procedural now
 *
 * Revision 1.14  2003/05/09 08:58:20  bzfsande
 * bugfixes
 *
 * Revision 1.13  2003/04/04 14:59:18  bzfsande
 * The base grid is now array-based
 *
 * Revision 1.12  2003/03/24 13:26:53  bzfsande
 * separated the Parametrization object into a base grid object connected to a standard Surface
 *
 * Revision 1.11  2002/10/28 12:18:17  bzfsande
 * fixes
 *
 * Revision 1.10  2002/10/08 13:00:29  bzfsande
 * mapping function with an explicit seed
 *
 * Revision 1.9  2002/10/02 15:22:39  bzfsande
 * Introduced the global #imagePos# array
 * - any data defined on the vertices of the original surface
 *   can now be queried
 * - The base grid vertices don't have to be at the same
 *   position as their homologues on the original surface anymore
 *
 * Revision 1.8  2002/07/25 13:25:15  bzfsande
 * access to the true surface normals
 *
 * Revision 1.7  2001/12/18 15:36:43  bzfsande
 * better smoothing across patch boundaries
 *
 * Revision 1.6  2001/12/02 18:13:24  bzfsande
 * fixes
 *
 * Revision 1.5  2001/11/02 15:28:27  bzflamec
 * compile on windows
 *
 * Revision 1.4  2001/10/05 10:41:19  bzfsande
 * fixes
 *
 * Revision 1.3  2001/09/26 09:33:01  bzfsande
 * edgeIterators added
 *
 * Revision 1.2  2001/09/21 08:06:35  bzfsande
 * .
 *
 * Revision 1.1  2001/09/13 12:39:52  bzfsande
 * initial version
 *
 * Revision 1.10  2001/09/09 15:59:38  bzfsande
 * added boundaryIds and restructured the file format
 *
 * Revision 1.8  2001/08/17 15:37:58  bzfsande
 * more fixes
 *
 * Revision 1.7  2001/08/13 15:34:48  bzfsande
 * more fixes
 *
 * Revision 1.6  2001/08/08 14:41:11  bzfsande
 * old bugs gone - new bugs added
 *
 * Revision 1.5  2001/04/27 13:04:38  bzfsande
 * new stuff:
 * - reader & writer for AmiraMesh
 * - A histogram plotting function for evaluation
 *   the parametrization quality
 * - first attempts on constructing a GLOBAL parametrization
 *
 * Revision 1.4  2001/04/06 12:24:49  bzfsande
 * a test for self-intersections and topology changes, simplification using a 
 * Hausdorff-type distance function and a brand-new internal structure
 *
 * Revision 1.3  2001/02/16 11:53:57  bzfsande
 * added:
 * - the Brown/Faigle point location algorithm
 * - guarantee for foldover-free parametrizations
 * - barycentric parametrization
 * - tons of bugfixes
 *
 * Revision 1.2  2001/01/31 16:20:30  bzfsande
 * a complete, but buggy, version of MAPS
 *
 * Revision 1.1  2000/12/05 13:38:48  bzfsande
 * MAPS
 *
 */
/////////////////////////////////////////////////////////////////
#ifndef PARAMETRIZATION_H
#define PARAMETRIZATION_H

#include "hxsurface/Surface.h"

#include <amiramesh/HxParamBundle.h>

#include "psurfaceAPI.h"

#include <psurface/McSurfaceBase.h>
#include <psurface/Domains.h>
#include <psurface/SurfacePathSet.h>
#include <psurface/GlobalNodeIdx.h>
#include <psurface/NodeBundle.h>

class AmiraMesh;

/** The parametrization of an arbitrary surface over a simple base grid */
class PSURFACE_API Parametrization : public McSurfaceBase<DomainVertex, DomainEdge, DomainTriangle>{

public:

    /** \brief An exception being thrown by code in this class */
    class ParamError {};

    /** This little class encapsulates all data about surface patches 
     */
    class Patch {
    public:

        Patch() : innerRegion(0), outerRegion(0), boundaryId(0) {}

        /// The number of the material on one side of the triangle patch
        int innerRegion;

        /// The number of the material on the other side of the patch
        int outerRegion;

        /// An additional ID storing some form of material property
        int boundaryId;
    };

    /// Default constructor
    Parametrization(HxParamBundle* bundle=NULL);

    /// Destructor
    virtual ~Parametrization();

    /** \brief Initializes the parametrization from a given second one.
     */
    void init(const Parametrization* other);

    /** \brief Empties the Parametrization object
     */
    void clear();

    /// Get box containing all vertices.
    virtual void getBoundingBox(McBox3f& bbox) const;

    /** \brief Sets up the internal data structures needed by the map() method.
     * 
     * The map() method uses the Brown/Faigle algorithm for quick point location
     * in plane triangulations.  This algorithm requires the data structure for
     * the plane triangulation to comply with certain requirements.  In particular,
     * the edges that leave a vertex have to be in a cyclic order.  Thus
     * before usage of the parametrization, this routine has to be called once.
     */
    void createPointLocationStructure();

    /// Returns the number of patches.
    int numPatches() const {
        return patches.size();
    }

    int getNumNodes() const;

    /** \brief Returns the number of node targets
     *
     * \deprecated
     */
    int numNodes() const {
        return iPos.size();
     }

    /// computes the number of true nodes
    int getNumTrueNodes();

    /** \brief Cleans up the internal data structure.  
     *
     * It does the same as McSurfaceBase::garbageCollection() but it also
     * updates the paths.
     */
    void garbageCollection();


    /// Writes the parametrization in AmiraMesh format
    static int writeAmiraMesh(Parametrization* par, const char* filename);

    /// Reads the parametrization from an AmiraMesh object
    static Parametrization* readAmiraMesh(AmiraMesh* am, const char* filename);

    /// AmiraMesh Reader using an existing AmiraMesh object (is used by derived classes)
    bool initFromAmiraMesh(AmiraMesh* am, const char* filename, Surface* surf);

    /** \brief Creates an explicit Surface object from the information implicitly
     * given by the parametrization.
     */
    void setupOriginalSurface();

    /** \brief Copies the path information from the SurfacePathSet to the
     * parameters.
     *
     * Copies the path information from the SurfacePathSet to the
     * parameters.  This routine is   called before writing a Parametrization
     * to disk, because on disk paths are stored together with the parameters.
     */
    void savePaths(HxParamBundle& parameters);

    /** \brief Reads paths from the parameters HxParamBundle into the
     * Parametrization's own SurfacePathSet.
     */
    void getPaths(const HxParamBundle& parameters);

    /** \brief Adds the triangular closure on each triangle.
     *
     * The image surface of the parametrization consists solely of triangles.
     * However, the preimage of such a triangle, restricted to a domain grid
     * triangle may be a quadrangle.  Since many operations on a base grid
     * triangle expect the graph there to be only triangles we triangulate
     * the graph using this routine.
     *
     * \bug When Ghost Nodes are present, more complex polygons than just
     * quadrangles can occur.  This routine is presently not able to handle them.
     */
    void insertExtraEdges();


    /// Removes the triangular closure on each triangle
    void removeExtraEdges();

    /// Returns a reference to the node specified by a GlobalNodeIdx
    Node& nodes(const GlobalNodeIdx& n) {return triangles(n.tri).nodes[n.idx];}

    /// Returns a const reference to the node specified by a GlobalNodeIdx
    const Node& nodes(const GlobalNodeIdx& n) const {return triangles(n.tri).nodes[n.idx];}
    /** \brief Procedural interface to the target Position in \f$R^3\f$ of a node.
     */
    McVec3f imagePos(const GlobalNodeIdx& n) const {
        return imagePos(n.tri, n.idx);
    }

    /** \brief Procedural interface to the target Position in \f$R^3\f$ of a node.
     */
    McVec3f imagePos(TriangleIdx tri, NodeIdx node) const {
        const Node& cN = triangles(tri).nodes[node];

        switch (cN.type) {
        case Node::GHOST_NODE: {
            const Surface::Triangle& cT = surface->triangles[cN.getNodeNumber()];
            return PlaneParam::linearInterpol(cN.dP, surface->points[cT.points[0]],
                                              surface->points[cT.points[1]], surface->points[cT.points[2]]);
        }
        case Node::INTERSECTION_NODE:
            return iPos[cN.getNodeNumber()];

        default:
            return surface->points[cN.getNodeNumber()];
        }

    }

    /** \brief Returns the local coordinates of the image of a node with
     * respect to a given image triangle.
     */
    McVec2f getLocalTargetCoords(const GlobalNodeIdx& n, int targetTri) const;

    /**@name Mapping functions */
    //@{

    /** \brief The main routine to evaluate the parametrization function.
     *
     * The routine evaluates the parametrization function embodied by the
     * Parametrization class.  It takes a point \f$x\f$ on the domain grid (specified
     * by giving a triangle number and barycentric coordinates on that triangle)
     * and returns the corresponding image point \f$\phi(x)\f$.  That image point is 
     * given via the three vertices of the triangle \f$t\f$ it's on and its position
     * on \f$t\f$ in barycentric coordinates.  This general way of returning \f$\phi(x)\f$
     * allows to correctly interpolate data given on the target surface.
     * If you just want the position of \f$\phi(x)\f$ in \f$R^3\f$ or the surface
     * normal at that point you can call the appropriate convenience functions.
     *
     * @return <tt>true</tt> if everything went correctly, <tt> false</tt> if not.
     *
     * \todo Make this routine GHOST_NODE-proof 
     * \todo Document seed!
     * 
     */
    int map(TriangleIdx tri,                ///< The triangle of the input point \f$x\f$
            McVec2f& p,                      ///< The barycentric coordinates of \f$x\f$ with respect to tri
            McVec3i& vertices,               ///< Return value: The three vertices of the triangle that \f$\phi(x)\f$ is on
            McVec2f& coords,                 ///< The barycentric coordinates of \f$\phi(x)\f$ wrt <tt>vertices</tt>
            int seed=-1                      
            ) const;


    /** \brief Convenience function for accessing the position of points on the target surface.
     *
     * Given a point \f$x\f$ on the base grid, this routine returns the position
     * of \f$\phi(x)\f$. 
     * 
     * Internally, this routine calls map().  Thus, it is not faster than calling
     * map() and computing the position oneself.
     *
     * @return <tt>true</tt> if everything went correctly, <tt> false</tt> if not.
     */
    int positionMap(TriangleIdx tri, McVec2f& p, McVec3f& result) const;

    /** \brief Convenience function for accessing the normals of the target surface.
     *
     * Given a point \f$x\f$ on the base grid, this routine returns the normal
     * of the target surface at \f$\phi(x)\f$.  It returns the direct triangle
     * normal, i.e., the normals are constant on each target triangle.
     * 
     * Internally, this routine calls map().  Thus, it is not faster than calling
     * map() and computing the normal oneself.
     *
     * @return <tt>true</tt> if everything went correctly, <tt> false</tt> if not.
     */
    int directNormalMap(TriangleIdx tri, McVec2f& p, McVec3f& result) const;

    //@}



    /** \brief Given three nodes that form a triangle \f$t\f$ on a base grid triangle,
     * this routine returns the preimage nodes of the target triangle \f$T\f$ that
     * \f$t\f$ is part of.
     */
    void getActualVertices(TriangleIdx tri,                     
                           const McSArray<NodeIdx, 3>& nds,
                           McSArray<GlobalNodeIdx, 3>& vertices   ///< The result nodes are returned here
                           ) const;

    /** \brief Given three nodes that form a triangle \f$t\f$ on a base grid triangle,
     * this routine returns the index of the target triangle \f$T\f$ that
     * \f$t\f$ is part of.
     *
     * \return A triangle index, -1 if the input data was invalid.
     */
    int getImageSurfaceTriangle(TriangleIdx tri,
                                const McSArray<NodeIdx, 3>& nds
                                ) const;

    /** \brief Returns the set of all target triangles that contain the image of a node.
     *
     * <b> Warning: </b> Shouldn't be called for an INTERSECTION_NODE
     * \todo Make this callable for INTERSECTION_NODEs.
     * <b> Warning: </b> The method assumes that Surface::computeTrianglesPerPoint has
     * been called before.
     */
    McSmallArray<int, 6> getTargetTrianglesPerNode(const GlobalNodeIdx& n) const;



    /** \brief Inverses the orientation of all triangles (of a given patch)
     *
     * Inverses the orientation of the base grid triangles while maintaining
     * the parametrization consistent.  If a patch number is given, only triangles
     * of that patch are inverted, if not all triangles are affected.
     *
     * @return Number of flipped triangles.
     */
    int invertTriangles(int patch=-1);


    /** \brief Finds the endnode of the edge an INTERSECTION_NODE is on.
     *
     * This routine returns different things depending on the type of the Node you
     * call it for.
     * <ul>
     * <li> If it is not an INTERSECTION_NODE it just returns the input arguments.
     * <li> If it is an INTERSECTION_NODE, the algorithm starts looking along the
     * edge the input node is on, until it finds a nonINTERSECTION_NODE.  That node
     * is the end of the edge in that direction.  The search direction is the one
     * away from the input triangle, without first traversing it.
     * </ul>
     */
    GlobalNodeIdx getOtherEndNode(TriangleIdx tri, NodeIdx cN) const;


    /** \brief Tests the object for consistency.
     *
     * This routine is intended for debugging purposes.  It performs an
     * internal consistency check.  If a part of the check fails it
     * prints a (somewhat) meaningfull error message, the string in
     * <tt>where</tt> to tell you where the check was called from
     * and then evokes a SIGABRT.
     *
     * If NDEBUG is set checkConsistency() is completely empty.
     */
    void checkConsistency(const char* where) const;

    NodeIdx addNode(TriangleIdx tri, const McVec3f& p);

    NodeIdx addInteriorNode(TriangleIdx tri, const McVec2f& dom, int nodeNumber);

    NodeIdx addGhostNode(TriangleIdx tri, int corner, int targetTri, const McVec2f& localTargetCoords);

    NodeIdx addCornerNode(TriangleIdx tri, int corner, int nodeNumber);

    /** \todo Sollte vielleicht ein Bundle zurückgeben */
    NodeIdx addIntersectionNodePair(TriangleIdx tri1, TriangleIdx tri2,
                                    const McVec2f& dP1, const McVec2f& dP2, 
                                    int edge1, int edge2, const McVec3f& range);

    NodeIdx addTouchingNode(TriangleIdx tri, const McVec2f& dP, int edge, int nodeNumber);

    /** \todo Sollte vielleicht ein Bundle zurückgeben */
    NodeIdx addTouchingNodePair(TriangleIdx tri1, TriangleIdx tri2,
                                const McVec2f& dP1, const McVec2f& dP2, 
                                int edge1, int edge2, int nodeNumber);

    void addParTriangle(TriangleIdx tri, const McVec3i& p);
    
    /** Not a fast implementation, but it works even if the edgePoint arrays
     * haven't been set up yet.
     */
    NodeBundle getNodeBundleAtVertex(VertexIdx v) const;

protected:

    /// This is a service routine only for getTargetTrianglesPerNode
    void getTrianglesPerEdge(int from, int to, McSmallArray<int, 6>& tris, int exception) const;

    /** \brief Internal routine used by map() 
     */
    void handleMapOnEdge(TriangleIdx tri, const McVec2f& p, const McVec2f& a, const McVec2f& b,
                         int edge, int edgePos, McSArray<GlobalNodeIdx, 3>& vertices, McVec2f& coords) const;

    /** \brief Internal routine used by setupOriginalSurface() */
    void appendTriangleToOriginalSurface(const McVec3i& v, int patch);


    /////////////////////////////////////////////////////
    // Data members

public:
    /** This flag needs to set if the the map()-method on a DomainTriangle is used
        It is set by calling the method createPointLocationStructure().
        You should not set it deliberately unless you know what you're doing.*/
    bool hasUpToDatePointLocationStructure;


    /// All base grid patches
    McDArray<Patch> patches;

    /// A set of arbitrary parameters
    HxParamBundle* params;

    bool hasOwnParamBundle;

    /// The image positions of all nodes \deprecated To be replaced by a procedural interface
    McDArray<McVec3f> iPos;

    /// The corresponding image surface
    Surface* surface;

    /// Surface paths on the base grid
    SurfacePathSet paths;

};


#endif


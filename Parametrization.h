/////////////////////////////////////////////////////////////////
/*
 * $Id: Parametrization.h,v 1.2 2007/10/18 16:05:58 sander Exp $
 */
/////////////////////////////////////////////////////////////////
#ifndef PARAMETRIZATION_H
#define PARAMETRIZATION_H

#ifdef PSURFACE_STANDALONE
#include "TargetSurface.h"
#else
#include "hxsurface/Surface.h"
#endif

#include <amiramesh/HxParamBundle.h>

#include "psurfaceAPI.h"

#include "McSurfaceBase.h"
#include "Domains.h"
#include "SurfacePathSet.h"
#include "GlobalNodeIdx.h"
#include "NodeBundle.h"

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
    virtual void getBoundingBox(Box<float,3>& bbox) const;

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
    McVec3f imagePos(int tri, NodeIdx node) const {
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
    int map(int tri,                ///< The triangle of the input point \f$x\f$
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
    int positionMap(int tri, McVec2f& p, McVec3f& result) const;

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
    int directNormalMap(int tri, McVec2f& p, McVec3f& result) const;

    //@}



    /** \brief Given three nodes that form a triangle \f$t\f$ on a base grid triangle,
     * this routine returns the preimage nodes of the target triangle \f$T\f$ that
     * \f$t\f$ is part of.
     */
    void getActualVertices(int tri,                     
                           const McSArray<NodeIdx, 3>& nds,
                           McSArray<GlobalNodeIdx, 3>& vertices   ///< The result nodes are returned here
                           ) const;

    /** \brief Given three nodes that form a triangle \f$t\f$ on a base grid triangle,
     * this routine returns the index of the target triangle \f$T\f$ that
     * \f$t\f$ is part of.
     *
     * \return A triangle index, -1 if the input data was invalid.
     */
    int getImageSurfaceTriangle(int tri,
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
    GlobalNodeIdx getOtherEndNode(int tri, NodeIdx cN) const;


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

    NodeIdx addNode(int tri, const McVec3f& p);

    NodeIdx addInteriorNode(int tri, const McVec2f& dom, int nodeNumber);

    NodeIdx addGhostNode(int tri, int corner, int targetTri, const McVec2f& localTargetCoords);

    NodeIdx addCornerNode(int tri, int corner, int nodeNumber);

    /** \todo Sollte vielleicht ein Bundle zurückgeben */
    NodeIdx addIntersectionNodePair(int tri1, int tri2,
                                    const McVec2f& dP1, const McVec2f& dP2, 
                                    int edge1, int edge2, const McVec3f& range);

    NodeIdx addTouchingNode(int tri, const McVec2f& dP, int edge, int nodeNumber);

    /** \todo Sollte vielleicht ein Bundle zurückgeben */
    NodeIdx addTouchingNodePair(int tri1, int tri2,
                                const McVec2f& dP1, const McVec2f& dP2, 
                                int edge1, int edge2, int nodeNumber);

    void addParTriangle(int tri, const McVec3i& p);
    
    /** Not a fast implementation, but it works even if the edgePoint arrays
     * haven't been set up yet.
     */
    NodeBundle getNodeBundleAtVertex(int v) const;

protected:

    /// This is a service routine only for getTargetTrianglesPerNode
    void getTrianglesPerEdge(int from, int to, McSmallArray<int, 6>& tris, int exception) const;

    /** \brief Internal routine used by map() 
     */
    void handleMapOnEdge(int tri, const McVec2f& p, const McVec2f& a, const McVec2f& b,
                         int edge, int edgePos, McSArray<GlobalNodeIdx, 3>& vertices, McVec2f& coords) const;

    /** \brief Internal routine used by setupOriginalSurface() */
    void appendTriangleToOriginalSurface(const McVec3i& v, int patch);


    /////////////////////////////////////////////////////
    // Data members
    /////////////////////////////////////////////////////

public:
    /** This flag needs to set if the the map()-method on a DomainTriangle is used
        It is set by calling the method createPointLocationStructure().
        You should not set it deliberately unless you know what you're doing.*/
    bool hasUpToDatePointLocationStructure;


    /// All base grid patches
    std::vector<Patch> patches;

    /// A set of arbitrary parameters
    HxParamBundle* params;

    bool hasOwnParamBundle;

    /// The image positions of all nodes \deprecated To be replaced by a procedural interface
    std::vector<McVec3f> iPos;

    /// The corresponding image surface
    Surface* surface;

    /// Surface paths on the base grid
    SurfacePathSet paths;

};


#endif


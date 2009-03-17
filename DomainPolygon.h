#ifndef DOMAIN_POLYGON_H
#define DOMAIN_POLYGON_H

#include <psurface/Domains.h>
#include <psurface/PlaneParam.h>
#include <psurface/Parametrization.h>
#include "psurfaceAPI.h"

class DomainTriangle;
class DomainVertex;
class DomainEdge;
class CircularPatch;


/** A polygon carrying a plane triangulation */
class PSURFACE_API DomainPolygon : public PlaneParam{
public:
    /// standard constructor
    DomainPolygon(Parametrization* _par) : par(_par) {};

    ///
    ~DomainPolygon() {};

    /// use of assignment operator and copy constructor is blocked
private:
    DomainPolygon& operator=(const DomainPolygon& rhs) {return *this;}

    DomainPolygon(const DomainPolygon& other) {}

public:
        
    ///
    void init(const DomainTriangle& tri, const McVec2f coords[3]);


    enum NodeLocation {IN_TRIANGLE, ON_SEGMENT, IN_POLYGON};
    ///
    NodeIdx cornerNode(int i) const {
        
        return edgePoints[i][0];
    }

    ///
    NodeIdx getNextEdgeNode(NodeIdx n) const {
        assert(!nodes[n].isINTERIOR_NODE());

        int edge    = nodes[n].getDomainEdge();
        int edgePos = nodes[n].getDomainEdgePosition();

        if (edgePos!=edgePoints[edge].size()-1) {
            return edgePoints[edge][edgePos+1];
        } else {
            int edge = (nodes[n].getDomainEdge()+1)%boundaryPoints.size();
            return edgePoints[edge][1];
        }
    }

    ///
    NodeIdx getPreviousEdgeNode(NodeIdx n) const {
       assert(!nodes[n].isINTERIOR_NODE());

        int edgePos = nodes[n].getDomainEdgePosition();

        if (edgePos!=0) {
            return edgePoints[nodes[n].getDomainEdge()][edgePos-1];
        } else {
            int edge = (nodes[n].getDomainEdge()+boundaryPoints.size()-1)%boundaryPoints.size();
            return edgePoints[edge][edgePoints[edge].size()-2];
        }
    }

    ///
    void unflipTriangles() {
        PlaneParam::unflipTriangles(par->iPos);
    }

    ///
    void applyParametrization() {
        PlaneParam::applyParametrization(par->iPos);
    }

    ///
    void garbageCollection() {
        McDArray<int> offArr;
        garbageCollection(offArr);
    }

    ///
    void garbageCollection(McDArray<int>& offArr);

    ///
    void createPointLocationStructure();

    ///
    void insertExtraEdges();
        
    /// merges the polygon with a triangle
    void mergeTriangle(int tri, McVec2f coords[3], int& newCenterNode,
                       McDArray<unsigned int>& nodeStack);

    /** \brief This uses a given triangulation to cut a DomainPolygon into a
     * set of DomainTriangles.
     *
     */
    bool triangulate(CircularPatch &fillIn, 
                     McDArray<unsigned int>& nodeStack);

    ///
    void cutParameterEdges(int boundaryIdx, NodeIdx startNode, NodeIdx lastNode,
                           McDArray<int>& nodeLocs,
                           DomainTriangle& cT,
                           const McSArray<McVec2f, 3>& newTriangleCoords,
                           McSmallArray<int, 2>& triNewEdgePoints,
                           McSmallArray<int, 2>& polyNewEdgePoints,
                           McDArray<unsigned int>& nodeStack);

    ///
    NodeIdx splitNode(NodeIdx cN, McDArray<int>& nodeLocs);

    ///
    unsigned int createNodePosition(McDArray<McVec3f>& nodePositions, McDArray<unsigned int>& nodeStack,
                                    const McVec3f& newImagePos);

    /// removes a vertex from the polygon
    void removeVertex(int point);

    /// does a cut from a given node to the first boundaryVertex
    void slice(int centerNode, int centerVertex, int bVertex);

    /** \brief Determines the intersection point of \f$(p1, p2) \f$ and (p3, p4) as the affine
     *  combination of p1 and p2 
     *
     * \todo Rewrite this */
    float computeIntersection(float &mu, const McVec2f &p1, const McVec2f &p2, 
                              const McVec2f &p3, const McVec2f &p4);

    /** \brief Transforms the domain positions into a world coordinate system
     *
     * Assuming that the domain positions of the nodes contained in the
     * PlaneParam underlying this DomainPolygon are given in barycentric 
     * coordinates, this function 
     * transforms them into a world coordinates.  The new coordinate system
     * is specified by supplying new coordinates for the three points
     * (1,0), (0,1), and (0, 0).  The transformation can be restricted
     * to nodes appearing in the local node array with a fixed index
     * or higher.
     *
     * \param newNodeIdx 
     * \param a, b, c : The new coordinates of the points \f$(1,0)\f$, 
     * \f$(0,1)\f$, and \f$(0,0)\f$, respectively.
     */
    void installWorldCoordinates(int newNodeIdx, const McVec2f &a, const McVec2f &b, const McVec2f &c){
        for (int i=newNodeIdx; i<nodes.size(); i++)
            nodes[i].setDomainPos(a*nodes[i].domainPos().x + b*nodes[i].domainPos().y + 
                c*(1-nodes[i].domainPos().x-nodes[i].domainPos().y));
    }

    void augmentNeighborIdx(int newNodeIdx, McSmallArray<int, 2> edgePoints[3]) {
        int i,d = newNodeIdx;

        for (i=newNodeIdx; i<nodes.size(); i++)
            for (int j=0; j<nodes[i].degree(); j++)
                nodes[i].neighbors(j) += d;

        for (i=0; i<3; i++)
            for (int j=0; j<edgePoints[i].size(); j++)
                edgePoints[i][j] += d;
    }

    void updateEdgePoints(McSmallArray<int, 2> edgePoints[3], int oldNode, int newNode) {
        int i;
        for (i=0; i<3; i++){
            if (edgePoints[i][0]==oldNode)
                edgePoints[i][0] = newNode;
            if (edgePoints[i].last()==oldNode)
                edgePoints[i].last() = newNode;
        }
    }


    /**@name debug code */
    //@{
    /// prints info about the triangle (to stdout)
    void print(bool showEdgePoints=false, bool showParamEdges=false, bool showNodes=false) const;

    /// 
    void checkConsistency(const char* where);
    //@}

    ////////////////////////////////////////////

    /// a list of the boundary points
    McDArray<int> boundaryPoints;

    /// for each edge of the polygon, a list of the mapping nodes that are right on this edge
    McDArray<McSmallArray<int, 2> > edgePoints;

    Parametrization* par;

};


#endif

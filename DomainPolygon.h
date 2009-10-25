#ifndef DOMAIN_POLYGON_H
#define DOMAIN_POLYGON_H

#include <psurface/Domains.h>
#include <psurface/PlaneParam.h>
#include <psurface/PSurface.h>

class DomainTriangle;
class DomainVertex;
class DomainEdge;
class CircularPatch;


/** A polygon carrying a plane triangulation */
class DomainPolygon : public PlaneParam<float> {
public:
    /// standard constructor
    DomainPolygon(PSurface<2,float>* _par) : par(_par) {};

    ///
    ~DomainPolygon() {};

    /// use of assignment operator and copy constructor is blocked
private:
    DomainPolygon& operator=(const DomainPolygon& rhs) {return *this;}

    DomainPolygon(const DomainPolygon& other) {}

public:
        
    ///
    void init(const DomainTriangle& tri, const StaticVector<float,2> coords[3]);


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
        PlaneParam<float>::unflipTriangles(par->iPos);
    }

    ///
    void applyParametrization() {
        PlaneParam<float>::applyParametrization(par->iPos);
    }

    ///
    void garbageCollection() {
        std::vector<int> offArr;
        garbageCollection(offArr);
    }

    ///
    void garbageCollection(std::vector<int>& offArr);

    ///
    void createPointLocationStructure();

    ///
    void insertExtraEdges();
        
    /// merges the polygon with a triangle
    void mergeTriangle(int tri, StaticVector<float,2> coords[3], int& newCenterNode,
                       std::vector<unsigned int>& nodeStack);

    /** \brief This uses a given triangulation to cut a DomainPolygon into a
     * set of DomainTriangles.
     *
     */
    bool triangulate(CircularPatch &fillIn, 
                     std::vector<unsigned int>& nodeStack);

    ///
    void cutParameterEdges(int boundaryIdx, NodeIdx startNode, NodeIdx lastNode,
                           std::vector<int>& nodeLocs,
                           DomainTriangle& cT,
                           const std::tr1::array<StaticVector<float,2>, 3>& newTriangleCoords,
                           std::vector<int>& triNewEdgePoints,
                           std::vector<int>& polyNewEdgePoints,
                           std::vector<unsigned int>& nodeStack);

    ///
    NodeIdx splitNode(NodeIdx cN, std::vector<int>& nodeLocs);

    ///
    unsigned int createNodePosition(std::vector<StaticVector<float,3> >& nodePositions, std::vector<unsigned int>& nodeStack,
                                    const StaticVector<float,3>& newImagePos);

    /// removes a vertex from the polygon
    void removeVertex(int point);

    /// does a cut from a given node to the first boundaryVertex
    void slice(int centerNode, int centerVertex, int bVertex);

    /** \brief Determines the intersection point of \f$(p1, p2) \f$ and (p3, p4) as the affine
     *  combination of p1 and p2 
     *
     * \todo Rewrite this */
    float computeIntersection(float &mu, const StaticVector<float,2> &p1, const StaticVector<float,2> &p2, 
                              const StaticVector<float,2> &p3, const StaticVector<float,2> &p4);

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
    void installWorldCoordinates(int newNodeIdx, const StaticVector<float,2> &a, const StaticVector<float,2> &b, const StaticVector<float,2> &c){
        for (int i=newNodeIdx; i<nodes.size(); i++)
            nodes[i].setDomainPos(a*nodes[i].domainPos()[0] + b*nodes[i].domainPos()[1] + 
                c*(1-nodes[i].domainPos()[0]-nodes[i].domainPos()[1]));
    }

    void augmentNeighborIdx(int newNodeIdx, std::tr1::array<std::vector<int>, 3>& edgePoints) {
        int i,d = newNodeIdx;

        for (i=newNodeIdx; i<nodes.size(); i++)
            for (int j=0; j<nodes[i].degree(); j++)
                nodes[i].neighbors(j) += d;

        for (i=0; i<3; i++)
            for (int j=0; j<edgePoints[i].size(); j++)
                edgePoints[i][j] += d;
    }

    void updateEdgePoints(std::tr1::array<std::vector<int>,3>& edgePoints, int oldNode, int newNode) {
        int i;
        for (i=0; i<3; i++){
            if (edgePoints[i][0]==oldNode)
                edgePoints[i][0] = newNode;
            if (edgePoints[i].back()==oldNode)
                edgePoints[i].back() = newNode;
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
    std::vector<int> boundaryPoints;

    /// for each edge of the polygon, a list of the mapping nodes that are right on this edge
    std::vector<std::vector<int> > edgePoints;

    PSurface<2,float>* par;

};


#endif

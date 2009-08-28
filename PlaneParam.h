/////////////////////////////////////////////////////////////////
/*
 * $Id: PlaneParam.h,v 1.2 2007/10/18 16:05:58 sander Exp $
 */
/////////////////////////////////////////////////////////////////
#ifndef PLANE_PARAM_H
#define PLANE_PARAM_H

#include <psurface/StaticVector.h>
#include <vector>
#include <algorithm>

#include "Node.h"

template<class T> class SparseMatrix;

class Parametrization;

/** The parametrization of a surface in space.  It is actually a planar graph
   having only triangular and unbounded faces.  The nodes of the graph store
   a 3D position, the rest is computed via linear interpolation. */
class PlaneParam{
public:

    class DirectedEdgeIterator {
    public:
        DirectedEdgeIterator() : fromNode(-1), neighborIdx(0) , nodes(0){}

        DirectedEdgeIterator(const std::vector<Node>& _nodes) :  fromNode(-1), neighborIdx(0) {
            nodes = &_nodes;
        }

        int from() const {return fromNode;}

        int to() const {return (*nodes)[fromNode].neighbors(neighborIdx);}

        bool isValid() const {
            return fromNode>=0 && fromNode<nodes->size();
        }

        ///
        DirectedEdgeIterator& operator++();

        ///
        void invert();

        ///
        DirectedEdgeIterator getONext() const {
            DirectedEdgeIterator dest = *this;
            dest.fromNode = fromNode;
            dest.neighborIdx = (neighborIdx+1)%(*nodes)[fromNode].degree();

            return dest;
        }

        ///
        DirectedEdgeIterator getDPrev() const {
            DirectedEdgeIterator dest = *this;
            
            dest.invert();
            dest.neighborIdx = (dest.neighborIdx + (*nodes)[dest.fromNode].degree()-1)%(*nodes)[dest.fromNode].degree();
            dest.invert();

            return dest;
        }

        int fromNode;
        int neighborIdx;

        const std::vector<Node>* nodes;

    };

    class UndirectedEdgeIterator {
    public:
        UndirectedEdgeIterator() : fromNode(-1), neighborIdx(0) , nodes(0){}

        UndirectedEdgeIterator(const std::vector<Node>& _nodes) :  fromNode(-1), neighborIdx(0) {
            nodes = &_nodes;
        }

        NodeIdx from() const {return fromNode;}

        int to() const {return (*nodes)[fromNode].neighbors(neighborIdx);}

        bool isValid() const {
            return fromNode>=0 && fromNode<nodes->size();
        }

        bool isRegularEdge() const {
            return (*nodes)[fromNode].neighbors(neighborIdx).isRegular();
        }

        UndirectedEdgeIterator& operator++();

        bool isCorrectlyOriented() const {return from() < to();}

        NodeIdx fromNode;
        int neighborIdx;

        const std::vector<Node>* nodes;

    };

    /** \brief An iterator over all triangles in a plane graph */
    class TriangleIterator {
    public:

        /** \brief Default constructor */
        TriangleIterator() {}

        TriangleIterator(const DirectedEdgeIterator& firstEdge);

        /** \brief Prefix increment */
        TriangleIterator& operator++();

        bool isValid() const {return cE.isValid();}

        /** \brief Access to the three vertices of the triangle.  

        It's read-only
         * because this class is an *iterator*, and not an actual triangle
         * object. */
        NodeIdx vertices(int i) const {
            assert(0<=i && i<3);
            if (i==0)
                return cE.from();
            else if (i==1)
                return cE.to();
            else
                return cE.getONext().to();
        }

        /** \brief Make a copy of the three vertex indices  */
        std::tr1::array<NodeIdx, 3> vertices() const {
            std::tr1::array<NodeIdx, 3> result;
            result[0] = vertices(0);
            result[1] = vertices(1);
            result[2] = vertices(2);
            return result;
        }

    protected:
        bool isCorrectlyOriented() const;

        DirectedEdgeIterator cE;
    };

    /**@name Access to different iterators */
    //@{

    TriangleIterator firstTriangle() const;

    UndirectedEdgeIterator firstUndirectedEdge() const {

        UndirectedEdgeIterator edge(nodes);

        //firstPseudoEdge(edge);
        if (nodes.size()==0){
            edge.fromNode = -1;
            return edge;
        } else
            edge.fromNode = 0;

        while (!nodes[edge.fromNode].degree()){
            edge.fromNode++;
            if (!edge.isValid())
                return edge;
        }

        edge.neighborIdx = 0;

        ////////////////////////////////
        // if (!edge.isValid(nodes))
//          return;

        while (!edge.isCorrectlyOriented()){
            //nextPseudoEdge(edge);
            if (edge.neighborIdx < nodes[edge.from()].degree()-1)
                edge.neighborIdx++;
            else {
                do{
                    edge.fromNode++;
                    if (!edge.isValid())
                        return edge;
                }while (!nodes[edge.fromNode].degree());
                
                edge.neighborIdx = 0;
            }
            
            /////////////////////////
            if (!edge.isValid())
                return edge;
        }
        
        return edge;
        
    }


    /** \brief Creates a directed edge iterator
     *
     * \param n Creates an edge iterator connected to this node
     */
    DirectedEdgeIterator firstDirectedEdge(NodeIdx n=0) const {

        DirectedEdgeIterator edge(nodes);

        if (n<0 || n>=nodes.size()){
            edge.fromNode = -1;
            return edge;
        } else
            edge.fromNode = n;

        while (!nodes[edge.fromNode].degree()){
            edge.fromNode++;
            if (!edge.isValid())
                return edge;
        }

        edge.neighborIdx = 0;

        return edge;
    }

    DirectedEdgeIterator getDirectedEdgeIterator(int from, int to) const {

        DirectedEdgeIterator edge(nodes);

        assert(from>=0 && from<nodes.size() && to>=0 && to<nodes.size());

        edge.fromNode = from;

        //edge.neighborIdx = mcSmallArray::index(nodes[from].nbs, to);
//         if (edge.neighborIdx == -1)
//             edge.fromNode = -1;

        std::vector<Node::NeighborReference>::const_iterator findIt = std::find(nodes[from].nbs.begin(), nodes[from].nbs.end(), to);

        if (findIt ==nodes[from].nbs.end()) {
            edge.neighborIdx = -1;
            edge.fromNode = -1;
        } else
            edge.neighborIdx = findIt-nodes[from].nbs.begin();
        

        return edge;
    }

    //@}
    /////////////////////////////////////////////////////
    // The actual parametrization

    ///
    PlaneParam() {}

    ///
    ~PlaneParam();

    ///
    void makeOneTriangle(int a, int b, int c)
    {
        nodes.resize(3);

        nodes[0].setValue(StaticVector<float,2>(1, 0), a, Node::CORNER_NODE);
        nodes[1].setValue(StaticVector<float,2>(0, 1), b, Node::CORNER_NODE);
        nodes[2].setValue(StaticVector<float,2>(0, 0), c, Node::CORNER_NODE);

        addEdge(0, 1);
        addEdge(1, 2);
        addEdge(2, 0);

    }

    /** \brief Adds an edge to the graph
     *
     * This routine adds an edge to the planar graph.  It does not check
     * whether the edge exists already.  Thus, multiedges are possible. */
    void addEdge(int from,    //!< First endnode of the new edge
                 int to,       //!< Second endnode of the new edge
                 bool triangularClosure=false   /**< This flag marks whether an edge actually 
                                                 * appears in the original surface or whether 
                                                 * it has just been inserted to turn the 
                                                 * graph into a triangulation. */
                 ){
        nodes[from].appendNeighbor(Node::NeighborReference(to, triangularClosure));
        nodes[to].appendNeighbor(Node::NeighborReference(from, triangularClosure));
    }

    ///
    void removeEdge(int from, int to){
        nodes[from].removeReferenceTo(to);
        nodes[to].removeReferenceTo(from);
    }

    ///
    unsigned int getNumEdges() const {
        int n=0;
        for (size_t i=0; i<nodes.size(); i++)
            n += nodes[i].degree();

        return n/2;
    }

    ///
    unsigned int getNumRegularEdges() const {
        int n=0;
        for (int i=0; i<nodes.size(); i++)
            for (int j=0; j<nodes[i].degree(); j++)
                if (nodes[i].neighbors(j).isRegular())
                    n++;

        return n/2;
    }


    
    ///
    void countNodes(int& intersectionNodes, int& touchingNodes, int& interiorNodes) const {
        intersectionNodes = touchingNodes = interiorNodes = 0;

        for (int i=0; i<nodes.size(); i++){
            switch (nodes[i].getType()) {
            case Node::INTERSECTION_NODE:
                intersectionNodes++;
                break;
            case Node::TOUCHING_NODE:
                touchingNodes++;
                break;
            case Node::INTERIOR_NODE:
                interiorNodes++;
                break;
            }
        }
    }

    ///
    void mergeNodes(int one, int other) {

        int i, j, k;

        // remove mutual references
        for (i=nodes[other].degree()-1; i>=0; i--)
            if (nodes[other].neighbors(i) == one)
                nodes[other].removeNeighbor(i);
            else
                nodes[nodes[other].neighbors(i)].replaceReferenceTo(other, one);
        
        for (i=nodes[one].degree()-1; i>=0; i--)
            if (nodes[one].neighbors(i) == other)
                nodes[one].removeNeighbor(i);
            
        for (i=0; i<nodes[other].degree(); i++)
            nodes[one].appendNeighbor(nodes[other].neighbors(i));


        // remove double edges
        for (i=nodes[one].degree()-1; i>=0; i--)
            for (j=i-1; j>=0; j--) 
                if (nodes[one].neighbors(i)==nodes[one].neighbors(j)) {
                    for (k=nodes[nodes[one].neighbors(i)].degree()-1; k>=0; k--)
                        if (nodes[nodes[one].neighbors(i)].neighbors(k)==one) {
                            nodes[nodes[one].neighbors(i)].removeNeighbor(k);
                            break;
                        }
                    
                    nodes[one].removeNeighbor(i);
                    break;
                }
        
            invalidate(other);

    }

    /**@name access methods */
    //@{
    ///
    int map(StaticVector<float,2>& domainCoord, std::tr1::array<NodeIdx, 3>& vertices, StaticVector<float,2>& coords,
            int seed=-1) const;
    //@}

    /** \brief Computes the orientation of a point relative to an oriented line
     *
     * \return
     *  -1 : clockwise
     *   0 : collinear
     *   1 : counterclockwise
     */
    signed char orientation(const DirectedEdgeIterator& cE, const StaticVector<float,2>& p) const {

        const StaticVector<float,2>& f = nodes[cE.from()].domainPos();
        const StaticVector<float,2>& t = nodes[cE.to()].domainPos();

        return orientation(f, t, p);
    }

    /** \brief Computes the orientation of three points
     *
     * \return
     *  -1 : clockwise
     *   0 : collinear
     *   1 : counterclockwise
     */
    static signed char orientation(const StaticVector<float,2>& a, const StaticVector<float,2>& b, const StaticVector<float,2>& c) {

        StaticVector<float,2> n = StaticVector<float,2>(a[1] - b[1], b[0] - a[0]);

        float scalarProd = n.dot(c-a);
        if (scalarProd > 0)
            return 1;
        else if (scalarProd < 0)
            return -1;
        else 
            return 0;
    }

    /** \brief Transforms the domain positions into a world coordinate system
     *
     * Assuming that the domain positions of the nodes contained in this
     * PlaneParam are given in barycentric coordinates, this function 
     * transforms them into a world coordinates.  The new coordinate system
     * is specified by supplying new coordinates for the three points
     * (1,0), (0,1), and (0, 0)
     */
    void installWorldCoordinates(const StaticVector<float,2> &a, const StaticVector<float,2> &b, const StaticVector<float,2> &c);

    /** assuming the domain coordinates are given as world coordinates
        this routines turns them into barycentric ones, given the vertices of a bounding triangle. */
    void installBarycentricCoordinates(const StaticVector<float,2> &a, const StaticVector<float,2> &b, const StaticVector<float,2> &c);

    ///
    void applyParametrization(const std::vector<StaticVector<float,3> >& nodePositions);

    ///
    void unflipTriangles(const std::vector<StaticVector<float,3> >& nodePositions);

    ///
    void augmentNeighborIdx(int d) {
        for (size_t i=0; i<nodes.size(); i++)
            for (int j=0; j<nodes[i].degree(); j++)
                nodes[i].neighbors(j) += d;

    }

    ///
    void invalidate(int i) {
        nodes[i].valid = false;
    }

    
    ///
    void makeCyclicBoundaryNode(Node& center, int next, int previous);

    void removeExtraEdges();
    
    void computeFloaterLambdas(SparseMatrix<float>& lambda_ij,
                               const std::vector<StaticVector<float,3> >& nodePositions);

    bool polarMap(const StaticVector<float,3>& center, const std::vector<StaticVector<float,3> > &threeDStarVertices, 
                  std::vector<StaticVector<float,2> >& flattenedCoords, std::vector<float>& theta);

public:
    /**@name debug code */
    //@{
    /// a debug consistency check
    void checkConsistency(const char* where) const;

    /// lists the contents
    void print(bool showNodes=false, bool showParamEdges=false, bool showExtraEdges=false) const;
    //@}

    ///
    static StaticVector<float,2> computeBarycentricCoords(const StaticVector<float,2> &p, const StaticVector<float,2> &a, const StaticVector<float,2> &b, const StaticVector<float,2> &c);

    /** \brief Computes the barycentric coordinates of a point.
     *
     * This routine computes the barycentric coordinates of a point in space 
     * with respect to a triangle in space.  It tacitly assumes that the point 
     * is coplanar with the triangle.
     */
    static StaticVector<float,2> computeBarycentricCoords(const StaticVector<float,3> &p, const StaticVector<float,3> &a, const StaticVector<float,3> &b, const StaticVector<float,3> &c);

    ///
    DirectedEdgeIterator BFLocate(const StaticVector<float,2> &p, int seed=-1) const;

    ///
    void makeCyclicGeometrically(Node& center);

    ///
    void makeCyclicInteriorNode(Node& center);

    ///
    //virtual void makeCyclicBoundaryNode(Node& center);
    
    ///
    bool DFSBoundaryVisit(const std::vector<Node::NeighborReference> &star, 
                          Node::NeighborReference u, int endNode,
                          std::vector<Node::NeighborReference> &outStar);

    ///
    bool DFSVisit(const std::vector<Node::NeighborReference> &star, const Node::NeighborReference& u, 
                  std::vector<Node::NeighborReference> &outStar);

public:
    ///////////////////////////////////////////////////////////////////////
    // These are interpolation methods for your programming convenience
    ///////////////////////////////////////////////////////////////////////

    template<class T>
    static T linearInterpol(const StaticVector<float,2>& p, const T& a, const T& b, const T& c){
        T result = p[0]*a + p[1]*b + (1-p[0]-p[1])*c;
        return result;
    }
        
    template<class T>
    static T linearInterpol(float lambda, const T& a, const T& b){
        T result = a + lambda*(b-a);
        return result;
    }

    //////////////////////////////////////////////////////////////

    ///
    std::vector<Node> nodes;

};

#endif

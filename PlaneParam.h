/////////////////////////////////////////////////////////////////
/*
 * $Id: PlaneParam.h,v 1.2 2007/10/18 16:05:58 sander Exp $
 *
 * $Log: PlaneParam.h,v $
 * Revision 1.2  2007/10/18 16:05:58  sander
 * stuff from 'contact' merged, renamed to psurface
 *
 * Revision 1.1  2007/10/17 13:16:55  sander
 * moved here from the ZIB server
 *
 * Revision 1.30  2006/12/20 09:19:15  bzfsande
 * fix for gcc-4.1: remove extra qualification
 * mailtoauthor: sander@zib.de
 *
 * Revision 1.29  2005/08/23 12:10:59  bzfsande
 * minor cleanup
 * mailtoauthor: sander@zib.de
 *
 * Revision 1.28  2005/02/03 17:27:33  bzfsande
 * documentation
 * mailtoauthor: sander@zib.de
 *
 * Revision 1.27  2003/08/27 12:22:19  bzfsande
 * fixes
 *
 * Revision 1.26  2003/06/30 12:44:44  bzfsande
 * generalizations needed for the contact library
 *
 * Revision 1.25  2003/06/05 13:01:33  bzfsande
 * introduces ghost nodes as a fifth node type
 * also, the access to the node domain positions is procedural now
 *
 * Revision 1.24  2003/05/09 08:58:20  bzfsande
 * bugfixes
 *
 * Revision 1.23  2003/04/04 14:59:18  bzfsande
 * The base grid is now array-based
 *
 * Revision 1.22  2003/03/24 13:26:53  bzfsande
 * separated the Parametrization object into a base grid object connected to a standard Surface
 *
 * Revision 1.21  2002/10/14 12:25:12  bzfsande
 * SUN fixes
 *
 * Revision 1.20  2002/10/08 13:00:29  bzfsande
 * mapping function with an explicit seed
 *
 * Revision 1.19  2002/10/02 15:22:39  bzfsande
 * Introduced the global #imagePos# array
 * - any data defined on the vertices of the original surface
 *   can now be queried
 * - The base grid vertices don't have to be at the same
 *   position as their homologues on the original surface anymore
 *
 */
/////////////////////////////////////////////////////////////////
#ifndef PLANE_PARAM_H
#define PLANE_PARAM_H

#include <mclib/McVec3f.h>
#include <mclib/McVec2f.h>
#include <mclib/McMat3f.h>
#include <mclib/McDArray.h>
#include <mclib/McSArray.h>
#include <mclib/McSmallArray.h>

#include "Node.h"

#include "psurfaceAPI.h"

template<class T, bool SYMMETRIC> class McSparseMatrix;
class Parametrization;

/** The parametrization of a surface in space.  It is actually a planar graph
   having only triangular and unbounded faces.  The nodes of the graph store
   a 3D position, the rest is computed via linear interpolation. */
class PSURFACE_API PlaneParam{
public:

    class DirectedEdgeIterator {
    public:
        DirectedEdgeIterator() : fromNode(-1), neighborIdx(0) , nodes(0){}

        DirectedEdgeIterator(const McDArray<Node>& _nodes) :  fromNode(-1), neighborIdx(0) {
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

        const McDArray<Node>* nodes;

    };

    class UndirectedEdgeIterator {
    public:
        UndirectedEdgeIterator() : fromNode(-1), neighborIdx(0) , nodes(0){}

        UndirectedEdgeIterator(const McDArray<Node>& _nodes) :  fromNode(-1), neighborIdx(0) {
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

        const McDArray<Node>* nodes;

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
        McSArray<NodeIdx, 3> vertices() const {
            McSArray<NodeIdx, 3> result;
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

        edge.neighborIdx = mcSmallArray::index(nodes[from].nbs, to);

        if (edge.neighborIdx == -1)
            edge.fromNode = -1;

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

        nodes[0].setValue(McVec2f(1, 0), a, Node::CORNER_NODE);
        nodes[1].setValue(McVec2f(0, 1), b, Node::CORNER_NODE);
        nodes[2].setValue(McVec2f(0, 0), c, Node::CORNER_NODE);

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
        for (int i=0; i<nodes.size(); i++)
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
    int map(McVec2f& domainCoord, McSArray<NodeIdx, 3>& vertices, McVec2f& coords,
            int seed=-1) const;
    //@}

    /** \brief Computes the orientation of a point relative to an oriented line
     *
     * \return
     *  -1 : clockwise
     *   0 : collinear
     *   1 : counterclockwise
     */
    signed char orientation(const DirectedEdgeIterator& cE, const McVec2f& p) const {

        const McVec2f& f = nodes[cE.from()].domainPos();
        const McVec2f& t = nodes[cE.to()].domainPos();

        return orientation(f, t, p);
    }

    /** \brief Computes the orientation of three points
     *
     * \return
     *  -1 : clockwise
     *   0 : collinear
     *   1 : counterclockwise
     */
    signed char orientation(const McVec2f& a, const McVec2f& b, const McVec2f& c) const {

        McVec2f n = McVec2f(a.y - b.y, b.x - a.x);

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
    void installWorldCoordinates(const McVec2f &a, const McVec2f &b, const McVec2f &c);

    /** assuming the domain coordinates are given as world coordinates
        this routines turns them into barycentric ones, given the vertices of a bounding triangle. */
    void installBarycentricCoordinates(const McVec2f &a, const McVec2f &b, const McVec2f &c);

    ///
    void applyParametrization(const McDArray<McVec3f>& nodePositions);

    ///
    void unflipTriangles(const McDArray<McVec3f>& nodePositions);

    ///
    void augmentNeighborIdx(int d) {
        for (int i=0; i<nodes.size(); i++)
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
    
    void computeFloaterLambdas(McSparseMatrix<float, false>& lambda_ij,
                               const McDArray<McVec3f>& nodePositions);

    bool polarMap(const McVec3f& center, const McSmallArray<McVec3f, 15> &threeDStarVertices, 
                  McSmallArray<McVec2f, 15>& flattenedCoords, McSmallArray<float, 15>& theta);

public:
    /**@name debug code */
    //@{
    /// a debug consistency check
    void checkConsistency(const char* where) const;

    /// lists the contents
    void print(bool showNodes=false, bool showParamEdges=false, bool showExtraEdges=false) const;
    //@}

    ///
    static McVec2f computeBarycentricCoords(const McVec2f &p, const McVec2f &a, const McVec2f &b, const McVec2f &c);

    /** \brief Computes the barycentric coordinates of a point.
     *
     * This routine computes the barycentric coordinates of a point in space 
     * with respect to a triangle in space.  It tacitly assumes that the point 
     * is coplanar with the triangle.
     */
    static McVec2f computeBarycentricCoords(const McVec3f &p, const McVec3f &a, const McVec3f &b, const McVec3f &c);

    ///
    DirectedEdgeIterator BFLocate(const McVec2f &p, int seed=-1) const;

    ///
    void makeCyclicGeometrically(Node& center);

    ///
    void makeCyclicInteriorNode(Node& center);

    ///
    //virtual void makeCyclicBoundaryNode(Node& center);
    
    ///
    bool DFSBoundaryVisit(const McSmallArray<Node::NeighborReference, 6> &star, 
                          const Node::NeighborReference& u, int endNode,
                          McSmallArray<Node::NeighborReference, 6> &outStar);

    ///
    bool DFSVisit(const McSmallArray<Node::NeighborReference, 6> &star, const Node::NeighborReference& u, 
                              McSmallArray<Node::NeighborReference, 6> &outStar);

public:
    ///////////////////////////////////////////////////////////////////////
    // These are interpolation methods for your programming convenience

    template<class T>
    static T linearInterpol(const McVec2f& p, const T& a, const T& b, const T& c){
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
    McDArray<Node> nodes;

};

#endif

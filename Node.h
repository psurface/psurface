#ifndef NODE_H
#define NODE_H

#include <iostream>
#include <algorithm>
#include <stdio.h>
#include "StaticVector.h"
#include <vector>

namespace psurface {

typedef int NodeIdx;


/** \brief The parametrization of a single point from the plane into space 
 *
 * There are five different types of nodes:
 * <ul>
 * <li> <b> Interior Nodes: </b> This is the mapping from a point on the interior of a base 
 *      grid triangle
 *      onto a vertex of the target surface.  Thus it exists once on the base grid, having
 *      unique barycentric coordinates on its triangle.  The target is specified using
 *      the number of the target vertex.
 * <li> <b> Touching Nodes: </b> The mapping from a point on the edge between two base grid
 *      triangles to a vertex on the target surface.
 * <li> <b> Corner Nodes: </b> The mapping from a corner of a base grid
 *      triangle to a vertex on the target surface.
 * <li> <b> Intersection Nodes: </b> When an edge of a target triangle leaves a base triangle
 *      through an edge, then the intersection points are mapped onto target edge points.
 * <li> <b> Ghost Nodes: </b> Each base vertex that is not a Corner node, becomes a ghost node.
 * <li> <b> Boundary Nodes: </b> If an edge leaves the image of the projection, e.g. when the domain boundary end
 *          a boundary node is added at the leaving edge. These nodes store the target vertex index of the point
 *          the edge is pointing to in the member "boundary".
 * </ul>

 \tparam ctype The type used for coordinates
 */
template <class ctype>
class Node {

public:

    /** \brief Encapsulates the neighborhood relationship between a node and others
     *
     * Encapsulates the neighborhood relationship between a node and others.  It's
     * basically just an index into the Node array, but it comes with an extra built-in
     * flag.  This flag specifies whether the graph edge embodied by the reference
     * is 'real', i.e. is (part of) the preimage of some edge in the target surface,
     * or whether it has just been added to turn the graph into a triangulation.
     * Those latter edges form the \em triangular \em closure.
     */
    class NeighborReference {
    public:
        /** \brief Default constructor.  Creates an invalid reference.
         */
        NeighborReference() : idx(-1), closure(false) {}

        
        NeighborReference(int _idx, bool _closure=false) {
            idx = _idx;
            closure = _closure;
        }

        //! Check if edge is real or from triangular closure.
        bool isRegular() const {return !closure;}

        /** \brief Assignment operator from an integer*/
        NeighborReference& operator =(int other) {
            idx = other;
            return *this;
        }

        /** \brief Standard assignment operator. IRIX seems to require this... */
        NeighborReference& operator =(const NeighborReference& other) {
            idx = other.idx;
            closure = other.closure;
            return *this;
        }


        NeighborReference& operator +=(int other) {
            idx += other;
            return *this;
        }

        NeighborReference& operator -=(int other) {
            idx -= other;
            return *this;
        }

        /** \brief Cast to integer */
        operator int() const {return idx;}

        int idx:31;
        bool closure:1;
    };
        
    enum NodeType {INTERIOR_NODE=0, 
                   INTERSECTION_NODE=1, 
                   CORNER_NODE=2, 
                   TOUCHING_NODE=3,
                   GHOST_NODE=4};

    ///
    Node() : valid(true), boundary(-1) {}
        
    /** \brief  Construct node from the local coords, the index and the node type. */
    Node(const StaticVector<ctype,2> &domain, int number, NodeType nodeType) : valid(true), boundary(-1) {
        setDomainPos(domain);
        nodeNumber = number;

        type = nodeType;
            
        clearNeighbors();
    }

        
    /** \brief Destructor. */
    ~Node() {}

    /** \brief Set local coordinates, node index and type. */
    void setValue(const StaticVector<ctype,2> &domain, int nN, NodeType nodeType, int targetIndex = -1) {
        setDomainPos(domain);
        nodeNumber = nN;
            
        type = nodeType;
        boundary = targetIndex;
    }

    void makeInteriorNode() {
        assert(!isINTERSECTION_NODE());
        type = INTERIOR_NODE;
    }

    void makeTouchingNode() {
        assert(!isINTERSECTION_NODE());
        type = TOUCHING_NODE;
    }

    void makeCornerNode() {
        assert(!isINTERSECTION_NODE());
        type = CORNER_NODE;
    }

    /** \brief Set node to be the corner node. */
    void makeCornerNode(int corner, int nN) {
        assert(corner==0 || corner==1 || corner==2);
        type = CORNER_NODE;
        nodeNumber = nN;
        edge = corner;
        if (corner==0)
            dP = StaticVector<ctype,2>(1,0);
        else if (corner==1)
            dP = StaticVector<ctype,2>(0,1);
        else
            dP = StaticVector<ctype,2>(0,0);
    }

    /** \brief Make a ghost node. */
    void makeGhostNode(int corner, int targetTri, const StaticVector<ctype,2>& localTargetCoords) {
        assert(corner==0 || corner==1 || corner==2);
        type = GHOST_NODE;
        nodeNumber = targetTri;
        edge = corner;
        dP = localTargetCoords;
    }

    /** \brief Set the only true neighbor reference if node is an intersection node. */
    void setNeighbor(int n) {
        assert(isINTERSECTION_NODE());
        nbs.resize(1);
        nbs[0] = n;
    }

    /** \brief Is boundary node, i.e. an intersection node on an edge that leaves the image of the projection. */
    bool isBoundary() const {
        return (boundary>=0);
    }

    /** \brief Check if the node is on a segment given by local coords of the endpoints. */
    bool isOnSegment(const StaticVector<ctype,2>& a, const StaticVector<ctype,2>& b, ctype eps) const {
        return ( (domainPos()-a).length() + (domainPos()-b).length()) / (a-b).length() <1.0 +eps;
    }

    /** \brief Check if node is connected to an other node. */
    bool isConnectedTo(const int other) const {
        for (int i=0; i<degree(); i++)
            if (neighbors(i) == other)
                return true;
            
        return false;
    }

    /** \brief The number of neighbors */
    int degree() const {
        return nbs.size(); 
    }

    /** \brief Remove all neighbors */
    void clearNeighbors() {
        nbs.clear();
    }

    /** \brief Erase i'th neighbor reference. */
    void removeNeighbor(int i){
        nbs.erase(nbs.begin() + i);
    }

    /** \brief Add new neighbor. */
    void appendNeighbor(const NeighborReference& newNeighbor){
        nbs.push_back(newNeighbor);
    }
        
    /** \brief Remove the reference to the neighbor which has the index 'other' */
    void removeReferenceTo(int other){
        for (int i=0; i<degree(); i++)
            if (neighbors(i) == other){
                removeNeighbor(i);
                return;
            }
    }
        
    /** \brief Replaces a reference (faster than remove + append). */
    bool replaceReferenceTo(int oldNeighbor, int newNeighbor){
        for (int i=0; i<degree(); i++)
            if (neighbors(i) == oldNeighbor){
                neighbors(i) = newNeighbor;
                return true;
            }

        return false;
    }

    /** \brief Check if node is valid. */
    bool isInvalid() const {return !valid;}

    /** \brief Get i'th neighbor reference. */
    NeighborReference& neighbors(int i) {return nbs[i];}

    /** \brief Get const i'th neighbor reference. */
    const NeighborReference& neighbors(int i) const {return nbs[i];}

    /** \brief Check if node is the first neighbor. */
    bool isFirstNeighbor(NodeIdx n) const {
        return (degree() && neighbors(0) == n);
    }

    /** \brief Check if node is the last neighbor. */
    bool isLastNeighbor(NodeIdx n) const {
        return (degree() && nbs.back() == n);
    }

    /** \brief If the node is an intersection node, then it must have
     *         one neighbor that is not due to the triangluar closure.
     */
    const NeighborReference& theInteriorNode() const {
        assert(isINTERSECTION_NODE());
            
        for (int i=0; i<degree(); i++)
            if (!nbs[i].closure)
                return nbs[i];

        assert(false);
        return nbs[0];
    }
    /** \brief Swap neighbor references. */
    void swapNeighbors(int i, int j) {
        assert(0<=i && i<degree());
        assert(0<=j && j<degree());

        std::swap(nbs[i], nbs[j]);
    }

    /** \brief Reverse the order of the neighbors. */
    void reverseNeighbors() {
        std::reverse(nbs.begin(),nbs.end());
    }
    /** \brief Get node type. */
    NodeType getType() const {return type;}

    /** \brief Get node index. */
    unsigned int getNodeNumber() const {
        return nodeNumber;
    }

    /** \brief If node lives on an edge, return the index of the edge. */
    unsigned int getDomainEdge() const {
        assert(!isINTERIOR_NODE());
        return edge;
    }

    /** \brief Returns a 1-dimensional coordinate \f$ \lambda \in [0,1] \f$ giving
     *         the position of an edge node on its edge.
     *         Is undefined if the node is interior.  Is also undefined if the node
     *         lives on a corner!!!
     */
    ctype getDomainEdgeCoord() const {
        assert(isINTERSECTION_NODE() || isTOUCHING_NODE());
        switch (getDomainEdge()) {
        case 0: return domainPos()[1];
        case 1: return 1-domainPos()[1];
        case 2: return domainPos()[0];
        }
        print();
        assert(false);
    }

    /** \brief Returns a 1-dimensional coordinate \f$ \lambda \in [0,1] \f$ giving
     *         the position of an edge node on its edge.
     *         Is undefined if the node is interior.  Is also undefined if the node
     *         lives on a corner!!!
     */
    ctype getDomainEdgeCoord(int edge) const {
        assert(!isINTERIOR_NODE());
        switch (edge) {
        case 0: return domainPos()[1];
        case 1: return 1-domainPos()[1];
        case 2: return domainPos()[0];
        }
        
        assert(false);
    }

    /** \brief Returns the index of an edge node resp. to its edgePoint array. */
    unsigned int getDomainEdgePosition() const {
        // this is how it should be
        //assert(isINTERSECTION_NODE() || isTOUCHING_NODE());
        return edgePosition;
    }



    /** If the node is of a type which is located at corners, this method
     * returns the corner.  If not, the result is undefined.
     */
    unsigned int getCorner() const {
        assert(isCORNER_NODE() || isGHOST_NODE());
        return edge;
    }

    /** \brief Get the triangle a ghost node points to.
     *
     * Ghost nodes do not point to vertices of the target surface but
     * to a point on a triangle.  Thus it makes sense to ask for that
     * triangle.  Calling this routine for non-ghost nodes results
     * in undefined, possibly fatal, behaviour.
     */
    int getTargetTriangle () const {
        assert(isGHOST_NODE());
        return nodeNumber;
    }
    /** \brief Set the edge index. */
    void setDomainEdge(int i) {
        assert(!isINTERIOR_NODE());
        edge = i;
    }
    /** \brief Set position of the edge node in the edge point array. */
    void setDomainEdgePosition(int i) {
        assert(!isINTERIOR_NODE());
        edgePosition = i;
    }
    /** \brief Check if node is a corner node. */
    bool isCORNER_NODE()       const {return type==CORNER_NODE;}
    /** \brief Check if node is an intersection node. */
    bool isINTERSECTION_NODE() const {return type==INTERSECTION_NODE;}
    /** \brief Check if node is a touching node. */
    bool isTOUCHING_NODE()     const {return type==TOUCHING_NODE;}
    /** \brief Check if node is an interior node. */
    bool isINTERIOR_NODE()     const {return type==INTERIOR_NODE;}
    /** \brief Check if node is a ghost node. */
    bool isGHOST_NODE()        const {return type==GHOST_NODE;}

    /** \brief Checks whether node is on a DomainTriangle edge.
     * 
     * This method returns true if the node is of a type that is by
     * definition located on the edge of the supporting base grid
     * triangle.  <b> Warning: </b> It returns false for nodes
     * on corners!
     */
    bool isOnEdge() const {return isTOUCHING_NODE() || isINTERSECTION_NODE();}

    /** \brief Checks whether node is on a DomainTriangle edge.
     *
     * Unlike the method isOnEdge without arguments this methods
     * also returns true for nodes on corners, given that the corner
     * belongs to the correct edge.
     */
    bool isOnEdge(unsigned int edge) const {
        if (isOnCorner()) 
            return (getCorner()==edge || getCorner()==((edge+1)%3));
        else if (isOnEdge())
            return getDomainEdge()==edge;

        return false;
    }

    /** \brief Check if node is on a corner. */
    bool isOnCorner() const {return isCORNER_NODE() || isGHOST_NODE();}

    /** \brief Print node information. */
    void print(bool showNeighbors=true) const {

        printf("dom (%f %f) ", domainPos()[0], domainPos()[1]);
        
        switch(type){
        case INTERIOR_NODE:
            printf("INTERIOR_NODE");
            break;
        case TOUCHING_NODE:
            printf("TOUCHING_NODE");
            break;
        case INTERSECTION_NODE:
            printf("INTERSECTION_NODE");
            break;
        case CORNER_NODE:
            printf("CORNER_NODE");
            break;
        case GHOST_NODE:
            printf("GHOST_NODE");
            break;
        }
        
        printf(" number %d", nodeNumber);
        printf(" is Boundary %d", boundary);
        
        if (isOnEdge())
            std::cout << "  edge: " << getDomainEdge() << "  edgePos " << getDomainEdgePosition() << std::endl;
        else if (isOnCorner())
            printf("  corner: %d\n", getCorner());
        else
        printf("\n");
        
        if (showNeighbors)
            for (int i=0; i<degree(); i++)
                printf("   %d %s\n", (int)nbs[i], nbs[i].isRegular() ? " " : "c");
        
    }

    /** \brief Get the barycentric coordinates of the node. */
    StaticVector<ctype,2> domainPos() const {
        if (isGHOST_NODE())
            switch (edge) {
            case 0: return StaticVector<ctype,2>(1, 0);
            case 1: return StaticVector<ctype,2>(0, 1);
            case 2: return StaticVector<ctype,2>(0, 0);
            }

        return dP;
    }

    /** \brief Set barycentric coordinates of the node. */
    void setDomainPos(const StaticVector<ctype,2>& p) {dP = p;}

    //private:
    //! Local barycentric coordinates of the node within the triangle.
    StaticVector<ctype,2> dP;

public:    
    //!
    int valid:1;

    //! Type of the node
    NodeType type:3;

    //! Index of target vertex
    unsigned int nodeNumber:28;

    //! This value is set to the target vertex index for boundary nodes and -1 else
    int boundary;

public:
    //! Vector containing all nodes that are connected to this one.
    std::vector<NeighborReference> nbs;
                
    ///////////////////////////////////////
    // This is only for nodes on the boundary of a base grid triangle
                
protected:
    //! If the node is on an edge, the index of the edge is stored
    unsigned int edge:8;
    //! Ordered position on that edge.(0 = 1.Corner, nEdgeNodes-1 = 2.Corner)
    unsigned int edgePosition:24;
                
};

} // namespace psurface

#endif

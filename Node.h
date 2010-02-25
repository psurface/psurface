/////////////////////////////////////////////////////////////////
/*
 * $Id: Node.h,v 1.2 2007/10/18 16:05:58 sander Exp $
 *
 * $Log: Node.h,v $
 * Revision 1.2  2007/10/18 16:05:58  sander
 * stuff from 'contact' merged, renamed to psurface
 *
 * Revision 1.1  2007/10/17 13:16:55  sander
 * moved here from the ZIB server
 *
 * Revision 1.3  2005/08/23 12:10:59  bzfsande
 * minor cleanup
 * mailtoauthor: sander@zib.de
 *
 * Revision 1.2  2004/10/01 13:08:37  bzfsande
 * using mclong instead of unsigned int in order to resolve overloading ambiguity
 * mailtoauthor: sander@zib.de
 *
 * Revision 1.1  2004/01/20 15:01:55  bzfsande
 * the Node class now gets its own file
 *
 */
/////////////////////////////////////////////////////////////////
#ifndef NODE_H
#define NODE_H

#include <iostream>
#include <stdio.h>
#include <psurface/StaticVector.h>
#include <vector>

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
 * <li> <b> Corner Nodes: </b> ddd
 * <li> <b> Intersection Nodes: </b> ddd
 * <li> <b> Ghost Nodes: </b> ddd
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
    Node() : valid(true) {}
        
    ///
    Node(const StaticVector<ctype,2> &domain, int number, NodeType nodeType) : valid(true) {
        setDomainPos(domain);
        nodeNumber = number;

        type = nodeType;
            
        clearNeighbors();
    }

        
    ///
    ~Node() {}

    ///
    void setValue(const StaticVector<ctype,2> &domain, int nN, NodeType nodeType) {
        setDomainPos(domain);
        nodeNumber = nN;
            
        type = nodeType;
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

    void makeGhostNode(int corner, int targetTri, const StaticVector<ctype,2>& localTargetCoords) {
        assert(corner==0 || corner==1 || corner==2);
        type = GHOST_NODE;
        nodeNumber = targetTri;
        edge = corner;
        dP = localTargetCoords;
    }

    void setNeighbor(int n) {
        assert(isINTERSECTION_NODE());
        nbs.resize(1);
        nbs[0] = n;
    }

    ///
    bool isOnSegment(const StaticVector<ctype,2>& a, const StaticVector<ctype,2>& b, ctype eps) const {
        return ( (domainPos()-a).length() + (domainPos()-b).length()) / (a-b).length() <1.0 +eps;
    }

    ///
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

    ///
    void removeNeighbor(int i){
        nbs.erase(nbs.begin() + i);
    }

    ///
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
        
    /// replaces a reference (faster than remove + append)
    bool replaceReferenceTo(int oldNeighbor, int newNeighbor){
        for (int i=0; i<degree(); i++)
            if (neighbors(i) == oldNeighbor){
                neighbors(i) = newNeighbor;
                return true;
            }

        return false;
    }
        
    bool isInvalid() const {return !valid;}

    ///
    NeighborReference& neighbors(int i) {return nbs[i];}

    ///
    const NeighborReference& neighbors(int i) const {return nbs[i];}

    ///
    bool isFirstNeighbor(NodeIdx n) const {
        return (degree() && neighbors(0) == n);
    }

    ///
    bool isLastNeighbor(NodeIdx n) const {
        return (degree() && nbs.back() == n);
    }

    const NeighborReference& theInteriorNode() const {
        assert(isINTERSECTION_NODE());
            
        for (int i=0; i<degree(); i++)
            if (!nbs[i].closure)
                return nbs[i];

        assert(false);
        return nbs[0];
    }

    void swapNeighbors(int i, int j) {
        assert(0<=i && i<degree());
        assert(0<=j && j<degree());

        std::swap(nbs[i], nbs[j]);
    }

    void reverseNeighbors() {
        std::reverse(nbs.begin(),nbs.end());
    }

    NodeType getType() const {return type;}

    unsigned int getNodeNumber() const {
        return nodeNumber;
    }

    unsigned int getDomainEdge() const {
        assert(!isINTERIOR_NODE());
        return edge;
    }

    /** Returns a 1-dimensional coordinate \f$ \lambda \in [0,1] \f$ giving
     * the position on an edge node on its edge.
     * Is undefined if the node is interior.  Is also undefined if the node
     * lives on a corner!!!     
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

    void setDomainEdge(int i) {
        assert(!isINTERIOR_NODE());
        edge = i;
    }

    void setDomainEdgePosition(int i) {
        assert(!isINTERIOR_NODE());
        edgePosition = i;
    }

    bool isCORNER_NODE()       const {return type==CORNER_NODE;}
    bool isINTERSECTION_NODE() const {return type==INTERSECTION_NODE;}
    bool isTOUCHING_NODE()     const {return type==TOUCHING_NODE;}
    bool isINTERIOR_NODE()     const {return type==INTERIOR_NODE;}
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

    ///
    bool isOnCorner() const {return isCORNER_NODE() || isGHOST_NODE();}

    ///
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

    /// query domain position
    StaticVector<ctype,2> domainPos() const {
        if (isGHOST_NODE())
            switch (edge) {
            case 0: return StaticVector<ctype,2>(1, 0);
            case 1: return StaticVector<ctype,2>(0, 1);
            case 2: return StaticVector<ctype,2>(0, 0);
            }

        return dP;
    }

    /// set domain position
    void setDomainPos(const StaticVector<ctype,2>& p) {dP = p;}

    //private:
    StaticVector<ctype,2> dP;

public:    
    ///
    int valid:1;

    ///
    NodeType type:3;

    /// 
    unsigned int nodeNumber:28;

    /// connectivity
public:
    std::vector<NeighborReference> nbs;
                
    ///////////////////////////////////////
    // This is only for nodes on the boundary of a base grid triangle
                
protected:
    unsigned int edge:8;

    unsigned int edgePosition:24;
                
};

#endif

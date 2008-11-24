#ifndef DOMAINS_H
#define DOMAINS_H

#include <mclib/McVec2f.h>

#include "McPointerSurfaceParts.h"
#include "PlaneParam.h"

#include "psurfaceAPI.h"

class DomainTriangle;
class DomainVertex;
class DomainEdge;


/** A vertex of the base domain */
class PSURFACE_API DomainVertex: public McVertex {
public:

    DomainVertex() {}

    DomainVertex(const McVec3f &a) : McVertex(a) {}

    ~DomainVertex() {}

    DomainVertex& operator=(const McVec3f& rhs) {
        x = rhs.x;
        y = rhs.y;
        z = rhs.z;
        return *this;
    }

};


/** A triangle containing a plane triangulation */
class PSURFACE_API DomainTriangle : public McTriangle,
                                           public PlaneParam
{
public:
    /// default constructor
    DomainTriangle()
    {
        edgePoints[0].clear();
        edgePoints[1].clear();
        edgePoints[2].clear();
    }

    /// creates a domain triangle with an empty parametrization
    DomainTriangle(VertexIdx vertexIdx[3]) : McTriangle(vertexIdx), PlaneParam() {}

    /// creates a domain triangle with an empty parametrization
    DomainTriangle(VertexIdx a, VertexIdx b, VertexIdx c) : McTriangle(a, b, c), PlaneParam() {}

    ~DomainTriangle() {}

public:

    /// creates the identical parametrization
    void makeOneTriangle(int a, int b, int c){
        PlaneParam::makeOneTriangle(a, b, c);

        edgePoints[0].resize(2);
        edgePoints[0][0] = 0;
        edgePoints[0][1] = 1;

        edgePoints[1].resize(2);
        edgePoints[1][0] = 1;
        edgePoints[1][1] = 2;

        edgePoints[2].resize(2);
        edgePoints[2][0] = 2;
        edgePoints[2][1] = 0;
    }   

    /**
     * <b> Warning:</b> This routine might not work properly if GHOST_NODEs are present!
     */
    void insertExtraEdges();

    /**
     *
     * \bug Crashes if edgePoint array contains less than two entries.
     */
    void createPointLocationStructure();



    /// inverses orientation
    void flip();

    /** \brief Turns one third
     *
     * \todo The transformation of the node-domainPositions is not efficient!
     */
    void rotate();

    void updateEdgePoints(int oldNode, int newNode);

    ///
    void adjustTouchingNodes();

    void augmentNeighborIdx(int d) {
        PlaneParam::augmentNeighborIdx(d);

        for (int i=0; i<3; i++)
            for (int j=0; j<edgePoints[i].size(); j++)
                edgePoints[i][j] += d;
    }

    ///
    int cornerNode(int i) const {
        //assert(edgePoints[i][0]->isCORNER_NODE());
        return edgePoints[i][0];
    }

    /** \brief Returns the index of an edge node resp. to its edgePoint array. */
    unsigned int getDomainEdgePosition(NodeIdx cN, size_t j) const {
        assert(!nodes[cN].isINTERIOR_NODE());
        if (nodes[cN].isTOUCHING_NODE() || nodes[cN].isINTERSECTION_NODE())
            return nodes[cN].getDomainEdgePosition();

        if (nodes[cN].getCorner()==j)
            return 0;
        else if (nodes[cN].getCorner() == ((j+1)%3))
            return edgePoints[j].size()-1;

        assert(false);
    }

    /** assuming the domain coordinates are given as world coordinates
        this routines turns them into barycentric ones. */
    void installBarycentricCoordinates(){

        const McVec2f a = nodes[cornerNode(0)].domainPos();
        const McVec2f b = nodes[cornerNode(1)].domainPos();
        const McVec2f c = nodes[cornerNode(2)].domainPos();

        PlaneParam::installBarycentricCoordinates(a, b, c);
    }



    /**@name debug code */
    //@{
    /// prints info about the triangle (to stdout)
    void print(bool showEdgePoints=false, bool showParamEdges=false, bool showNodes=false) const; 

    /// checks the triangle for internal consistency
    void checkConsistency(const char* where) const;

    //@}


    /// a list of all nodes that are located exactly on the boundary of the triangle
    McSArray<McSmallArray<NodeIdx, 2>, 3> edgePoints;

    /// the patch number
    int patch;

};




//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////


/** An edge in the base domain */
class PSURFACE_API DomainEdge: public McEdge<DomainVertex> {
public:
    /// default constructor
    DomainEdge() : McEdge<DomainVertex>() {}
    
    ///
    DomainEdge(VertexIdx a, VertexIdx b) : McEdge<DomainVertex>(a, b) {}

    ~DomainEdge() {}

};

#endif

#ifndef DOMAINS_H
#define DOMAINS_H

#include "StaticVector.h"

#include "SurfaceParts.h"
#include "PlaneParam.h"

#include "psurfaceAPI.h"

namespace psurface {

/** A triangle containing a plane triangulation */
template <class ctype>
class PSURFACE_API DomainTriangle : public McTriangle,
                       public PlaneParam<ctype>
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
    DomainTriangle(int vertexIdx[3]) : McTriangle(vertexIdx), PlaneParam<ctype>() {}

    /// creates a domain triangle with an empty parametrization
    DomainTriangle(int a, int b, int c) : McTriangle(a, b, c), PlaneParam<ctype>() {}

public:

    /// creates the identical parametrization
    void makeOneTriangle(int a, int b, int c){
        PlaneParam<ctype>::makeOneTriangle(a, b, c);

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

    /** \brief Cyclically permute a std::vector */
    template <class T>
    static void rotate(std::vector<T>& vec, int offset) {
        int i, s = vec.size();
        T* data = &vec[0];
        if (offset<0) {
            int n = -offset;
            T* tmp = (T*)alloca(n*sizeof(T));
            for (i=0; i<n; i++)
                tmp[i] = data[i];
            for (i=0; i<s-n; i++)
                data[i] = data[i+n];
            for (i=0; i<n; i++)
                data[i+s-n] = tmp[i];
        } else if (offset>0) {
            int n = offset;
            T* tmp = (T*)alloca(n*sizeof(T));
            for (i=0; i<n; i++)
                tmp[i] = data[s-n+i];
            for (i=s-n-1; i>=0; i--)
                data[n+i] = data[i];
            for (i=0; i<n; i++)
                data[i] = tmp[i];
        }
    }

    void updateEdgePoints(int oldNode, int newNode);

    ///
    void adjustTouchingNodes();

    void augmentNeighborIdx(int d) {
        PlaneParam<ctype>::augmentNeighborIdx(d);

        for (int i=0; i<3; i++)
            for (size_t j=0; j<edgePoints[i].size(); j++)
                edgePoints[i][j] += d;
    }

    ///
    int cornerNode(int i) const {
        //assert(edgePoints[i][0]->isCORNER_NODE());
        return edgePoints[i][0];
    }

    /** \brief Returns the index of an edge node resp. to its edgePoint array. */
    unsigned int getDomainEdgePosition(NodeIdx cN, size_t j) const {
        assert(!this->nodes[cN].isINTERIOR_NODE());
        if (this->nodes[cN].isTOUCHING_NODE() || this->nodes[cN].isINTERSECTION_NODE())
            return this->nodes[cN].getDomainEdgePosition();

        if (this->nodes[cN].getCorner()==j)
            return 0;
        else if (this->nodes[cN].getCorner() == ((j+1)%3))
            return edgePoints[j].size()-1;

        assert(false);
    }

    /** assuming the domain coordinates are given as world coordinates
        this routines turns them into barycentric ones. */
    void installBarycentricCoordinates(){

        const StaticVector<ctype,2> a = this->nodes[cornerNode(0)].domainPos();
        const StaticVector<ctype,2> b = this->nodes[cornerNode(1)].domainPos();
        const StaticVector<ctype,2> c = this->nodes[cornerNode(2)].domainPos();

        PlaneParam<ctype>::installBarycentricCoordinates(a, b, c);
    }



    /**@name debug code */
    //@{
    /// prints info about the triangle (to stdout)
    void print(bool showEdgePoints=false, bool showParamEdges=false, bool showNodes=false) const; 

    /// checks the triangle for internal consistency
    void checkConsistency(const char* where) const;

    //@}


    /// a list of all nodes that are located exactly on the boundary of the triangle
    std::tr1::array<std::vector<NodeIdx>, 3> edgePoints;

    /// the patch number
    int patch;

};

} // namespace psurface

#endif

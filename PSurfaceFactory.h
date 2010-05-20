#ifndef PSURFACE_FACTORY_H
#define PSURFACE_FACTORY_H

#include <tr1/array>

template <int dim, class ctype> class PSurface;
class Surface;
template <class ctype, int dim> class StaticVector;

template <int dim, class ctype>
class PSurfaceFactory
{
public:

    PSurfaceFactory(PSurface<dim,ctype>* psurface)
        : psurface_(psurface)
    {}

    void setTargetSurface(Surface* surface);

    void insertVertex(const StaticVector<ctype,dim+1>& position);

    /** \brief Insert a new domain triangle

    \return The index of the newly inserted triangle
    */
    unsigned int insertSimplex(const std::tr1::array<unsigned int, dim+1>& v);

    /**
       \param domainVertex if the normal projection hits a base grid vertex, this is the vertex
    */
    void insertTargetVertexMapping(unsigned int targetVertex, 
                                   unsigned int domainTriangle, 
                                   const StaticVector<ctype,dim>& domainLocalPosition,
                                   NodeBundle& projectedTo,
                                   int& domainVertex);

    void insertGhostNode(unsigned int domainVertex,
                         unsigned int targetTriangle,
                         const StaticVector<ctype,dim>& targetLocalPosition);

    void insertEdge();

    NodeIdx addInteriorNode(int tri, const StaticVector<ctype,2>& dom, int nodeNumber);

    NodeIdx addGhostNode(int tri, int corner, int targetTri, const StaticVector<ctype,2>& localTargetCoords);

    NodeIdx addCornerNode(int tri, int corner, int nodeNumber);

    NodeBundle addIntersectionNodePair(int tri1, int tri2,
                                    const StaticVector<ctype,2>& dP1, const StaticVector<ctype,2>& dP2, 
                                    int edge1, int edge2, const StaticVector<ctype,3>& range);

    NodeIdx addTouchingNode(int tri, const StaticVector<ctype,2>& dP, int edge, int nodeNumber);

    /** \todo Sollte vielleicht ein Bundle zurückgeben */
    NodeIdx addTouchingNodePair(int tri1, int tri2,
                                const StaticVector<ctype,2>& dP1, const StaticVector<ctype,2>& dP2, 
                                int edge1, int edge2, int nodeNumber);

    void addParTriangle(int tri, const std::tr1::array<int,3>& p);
    


protected:

    void addCornerNodeBundle(int domainVertex, int targetVertex);

private:

    PSurface<dim,ctype>* psurface_;

};

#endif

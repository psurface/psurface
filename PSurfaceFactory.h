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

    void insertTargetVertexMapping(unsigned int targetVertex, 
                                   unsigned int domainTriangle, 
                                   const StaticVector<ctype,dim>& domainLocalPosition);

    void insertGhostNode(unsigned int domainVertex,
                         unsigned int targetTriangle,
                         const StaticVector<ctype,dim>& targetLocalPosition);

    void insertEdge();

private:

    PSurface<dim,ctype>* psurface_;

};

#endif

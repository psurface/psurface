#include <vector>

#include <psurface/PSurface.h>
#include <psurface/PSurfaceFactory.h>

template <int dim, class ctype>
void PSurfaceFactory<dim,ctype>::setTargetSurface(Surface* surface)
{
    psurface_->surface = surface;
}

template <int dim, class ctype>
void PSurfaceFactory<dim,ctype>::insertVertex(const StaticVector<ctype,dim+1>& position)
{
    psurface_->newVertex(position);
}

template <int dim, class ctype>
unsigned int  PSurfaceFactory<dim,ctype>::insertSimplex(const std::tr1::array<unsigned int, dim+1>& v)
{
    unsigned int idx = psurface_->createSpaceForTriangle(v[0], v[1], v[2]);
    psurface_->integrateTriangle(idx);
    return idx;
}

template <int dim, class ctype>
void PSurfaceFactory<dim,ctype>::insertTargetVertexMapping(unsigned int targetVertex, 
                                                           unsigned int domainTriangle, 
                                                           const StaticVector<ctype,dim>& domainLocalPosition)
{
}

template <int dim, class ctype>
void PSurfaceFactory<dim,ctype>::insertGhostNode(unsigned int domainVertex,
                                                 unsigned int targetTriangle,
                                                 const StaticVector<ctype,dim>& targetLocalPosition)
{
    std::vector<int> neighbors = psurface_->getTrianglesPerVertex(domainVertex);
    
    for (int i=0; i<neighbors.size(); i++) {
        
        int corner = psurface_->triangles(neighbors[i]).getCorner(domainVertex);
        psurface_->addGhostNode(neighbors[i], corner, targetTriangle, targetLocalPosition);
        
    }
}

template <int dim, class ctype>
void PSurfaceFactory<dim,ctype>::insertEdge()
{
}

// ///////////////////////////////////////////////////////////////////////
//   Explicitly instantiate 'float' and 'double' versions of this code
// ///////////////////////////////////////////////////////////////////////

//template class PSurfaceFactory<1,float>;
//template class PSurfaceFactory<1,double>;

template class PSurfaceFactory<2,float>;
template class PSurfaceFactory<2,double>;

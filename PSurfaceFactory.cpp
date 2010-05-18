#include <psurface/PSurfaceFactory.h>

template <int dim, class ctype>
unsigned int  PSurfaceFactory<dim,ctype>::insertSimplex(const std::tr1::array<unsigned int, dim+1>& v)
{
    assert(dim==2);
    unsigned int idx = psurface_->createSpaceForTriangle(v[0], v[1], v[2]);
    psurface_->integrateTriangle(idx);
    return idx;
}

// ///////////////////////////////////////////////////////////////////////
//   Explicitly instantiate 'float' and 'double' versions of this code
// ///////////////////////////////////////////////////////////////////////

//template class PSurfaceFactory<1,float>;
//template class PSurfaceFactory<1,double>;

template class PSurfaceFactory<2,float>;
template class PSurfaceFactory<2,double>;

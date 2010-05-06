#ifndef PSURFACE_FACTORY_H
#define PSURFACE_FACTORY_H

template <int dim, class ctype>
class PSurfaceFactory
{
public:

    PSurfaceFactory(PSurface<dim,ctype>* psurface)
        : psurface_(psurface)
    {}

    void insertVertex(const StaticVector<ctype,dim+1>& position)
    {
        psurface_->newVertex(position);
    }

    /** \brief Insert a new domain triangle

    \return The index of the newly inserted triangle
    */
    unsigned int insertSimplex(const std::tr1::array<unsigned int, dim+1>& v)
    {
        assert(dim==2);
        unsigned int idx = psurface_->createSpaceForTriangle(v[0], v[1], v[2]);
        psurface_->integrateTriangle(idx);
        return idx;
    }

private:

    PSurface<dim,ctype>* psurface_;

};

#endif

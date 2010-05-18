#ifndef PSURFACE_FACTORY_H
#define PSURFACE_FACTORY_H

template <int dim, class ctype> class PSurface;
class Surface;

template <int dim, class ctype>
class PSurfaceFactory
{
public:

    PSurfaceFactory(PSurface<dim,ctype>* psurface)
        : psurface_(psurface)
    {}

    void setTargetSurface(const Surface* surface);

    void insertVertex(const StaticVector<ctype,dim+1>& position);

    /** \brief Insert a new domain triangle

    \return The index of the newly inserted triangle
    */
    unsigned int insertSimplex(const std::tr1::array<unsigned int, dim+1>& v);

private:

    PSurface<dim,ctype>* psurface_;

};

#endif

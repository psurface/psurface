#ifndef SURFACE_PATH_SET
#define SURFACE_PATH_SET

#include <psurface/SurfacePath.h>


class PARAMETRIZATION_API SurfacePathSet : public McDArray<SurfacePath> {

public:

    void removePoint(VertexIdx p)
    {
        for (int i=size()-1; i>=0; i--)
            (*this)[i].removePoint(p);
    }

};


#endif

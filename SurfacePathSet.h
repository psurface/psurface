#ifndef SURFACE_PATH_SET
#define SURFACE_PATH_SET

#include <vector>

#include <psurface/SurfacePath.h>


class PSURFACE_API SurfacePathSet : public std::vector<SurfacePath> {

public:

    void removePoint(int p)
    {
        for (int i=size()-1; i>=0; i--)
            (*this)[i].removePoint(p);
    }

};


#endif

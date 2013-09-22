#ifndef GMSHIO_H
#define GMSHIO_H
#include <vector>
#include <string>

#include "psurfaceAPI.h"

namespace psurface{

template<class ctype,int dim>
class GmshIO{
    public:
    PSurface<dim,ctype>* par;
    public:
    GmshIO(PSurface<dim,ctype>* psurface);
    ///read psurface_convert from Gmsh file
    void readGmsh(Surface* surf, const std::string&  filename);

    private:
    void    readfile(FILE * file, int cnt, const char* format,
                  void* t1, void* t2, void* t3, void* t4,
                  void* t5 , void* t6, void* , void* t8,
                  void* t9 , void* t10 );
    void    skipline(FILE * file);
};
}
#endif

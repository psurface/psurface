#ifndef GMSHIO_H
#define GMSHIO_H

#include <string>

namespace psurface{

template<class ctype,int dim>
class GmshIO{
    static void    readfile(FILE * file, int cnt, const char* format,
             void* t1, void* t2 = 0, void* t3 = 0, void* t4 = 0,
             void* t5 = 0, void* t6 = 0, void* t7 = 0, void* t8 = 0,
             void* t9 = 0, void* t10 = 0);
    static void    skipline(FILE * file);

    public:
    /// Reads the parametrization of psurface from Gmsh object
    static PSurface<dim, ctype>* readGmsh(const std::string& filename);
};
} // namespace psurface
#endif

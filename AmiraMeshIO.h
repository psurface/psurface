#ifndef AMIRAMESH_IO_H
#define AMIRAMESH_IO_H

#include "psurfaceAPI.h"

// forward declaration for stuff outside of namespace 'psurface'
class AmiraMesh;
#ifdef PSURFACE_STANDALONE
namespace psurface { class Surface; }
#else
class Surface;
#endif

namespace psurface {

// further forward declarations
template <int dim, class ctype> class PSurface;

template <class ctype>
class PSURFACE_API AmiraMeshIO
{
public:

#if defined HAVE_AMIRAMESH || !defined PSURFACE_STANDALONE
    /// Writes the parametrization in AmiraMesh format
    static int writeAmiraMesh(PSurface<2,ctype>* psurface, const char* filename);

    /** \brief Reads the parametrization from an AmiraMesh object
        \todo The return value should be PSurface<dim,ctype>*
    */
    static void* readAmiraMesh(AmiraMesh* am, const char* filename);

    /// AmiraMesh Reader using an existing AmiraMesh object (is used by derived classes)
    static bool initFromAmiraMesh(PSurface<2,ctype>* psurface, AmiraMesh* am, const char* filename, Surface* surf);
#endif

};

} // namespace psurface

#endif

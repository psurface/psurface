#ifndef PSURFACE_SMOOTHER_H
#define PSURFACE_SMOOTHER_H

#include "HxParamToolBox.h"
#include "DomainPolygon.h"

#include "psurfaceAPI.h"

namespace psurface {

template <int dim, class ctype> class PSurface;

template <class ctype>
class PSURFACE_API PSurfaceSmoother
{
    
public:

    /** \brief Smooth the parametrization graph across a specific base grid edge
     *
     * \param psurface The psurface object to be smoothed
     * \param edge The edge to smooth
     * \param keepPatch If this is true, then smoothing on patch boundaries will only be tangential
     * \param nodeStack Source of temporary memory.  Needs to be given from the outside, to be
     *          persistent across calls
     */
    static void applyEdgeRelaxation(PSurface<2,ctype>* psurface, int edge, 
                                    bool keepPatches, std::vector<unsigned int>& nodeStack);

    ///
    void applyVertexRelaxation();

    ///
    static void applyHorizontalRelaxation(DomainPolygon& quadri, PSurface<2,ctype>* psurface);
    
private:

    static void moveSubGraph(int startingNode, DomainPolygon& from, int centerNode);

};

}

#endif

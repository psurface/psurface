#ifndef PSURFACE_SMOOTHER_H
#define PSURFACE_SMOOTHER_H

#include <hxpsurface/HxParamToolBox.h>

template <int dim, class ctype> class PSurface;

template <class ctype>
class PSurfaceSmoother
{
    
public:

    ///
    static void applyEdgeRelaxation(PSurface<2,ctype>* psurface, int edge, 
                                    bool keepPatches, std::vector<unsigned int>& nodeStack);

    ///
    void applyVertexRelaxation();

    ///
    void applyHorizontalRelaxation(DomainPolygon& quadri, PSurface<2,ctype>* psurface);
    
private:

    static void moveSubGraph(int startingNode, DomainPolygon& from, int centerNode);

};

#endif

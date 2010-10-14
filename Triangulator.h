#ifndef MY_TRIANGULATOR
#define MY_TRIANGULATOR

#include <vector>

#include "Domains.h"
#include "MultiDimOctree.h"

#include "EdgeIntersectionFunctor.h"
#include "VertexHeap.h"

template <class ctype> class CircularPatch;
template <int dim, class ctype> class PSurface;
struct QualityRequest;

extern int debugCounter;
/** This class encapsulates algorithms that triangulates small circular
    and semi-circular patches. */
namespace Triangulator {
    
    /// Return orientation (-1 -> clockwise, 0 -> collinear, 1 -> counterclockwise).
    signed char orientation(const StaticVector<float,2>& a, const StaticVector<float,2>& b, const StaticVector<float,2>& c, const float eps=0.0);

    /// performs a flattening described in MAPS (SIGGRAPH 98) and a constrained Delaunay triangulation.
    void triangulateStar(const std::vector<int> &border, int center, 
                                          CircularPatch<float>& resultPatch, std::vector<StaticVector<float,2> >& flatBorder,
                         PSurface<2,float>* par);

    ///
    void triangulateHalfStar(const std::vector<int> &border, 
                                              int center, 
                                              CircularPatch<float>& resultPatch, 
                                              std::vector<StaticVector<float,2> >& flatBorder,
                             PSurface<2,float>* par);

    /// performs a flattening described in MAPS (SIGGRAPH 98) and a constrained Delaunay triangulation.
    void estimateStarError(const std::vector<int>& border, int center, 
                                            const QualityRequest &quality, 
                                            const std::vector<int> &fullStar, 
                                            VertexHeap::ErrorValue& qualityValue,
                                            MultiDimOctree<McEdge, EdgeIntersectionFunctor, float, 3, true>& edgeOctree, 
                           PSurface<2,float>* par); 

    ///
    void estimateHalfStarError(const std::vector<int> &border, int center, 
                                                const QualityRequest &quality,
                                                const std::vector<int> &fullStar, 
                                                VertexHeap::ErrorValue& qualityValue,
                                                MultiDimOctree<McEdge, EdgeIntersectionFunctor, float, 3, true>& edgeOctree, 
                               PSurface<2,float>* par); 


    void planeCDT(const std::vector<StaticVector<float,2> >& flatBorder, const std::vector<int>& border,
                  CircularPatch<float>& result, PSurface<2,float>* par);

    bool isLegalEdge(const StaticVector<float,2>& a, const StaticVector<float,2>& b, const StaticVector<float,2>& c, 
                                      const std::vector<StaticVector<float,2> > &polygon);


    inline float computeAspectRatio(const StaticVector<float,3>& x1, const StaticVector<float,3>& x2, const StaticVector<float,3>& x3) {
        const float a = (x1 - x2).length();
        const float b = (x2 - x3).length();
        const float c = (x3 - x1).length();

        const float aR = 2*a*b*c/((-a+b+c)*(a-b+c)*(a+b-c));
        
        return fabs(aR);
    }

    /// makes an overall evaluation of a CircularPatch using the given QualityRequest object
    void evaluate(const CircularPatch<float>* cP, int removedVertex, 
                                   const QualityRequest &quality, 
                                   VertexHeap::ErrorValue& qualityValue, 
                                   const std::vector<int> &fullStar, 
                                   MultiDimOctree<McEdge, EdgeIntersectionFunctor, float, 3, true>& edgeOctree, 
                  const PSurface<2,float>* par);

};

#endif



#ifndef PARAM_TOOLBOX_H
#define PARAM_TOOLBOX_H

#include <vector>

#include "StaticVector.h"
#include "MultiDimOctree.h"
#include "PlaneParam.h"
#include "Domains.h"

#include "EdgeIntersectionFunctor.h"

#include "psurfaceAPI.h"

namespace psurface {

class DomainPolygon;

template <int dim, class ctype> class PSurface;
struct QualityRequest;

namespace ParamToolBox {

    enum PointType {FEATURE_POINT=-1, REGULAR_POINT=0};
    
    /// Return orientation (-1 -> clockwise, 0 -> collinear, 1 -> counterclockwise).
    signed char PSURFACE_API orientation(const StaticVector<float,2>& a, const StaticVector<float,2>& b, const StaticVector<float,2>& c, const float eps=0.0);

    ///
    int PSURFACE_API computeFeatureStatus(const PSurface<2,float>* par, int v, int& featureEdgeA, 
                                          int& featureEdgeB);

    ///
    bool PSURFACE_API findAllHalfStars(int centerIdx, 
                                           int featureEdgeA, int featureEdgeB, 
                                           std::vector<std::vector<int> > &halfStarVertices, 
                                           std::vector<std::vector<int> > &halfStarTris,
                          std::vector<int>& patches, const PSurface<2,float>* par);
    
    ///
    void PSURFACE_API makeFullStarOutOfHalfStars(std::vector<int>& halfStarVerticesA, 
                                                     std::vector<int>& halfStarTrisA, 
                                                     std::vector<int>& halfStarVerticesB, 
                                                     std::vector<int>& halfStarTrisB, 
                                                     std::vector<int>& fullStarVertices, 
                                                     std::vector<int>& fullStarTris);
    

    
    /// \todo The const_casts used in the routine should be  unnessecary
    void PSURFACE_API mergeTwoTrianglesIntoQuadrangle(int tri1, int tri2,
                                                          DomainPolygon& quadri, bool& flipped, 
                                                          StaticVector<float,2> quadriCoords[4],
                                                          std::vector<unsigned int>& nodeStack,
                                         const PSurface<2,float>* par);
    
    ///
    bool PSURFACE_API mergeStarIntoPolygon(int centerIdx, DomainPolygon& fullStar,
                                               std::vector<int>& fullStarTris, int& newCenterNode,
                                               std::vector<unsigned int>& nodeStack,
                              PSurface<2,float>* par);
    

    
    ///
    void PSURFACE_API flattenStar(int center, const std::vector<int> &threeDStarVertices,
                     std::vector<StaticVector<float,2> >& twoDVertexPos, const PSurface<2,float>* par);
    
    ///
    void PSURFACE_API flattenHalfStar(int center, const std::vector<int>& threeDStarVertices,
                         std::vector<StaticVector<float,2> >& twoDVertexPos, const PSurface<2,float>* par);

    ///
    void PSURFACE_API pizzaCutter(DomainPolygon& fullStar, int newCenterNode,
                                      int& newVertex, std::vector<int>& newTriangles);

    ///
    void PSURFACE_API moveSubGraph(int startingNode, 
                                       DomainPolygon& from, 
                                       std::vector<int>& nodeLocs,
                                       int centerNode);

    ///
    bool PSURFACE_API removeRegularPoint(PSurface<2,float>* par, 
                                             int centerPoint, 
                                             const QualityRequest &quality,
                                             MultiDimOctree<Edge, EdgeIntersectionFunctor, float, 3>* edgeOctree
                                             );

    ///
    bool PSURFACE_API removeFeatureLinePoint(PSurface<2,float>* par, 
                                                 int centerPoint, 
                                                 const QualityRequest &quality,
                                                 int numHalfStars,
                                                 int featureEdgeA,
                                                 int featureEdgeB,
                                                 MultiDimOctree<Edge, EdgeIntersectionFunctor, float, 3>* edgeOctree
                                                 );

    ///
    void PSURFACE_API convexify(std::vector<StaticVector<float,2> >& Coords);

    ///
    void PSURFACE_API convexifyHalfStar(std::vector<StaticVector<float,2> >& coords);

    ///
    bool PSURFACE_API singleTetrahedronTest(const PSurface<2,float>* par, 
                                                const std::vector<int>& fullStarVertices);

    ////////////////////////////////////////////////////////////
    // debug stuff

    ///
    int PSURFACE_API hasIntersections(const PSurface<2,float>* par);

    ///
    float PSURFACE_API longestEdge(const PSurface<2,float>* par); 

    ///
    void PSURFACE_API checkEdge(const PSurface<2,float>* par);
    
};

} // namespace psurface

#endif

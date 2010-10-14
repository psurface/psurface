#ifndef PARAM_TOOLBOX_H
#define PARAM_TOOLBOX_H

#include <vector>

#include "StaticVector.h"
#include "MultiDimOctree.h"
#include "PlaneParam.h"
#include "Domains.h"

#include "EdgeIntersectionFunctor.h"

class DomainPolygon;

template <int dim, class ctype> class PSurface;
class QualityRequest;

namespace ParamToolBox {

    enum PointType {FEATURE_POINT=-1, REGULAR_POINT=0};
    
    /// Return orientation (-1 -> clockwise, 0 -> collinear, 1 -> counterclockwise).
    signed char orientation(const StaticVector<float,2>& a, const StaticVector<float,2>& b, const StaticVector<float,2>& c, const float eps=0.0);

    ///
    int computeFeatureStatus(const PSurface<2,float>* par, int v, int& featureEdgeA, 
                                              int& featureEdgeB, 
                                              const std::vector<bool>* featureVertices=NULL, 
                                              const std::vector<int>* featureEdges=NULL
                                              );

    ///
    bool findAllHalfStars(int centerIdx, 
                                           int featureEdgeA, int featureEdgeB, 
                                           std::vector<std::vector<int> > &halfStarVertices, 
                                           std::vector<std::vector<int> > &halfStarTris,
                          std::vector<int>& patches, const PSurface<2,float>* par);
    
    ///
    void makeFullStarOutOfHalfStars(std::vector<int>& halfStarVerticesA, 
                                                     std::vector<int>& halfStarTrisA, 
                                                     std::vector<int>& halfStarVerticesB, 
                                                     std::vector<int>& halfStarTrisB, 
                                                     std::vector<int>& fullStarVertices, 
                                                     std::vector<int>& fullStarTris);
    

    
    /// \todo The const_casts used in the routine should be  unnessecary
    void mergeTwoTrianglesIntoQuadrangle(int tri1, int tri2,
                                                          DomainPolygon& quadri, bool& flipped, 
                                                          StaticVector<float,2> quadriCoords[4],
                                                          std::vector<unsigned int>& nodeStack,
                                         const PSurface<2,float>* par);
    
    ///
    bool mergeStarIntoPolygon(int centerIdx, DomainPolygon& fullStar,
                                               std::vector<int>& fullStarTris, int& newCenterNode,
                                               std::vector<unsigned int>& nodeStack,
                              PSurface<2,float>* par);
    

    
    ///
    void flattenStar(int center, const std::vector<int> &threeDStarVertices,
                     std::vector<StaticVector<float,2> >& twoDVertexPos, const PSurface<2,float>* par);
    
    ///
    void flattenHalfStar(int center, const std::vector<int>& threeDStarVertices,
                         std::vector<StaticVector<float,2> >& twoDVertexPos, const PSurface<2,float>* par);

    ///
    void pizzaCutter(DomainPolygon& fullStar, int newCenterNode,
                                      int& newVertex, std::vector<int>& newTriangles);

    ///
    void moveSubGraph(int startingNode, 
                                       DomainPolygon& from, 
                                       std::vector<int>& nodeLocs,
                                       int centerNode);

    ///
    bool removeRegularPoint(PSurface<2,float>* par, 
                                             int centerPoint, 
                                             const QualityRequest &quality,
                                             MultiDimOctree<McEdge, EdgeIntersectionFunctor, float, 3, true>* edgeOctree
                                             );

    ///
    bool removeFeatureLinePoint(PSurface<2,float>* par, 
                                                 int centerPoint, 
                                                 const QualityRequest &quality,
                                                 int numHalfStars,
                                                 int featureEdgeA,
                                                 int featureEdgeB,
                                                 MultiDimOctree<McEdge, EdgeIntersectionFunctor, float, 3, true>* edgeOctree,
                                                 std::vector<int>* featureEdges
                                                 );

    ///
    void convexify(std::vector<StaticVector<float,2> >& Coords);

    ///
    void convexifyHalfStar(std::vector<StaticVector<float,2> >& coords);

    ///
    bool singleTetrahedronTest(const PSurface<2,float>* par, 
                                                const std::vector<int>& fullStarVertices);

    ////////////////////////////////////////////////////////////
    // debug stuff

    ///
    void display(const DomainPolygon& pol, int idx);
    
    ///
    void display(const DomainTriangle<float>& tri, int idx);

    ///
    int hasIntersections(const PSurface<2,float>* par);

    ///
    float longestEdge(const PSurface<2,float>* par); 

    ///
    void checkEdge(const PSurface<2,float>* par);
    
};


#endif

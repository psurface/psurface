#include "Triangulator.h"
#include "PSurface.h"
#include "QualityRequest.h"
#include "HxParamToolBox.h"
#include "CircularPatch.h"



#ifdef _WIN32
inline int isnan(double x) {return _isnan(x);}
#endif

using namespace psurface;

signed char Triangulator::orientation(const StaticVector<float,2>& a, const StaticVector<float,2>& b, const StaticVector<float,2>& c, const float eps)
{

    float det = a[0] * (b[1]-c[1]) - b[0] * (a[1] - c[1]) + c[0] * (a[1] - b[1]);

    if (det>eps)
        return 1;
    else if (det<-eps)
        return -1;

    return 0;
    
}

void Triangulator::triangulateStar(const std::vector<int> &border, int center, 
                                   CircularPatch<float>& resultPatch, 
                                   std::vector<StaticVector<float,2> >& flatBorder,
                                   PSurface<2,float>* par)
{
    /////////////////////////////////////
    // computes the flattened coordinates
    ParamToolBox::flattenStar(center, border, flatBorder, par);

    for (int j=0; j<flatBorder.size(); j++)
        assert(!isnan(flatBorder[j][0]) &&!isnan(flatBorder[j][1]));

    ///////////////////////////////////////////
    // do a constrained Delaunay triangulation
    planeCDT(flatBorder, border, resultPatch, par);

}

void Triangulator::estimateStarError(const std::vector<int> &border, int center, 
                                     const QualityRequest &quality, const std::vector<int> &fullStar, 
                                     VertexHeap::ErrorValue& qualityValue,
                                     MultiDimOctree<Edge, EdgeIntersectionFunctor, float, 3>& edgeOctree, 
                                     PSurface<2,float>* par)
{
    /////////////////////////////////////
    // computes the flattened coordinates
    std::vector<StaticVector<float,2> > flatBorder;
    
    ParamToolBox::flattenStar(center, border, flatBorder, par);
    
    for (int j=0; j<flatBorder.size(); j++)
        assert(!isnan(flatBorder[j][0]) &&!isnan(flatBorder[j][1]));
    
    ///////////////////////////////////////////
    // do a constrained Delaunay triangulation
    CircularPatch<float> resultPatch(border.size()-2, par);
    
    planeCDT(flatBorder, border, resultPatch, par);
    
    //////////////////////////////////////////
    // evaluate triangulation, 
    evaluate(&resultPatch, center, quality, qualityValue, fullStar, edgeOctree, par);
    
    resultPatch.killAll();
    
}

// same thing for a half star
void Triangulator::triangulateHalfStar(const std::vector<int> &border, int center, 
                                       CircularPatch<float>& resultPatch, std::vector<StaticVector<float,2> >& flatBorder,
                                       PSurface<2,float>* par)
{
    /////////////////////////////////////
    // computes the flattened coordinates
    ParamToolBox::flattenHalfStar(center, border, flatBorder, par);


    ///////////////////////////////////////////
    // do a constrained Delaunay triangulation

    planeCDT(flatBorder, border, resultPatch, par);
}


// same thing for a half star
void Triangulator::estimateHalfStarError(const std::vector<int> &border, int center, 
                                         const QualityRequest &quality, const std::vector<int> &fullStar, 
                                         VertexHeap::ErrorValue& qualityValue,
                                         MultiDimOctree<Edge, EdgeIntersectionFunctor, float, 3>& edgeOctree, 
                                         PSurface<2,float>* par)
{
    /////////////////////////////////////
    // computes the flattened coordinates

    std::vector<StaticVector<float,2> > flatBorder;
    ParamToolBox::flattenHalfStar(center, border, flatBorder, par);


    ///////////////////////////////////////////
    // do a constrained Delaunay triangulation

    CircularPatch<float> resultPatch(border.size()-2, par);
    planeCDT(flatBorder, border, resultPatch, par);


    //////////////////////////////////////////
    // evaluate triangulation
    evaluate(&resultPatch, center, quality, qualityValue, fullStar, edgeOctree, par);

    resultPatch.killAll();
                
}



////////////////////////////////////////////////////////////////////////////
// the plane Delaunay triangulation of a plane star-shaped polygon
void Triangulator::planeCDT(const std::vector<StaticVector<float,2> >& flatBorder, const std::vector<int>& border,
                            CircularPatch<float>& result, PSurface<2,float>* par)
{

    int K = border.size();

    if (K==3){
        result[0] = par->createSpaceForTriangle(border[0], border[1], border[2]);
        return;
    }

    const signed char counterclockwise = 1;

    int k = 0;
    int idx = 0;
    int edgeIdx = 0;

    std::vector<int> tmpVertices   = border;
    std::vector<StaticVector<float,2> > tmpCoords   = flatBorder;

    while (K>4){

        int bestEdge = -1;
        int bestDelaunayEdge = -1;

        float bestAspectRatio = std::numeric_limits<float>::max();
        float bestDelaunayAspectRatio = std::numeric_limits<float>::max();

        ////////////////////////////////////////
        // look for best edge
        for (k=0; k<K; k++) {

            if (orientation(tmpCoords[k%K], tmpCoords[(k+1)%K], tmpCoords[(k+2)%K])==counterclockwise) {

                const float aspectRatio = computeAspectRatio(par->vertices(tmpVertices[k%K]), 
                                                             par->vertices(tmpVertices[(k+1)%K]), 
                                                             par->vertices(tmpVertices[(k+2)%K]));

                if (aspectRatio<bestAspectRatio) {
                    bestAspectRatio = aspectRatio;
                    bestEdge = k;
                }

                if (isLegalEdge(tmpCoords[k%K], tmpCoords[(k+1)%K], tmpCoords[(k+2)%K], flatBorder) &&
                    aspectRatio<bestDelaunayAspectRatio) {
                    bestDelaunayAspectRatio = aspectRatio;
                    bestDelaunayEdge = k;
                }
            }
        }
        
        assert(bestEdge!=-1);
   

        // /////////////////////////////////////
        // cut off that edge
        
        const int chosenEdge = (bestDelaunayEdge!=-1) ? bestDelaunayEdge : bestEdge;
        
        result[idx++] = par->createSpaceForTriangle(tmpVertices[chosenEdge%K],
                                                    tmpVertices[(chosenEdge+1)%K], 
                                                    tmpVertices[(chosenEdge+2)%K]);
        
        result.innerEdges[edgeIdx][0] = tmpVertices[chosenEdge%K];
        result.innerEdges[edgeIdx][1] = tmpVertices[(chosenEdge+2)%K];
        edgeIdx++;
        tmpVertices.erase(tmpVertices.begin()+((chosenEdge+1)%K));
        tmpCoords.erase(tmpCoords.begin()+((chosenEdge+1)%K)); 
        K--;

    }

    // /////////////////////////////////////
    // the case K=4
    // choose the triangulation which yields the lowest max aspect ratio
    const float aRa1 = computeAspectRatio(par->vertices(tmpVertices[0]), 
                                          par->vertices(tmpVertices[1]), 
                                          par->vertices(tmpVertices[2]));
    const float aRa2 = computeAspectRatio(par->vertices(tmpVertices[2]), 
                                          par->vertices(tmpVertices[3]), 
                                          par->vertices(tmpVertices[0]));

    const float aRb1 = computeAspectRatio(par->vertices(tmpVertices[1]), 
                                          par->vertices(tmpVertices[2]), 
                                          par->vertices(tmpVertices[3]));
    const float aRb2 = computeAspectRatio(par->vertices(tmpVertices[3]), 
                                          par->vertices(tmpVertices[0]), 
                                          par->vertices(tmpVertices[1]));

    const float maxA = (aRa1>aRa2) ? aRa1 : aRa2;
    const float maxB = (aRb1>aRb2) ? aRb1 : aRb2;

    if (maxA<maxB)  {
        result[idx++] = par->createSpaceForTriangle(tmpVertices[0], tmpVertices[1], tmpVertices[2]);
        result[idx++] = par->createSpaceForTriangle(tmpVertices[2], tmpVertices[3], tmpVertices[0]);
        result.innerEdges[edgeIdx][0] = tmpVertices[0];
        result.innerEdges[edgeIdx][1] = tmpVertices[2];
    } else {
        result[idx++] = par->createSpaceForTriangle(tmpVertices[1], tmpVertices[2], tmpVertices[3]);
        result[idx++] = par->createSpaceForTriangle(tmpVertices[3], tmpVertices[0], tmpVertices[1]);
        result.innerEdges[edgeIdx][0] = tmpVertices[1];
        result.innerEdges[edgeIdx][1] = tmpVertices[3];
    }

}


bool Triangulator::isLegalEdge(const StaticVector<float,2>& a, const StaticVector<float,2>& b, const StaticVector<float,2>& c,
                               const std::vector<StaticVector<float,2> > &polygon)
{
    //compute the circumcirle of the points a, b, c
    // code taken from 'GraphicsGems I'
    double d1 = (c - a).dot(b - a);
    double d2 = (c - b).dot(a - b);
    double d3 = (a - c).dot(b - c);
    double c1 = d2*d3, c2 = d3*d1, c3 = d1*d2;
    double c123 = c1 + c2 + c3;

    // test for degeneracy

    if (c123==0) return false;

    float radius = 0.5*sqrt((d1+d2) * (d2+d3) * (d3+d1)/c123);

    if (isnan(radius))
        return false;

    StaticVector<float,2> center = ((c2+c3)*a + (c3+c1)*b + (c1+c2)*c)/(2*c123);

    for (int i=0; i<polygon.size(); i++)
        if (polygon[i]!=a && polygon[i]!=b && polygon[i]!=c &&
            (polygon[i]-center).length()<radius)
            return false;
 
    return true;
}


void Triangulator::evaluate(const CircularPatch<float>* cP, int removedVertex, 
                            const QualityRequest &quality, VertexHeap::ErrorValue& error,
                            const std::vector<int> &fullStar, 
                            MultiDimOctree<Edge, EdgeIntersectionFunctor, float, 3>& edgeOctree, 
                            const PSurface<2,float>* par)
{
    error.unblock();
    
    for (int i=0; i<cP->size(); i++)
        for (int j=0; j<3; j++)
            assert( par->triangles((*cP)[i]).vertices[j] != -1);

    std::vector<int> closeEdges(0);

    if (quality.intersections){

        // get close edges that might intersect the new patch
        Box<float,3> resultBox;
        cP->getBoundingBox(resultBox);

        std::vector<Edge*> tmpCloseEdges;
        edgeOctree.lookup(resultBox, tmpCloseEdges);

        for (size_t i=0; i<tmpCloseEdges.size(); i++) {

            if (tmpCloseEdges[i]->isConnectedTo(removedVertex))
                continue;

            int newEdge = tmpCloseEdges[i] - &par->edges(0);

            if (std::find(closeEdges.begin(), closeEdges.end(), newEdge)==closeEdges.end())
                closeEdges.push_back(newEdge);

        }

        if (cP->intersectsParametrization(closeEdges) || cP->hasSelfintersections()){
            error.block();
            return;
        }

    }

    if (quality.maxEdgeLength>=0) {

        for (size_t i=0; i<cP->innerEdges.size(); i++)
            if ( (par->vertices(cP->innerEdges[i][0]) - par->vertices(cP->innerEdges[i][1])).length() > quality.maxEdgeLength){
                error.block();
                return;
            }

    }
    
    if (cP->inducesTopologyChange()){
        //printf("Induces TopChange\n");
        error.block();
        return;
    }
    
    // compute the Hausdorff distance
    float HausdorffDistance=0;
    if (quality.hausdorffDistance > 0.01) {
        //printf("ev 13\n");
        int nNodes = 0;
        
        for (int i=0; i<fullStar.size(); i++){

            const DomainTriangle<float>& cT = par->triangles(fullStar[i]);
            
            for (size_t cN=0; cN<cT.nodes.size(); cN++) {
                //printf("ev 13.3\n");
                if (cT.nodes[cN].isINTERIOR_NODE() || 
                    cT.nodes[cN].isTOUCHING_NODE()){
                    //printf("ev 13.4\n");
                    nNodes++;
                    
                    HausdorffDistance += cP->distanceTo(par->imagePos(fullStar[i], cN));
                    //printf("ev 13.5\n");
                }
            }
        }
        
        HausdorffDistance += cP->distanceTo(par->vertices(removedVertex));
        HausdorffDistance /= nNodes+1;
    } 
    //printf("ev 14\n");
    float aspectRatioImprovement =  0;

    if (quality.aspectRatio > 0.01) {
        //printf("ev 15\n");
        // compute max aspect ratio of the star at the vertex
        float oldMaxAspectRatio = 0;
        for (int i=0; i<fullStar.size(); i++){
            const float thisAspectRatio = par->aspectRatio(fullStar[i]);
            if (thisAspectRatio>oldMaxAspectRatio)
                oldMaxAspectRatio = thisAspectRatio;
        }
            
        // compute max aspect ratio of retriangulation
        const float newMaxAspectRatio = cP->maxAspectRatio();
        aspectRatioImprovement = newMaxAspectRatio - oldMaxAspectRatio;
    }
    
    error.unblock();
    error.value   = HausdorffDistance*quality.hausdorffDistance + (aspectRatioImprovement)*quality.aspectRatio;
    //printf("ev 16\n");
    return;
}

       


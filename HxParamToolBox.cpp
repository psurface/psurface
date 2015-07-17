#include "config.h"

#include <vector>

#include "HxParamToolBox.h"
#include "Domains.h"
#include "DomainPolygon.h"
#include "CircularPatch.h"
#include "Triangulator.h"
#include "PSurface.h"
#include "QualityRequest.h"

using namespace psurface;

signed char ParamToolBox::orientation(const StaticVector<float,2>& a, const StaticVector<float,2>& b, const StaticVector<float,2>& c, const float eps)
{

    float det = a[0] * (b[1]-c[1]) - b[0] * (a[1] - c[1]) + c[0] * (a[1] - b[1]);

    if (det>eps)
        return 1;
    else if (det<-eps)
        return -1;

    return 0;
}

int ParamToolBox::computeFeatureStatus(const PSurface<2,float>* par, int v,
                                       int& featureEdgeA, int& featureEdgeB)
{
    std::vector<int> featureEdges(0);

    for (int i=0; i<par->vertices(v).degree(); i++){
        const std::vector<int>& cFan = par->edges(par->vertices(v).edges[i]).triangles;

        if (cFan.size() != 2 || par->triangles(cFan[0]).patch != par->triangles(cFan[1]).patch)
            featureEdges.push_back(par->vertices(v).edges[i]);
    }

    if (featureEdges.size() == 0)
        return ParamToolBox::REGULAR_POINT;

    if (featureEdges.size() == 1 || featureEdges.size() > 2)
        return ParamToolBox::FEATURE_POINT;

    // #point# is connected to exactly two feature edges

    if (par->edges(featureEdges[0]).triangles.size() != par->edges(featureEdges[1]).triangles.size())
        return ParamToolBox::FEATURE_POINT;

    featureEdgeA = featureEdges[0];
    featureEdgeB = featureEdges[1];

    return par->edges(featureEdgeA).triangles.size();
}

//////////////////////////////////////////////////////////////////////////////////
// This routines looks for all the half stars that meet at a vertex.
//////////////////////////////////////////////////////////////////////////////////

bool ParamToolBox::findAllHalfStars(int centerIdx,
                                    int featureEdgeA, int featureEdgeB,
                                    std::vector<std::vector<int> > &halfStarVertices,
                                    std::vector<std::vector<int> > &halfStarTris,
                                    std::vector<int> &patches,
                                    const PSurface<2,float>* par)
{
    int i, j;

    int nextCurrentEdge;

    const std::vector<int> &halfStarStarts = par->edges(featureEdgeA).triangles;

    const int numHalfStars = halfStarStarts.size();

    halfStarVertices.resize(numHalfStars);
    halfStarTris.resize(numHalfStars);
    patches.resize(numHalfStars);

    // loop over all halfstars:
    for (i=0; i<numHalfStars; i++){

        halfStarVertices[i].clear();
        halfStarTris[i].clear();

        int currentEdge = featureEdgeA;
        int       currentTri  = halfStarStarts[i];
        patches[i] = par->triangles(currentTri).patch;

        while (currentEdge != featureEdgeB){

            assert(par->triangles(currentTri).isConnectedTo(centerIdx));

            halfStarVertices[i].push_back(par->edges(currentEdge).theOtherVertex(centerIdx));
            halfStarTris[i].push_back(currentTri);

            // look for next edge
            int j;
            for (j=0; j<3; j++){

                nextCurrentEdge = par->triangles(currentTri).edges[j];
                if (nextCurrentEdge != currentEdge &&
                    (par->edges(nextCurrentEdge).isConnectedTo(centerIdx)))
                    break;
            }

            assert(j<3);

            currentEdge = nextCurrentEdge;

            // there is a loop in the face graph which doesn't visit featureEdgeB
            if (currentEdge==featureEdgeA)
                return false;

            // look for next triangle
            if (currentEdge != featureEdgeB){
                // here we assume that currentEdge is NOT a feature edge
                // i.e. it is bordered by exactly two triangles
                if (par->edges(currentEdge).triangles[0] == currentTri)
                    currentTri = par->edges(currentEdge).triangles[1];
                else
                    currentTri = par->edges(currentEdge).triangles[0];

                assert(par->triangles(currentTri).patch==patches[i]);
            }

        }

        halfStarVertices[i].push_back(par->edges(featureEdgeB).theOtherVertex(centerIdx));
    }

    for (i=0; i<halfStarTris.size(); i++)
        for (j=0; j<halfStarTris[i].size(); j++)
            assert(par->triangles(halfStarTris[i][j]).patch == patches[i]);

    // have all edges been visited??
    int vertexCounter = 2;
    for (i=0; i<halfStarVertices.size(); i++)
        vertexCounter += (halfStarVertices[i].size()-2);

    if (par->vertices(centerIdx).degree() != vertexCounter){
        printf("two touching sheets found!\n");
        return false;
    }

    return true;
}

void ParamToolBox::makeFullStarOutOfHalfStars(std::vector<int>& halfStarVerticesA, std::vector<int>& halfStarTrisA,
                                              std::vector<int>& halfStarVerticesB, std::vector<int>& halfStarTrisB,
                                              std::vector<int>& fullStarVertices, std::vector<int>& fullStarTris)
{
        fullStarVertices = halfStarVerticesA;
        fullStarVertices.pop_back();
        std::reverse(halfStarVerticesB.begin(),halfStarVerticesB.end());
        fullStarVertices.insert(fullStarVertices.end(), halfStarVerticesB.begin(), halfStarVerticesB.end());
        fullStarVertices.pop_back();

        fullStarTris = halfStarTrisA;
        std::reverse(halfStarTrisB.begin(),halfStarTrisB.end());
        fullStarTris.insert(fullStarTris.end(), halfStarTrisB.begin(), halfStarTrisB.end());
}

void ParamToolBox::mergeTwoTrianglesIntoQuadrangle(int tri1Idx, int tri2Idx,
                                                   DomainPolygon& quadri, bool& flipped, StaticVector<float,2> quadCoords[4],
                                                   std::vector<unsigned int>& nodeStack,
                                                   const PSurface<2,float>* par)
{
    int i, j;
    DomainTriangle<float>& tri1 = const_cast<DomainTriangle<float>&>(par->triangles(tri1Idx));
    DomainTriangle<float>& tri2 = const_cast<DomainTriangle<float>&>(par->triangles(tri2Idx));

    // find the two common points
    int commonVertices[3] = {-1, -1, -1};
    int count=0;

    for (i=0; i<3; i++)
        for (j=0; j<3; j++)
            if (tri1.vertices[i]==tri2.vertices[j])
                commonVertices[count++] = tri1.vertices[i];

    assert(commonVertices[0]>=0 && commonVertices[1]>=0);
    assert(commonVertices[0]!=commonVertices[1]);

    if (commonVertices[2] != -1){
        printf("Coplanar Triangle found!\n");
        assert(false);
    }

    const Vertex<float> &thirdPoint1 = par->vertices(tri1.getThirdVertex(commonVertices[0], commonVertices[1]));
    const Vertex<float> &thirdPoint2 = par->vertices(tri2.getThirdVertex(commonVertices[0], commonVertices[1]));

    while (tri1.vertices[1] != tri1.getThirdVertex(commonVertices[0], commonVertices[1])){
        tri1.rotate();
    }



    while (tri2.vertices[0] != tri1.vertices[2]){
        tri2.rotate();
    }

    flipped = false;
    if (tri2.vertices[2] != tri1.vertices[0]){
        flipped = true;
        tri2.flip();
    }

    assert(tri1.vertices[0]==tri2.vertices[2] && tri1.vertices[2]==tri2.vertices[0]);

    commonVertices[0] = tri1.vertices[0];
    commonVertices[1] = tri1.vertices[2];

    // 'folds up' the two triangles along their common edge
    // that's obviously the most distortion-free mapping into R2

    const Vertex<float> &from = par->vertices(commonVertices[0]);
    const Vertex<float> &to   = par->vertices(commonVertices[1]);

    //  move the first triangle into the real plane
    float l  = (from-to).length();
    float l1 = (thirdPoint1-from).length();
    float l2 = (thirdPoint1-to).length();

    float alpha = (l2*l2 - l1*l1 - l*l);
    float x1    = alpha / (-2*l);
    float y1    = -sqrt(l1*l1 - alpha * alpha / (4*l*l));

    // same for the second triangle
    l1 = (thirdPoint2-from).length();
    l2 = (thirdPoint2-to).length();

    alpha    = (l2*l2 - l1*l1 - l*l);
    float x2 = alpha / (-2*l);
    float y2 = sqrt(l1*l1 - alpha * alpha / (4*l*l));

    quadCoords[0] = StaticVector<float,2>(0, 0);
    quadCoords[1] = StaticVector<float,2>(x1, y1);
    quadCoords[2] = StaticVector<float,2>(l, 0);
    quadCoords[3] = StaticVector<float,2>(x2, y2);


    // if the quadilateral is not convex (which is exceptional), we move its
    // vertices a little
    if (orientation(quadCoords[1], quadCoords[2], quadCoords[3])!=1) {
        float newL = x1 - y1*(x2-x1)/(y2-y1);
        newL *= 1.1;

        quadCoords[2][0] = newL;

    } else if (orientation(quadCoords[3], quadCoords[0], quadCoords[1])!=1) {
        float newL = x1 - y1*(x2-x1)/(y2-y1);
        newL -= 0.001;
        newL *= 1.1;

        quadCoords[0][0] = newL;
    }

    assert((orientation(quadCoords[1], quadCoords[2], quadCoords[3])!=-1 &&
           orientation(quadCoords[3], quadCoords[0], quadCoords[1])!=-1));

    // merge the two triangles into one quadrilateral
    StaticVector<float,2> secondTriCoords[3];
    secondTriCoords[0][0] = quadCoords[2][0];
    secondTriCoords[0][1] = quadCoords[2][1];
    secondTriCoords[1][0] = quadCoords[3][0];
    secondTriCoords[1][1] = quadCoords[3][1];
    secondTriCoords[2][0] = quadCoords[0][0];
    secondTriCoords[2][1] = quadCoords[0][1];

    int dummy;

    // stupid: copy from StaticVector<float,2> to StaticVector
    StaticVector<float,2> psurfaceQuadCoords[4];
    for (int i=0; i<4; i++) {
        psurfaceQuadCoords[i][0] = quadCoords[i][0];
        psurfaceQuadCoords[i][1] = quadCoords[i][1];
    }

    quadri.init(tri1, psurfaceQuadCoords);

    tri2.checkConsistency("illegal index");
    quadri.mergeTriangle(tri2Idx, secondTriCoords, dummy, nodeStack);

}

////////////////////////////////////////////////////////////////////////
// this routine turns a star shaped polygon into a convex one while
// deforming it as little as possible
void ParamToolBox::convexify(std::vector<StaticVector<float,2> >& coords)
{
    for (int j=0; j<coords.size(); j++)
        coords[j].normalize();

#if 0 // The actual code, needs to be debugged
    int j;
    const signed char counterclockwise = 1;
    const int N = flatCoords.size();

    for (j=0; j<N; j++)
        if (StaticVector<float,2>::orientation(flatCoords[j], flatCoords[(j+1)%N], flatCoords[(j+2)%N])!=counterclockwise){
            const StaticVector<float,2>& a = flatCoords[j];
            const StaticVector<float,2>& b = flatCoords[(j+2)%N];
            StaticVector<float,2>&       c = flatCoords[(j+1)%N];
            float lambda = (a.y*(b.x-a.x) - a.x*(b.y-a.y)) / (c.y*(b.x-a.x) - c.x*(b.y-a.y));
            lambda *= 1.1f;
            c *= lambda;

            j=0;
        }

    for (j=0; j<N; j++)
        assert(StaticVector<float,2>::orientation(flatCoords[j], flatCoords[(j+1)%N], flatCoords[(j+2)%N])==counterclockwise);
#endif
}


void ParamToolBox::convexifyHalfStar(std::vector<StaticVector<float,2> >& coords)
{

//     printf("convexify before:\n");
//     for (int i=0; i<coords.size(); i++)
//         printf("   (%f %f)\n", coords[i].x, coords[i].y);

    StaticVector<float,2> c = 0.5*(coords[0]+coords.back());

    for (int i=1; i<coords.size()-1; i++) {

        StaticVector<float,2>& a = coords[i];
        float radius = 0.5*(coords[0]-coords.back()).length();
        float radiusSquared = radius*radius;

        float squareRoot = 4*a.dot(c)*a.dot(c) - 4*a.length2()*(c.length2()-radiusSquared);

        float lambda1 = (2*a.dot(c) + sqrt(squareRoot)) / (2*a.length2());
//         float lambda2 = (2*a.dot(c) - sqrt(squareRoot)) / 2*a.length2();

//         printf("c (%f %f)   radius %f\n", c.x, c.y, radius);

//         printf("squareRoot %f,  lambda1 %f,  lambda2 %f\n", squareRoot, lambda1, lambda2);

//         printf("2normar %f    sqrt %f\n", 2*a.length()*radius, sqrt(squareRoot));

//         printf("r/norma %f\n", radius / a.length());

        a = lambda1*a;

    }

}


bool ParamToolBox::singleTetrahedronTest(const PSurface<2,float>* par, const std::vector<int>& fullStarVertices)
{
    if (fullStarVertices.size() != 3)
        return false;

    return (par->findTriangle(fullStarVertices[0], fullStarVertices[1], fullStarVertices[2]) != -1);
}

bool ParamToolBox::mergeStarIntoPolygon(int centerIdx, DomainPolygon& fullStar,
                                        std::vector<int>& fullStarTris, int& newCenterNode,
                                        std::vector<unsigned int>& nodeStack,
                                        PSurface<2,float>* par)
{
    int i, j;

    Vertex<float>& centerVertex = par->vertices(centerIdx);

    std::vector<std::vector<int> >   halfStarTris(0);
    std::vector<std::vector<int> >     halfStarVertices(0);

    std::vector<StaticVector<float,2> >    flatCoords;

    std::vector<int> fullStarVertices(0);

    StaticVector<float,2> triCoords[3];

    int featureEdgeA, featureEdgeB;

    bool flipped;

    int featureStatus = ParamToolBox::computeFeatureStatus(par, centerIdx, featureEdgeA, featureEdgeB);

    if (featureStatus != ParamToolBox::REGULAR_POINT)
        return false;

    // compute the star around this vertex:
    // any two edges will do here:
    featureEdgeA = centerVertex.edges[0];
    featureEdgeB = centerVertex.edges[1];

    // finds the two halfstars that make up the full star
    std::vector<int> patches;
    if (!findAllHalfStars(centerIdx, featureEdgeA, featureEdgeB, halfStarVertices, halfStarTris, patches, par))
        return false;

    makeFullStarOutOfHalfStars(halfStarVertices[0], halfStarTris[0],
                               halfStarVertices[1], halfStarTris[1],
                               fullStarVertices, fullStarTris);

    assert(fullStarTris.size()>=3);

    for (i=0; i<fullStarTris.size(); i++)
        par->triangles(fullStarTris[i]).checkConsistency("PreRelaxation\n");

    //flatten the star
    flattenStar(centerIdx, fullStarVertices, flatCoords, par);

    //convexify, if needed
    convexify(flatCoords);

    // merge the flattened star into one domainPolygon

    // turn the first triangle such that its points coincide with the first
    // three points in #fullStarVertices#

    DomainTriangle<float>& firstTri = par->triangles(fullStarTris[0]);

    while (firstTri.vertices[0]!=centerIdx){
        // turn triangle
        firstTri.rotate();
    }

    // flip star, if necessary
    flipped = false;
    if (firstTri.vertices[1] != fullStarVertices[0]){
        flipped = true;

        std::reverse(fullStarVertices.begin(),fullStarVertices.end());
        std::rotate(fullStarVertices.begin(),fullStarVertices.end()-2,fullStarVertices.end());

        std::reverse(fullStarTris.begin(), fullStarVertices.end());
        std::rotate(fullStarTris.begin(),fullStarTris.end()-1,fullStarTris.end());

    }

    ///////////////////
    // merge the star
    for (j=0; j<3; j++)
        if (par->triangles(fullStarTris[0]).vertices[j]==centerIdx){
            triCoords[j][0] = 0;
            triCoords[j][1] = 0;
            break;
        }

    triCoords[(j+1)%3][0] = (par->triangles(fullStarTris[0]).vertices[(j+1)%3] == fullStarVertices[0]) ? flatCoords[0][0] : flatCoords[1][0];
    triCoords[(j+1)%3][1] = (par->triangles(fullStarTris[0]).vertices[(j+1)%3] == fullStarVertices[0]) ? flatCoords[0][1] : flatCoords[1][1];

    triCoords[(j+2)%3][0] = (par->triangles(fullStarTris[0]).vertices[(j+2)%3] == fullStarVertices[0]) ? flatCoords[0][0] : flatCoords[1][0];
    triCoords[(j+2)%3][1] = (par->triangles(fullStarTris[0]).vertices[(j+2)%3] == fullStarVertices[0]) ? flatCoords[0][1] : flatCoords[1][1];

    fullStar.init(par->triangles(fullStarTris[0]), triCoords);

    for (i=1; i<fullStarTris.size(); i++){

        const std::array<int, 3>& currentTriPoints = par->triangles(fullStarTris[i]).vertices;

        for (j=0; j<3; j++)
            if (currentTriPoints[j] == centerIdx){
                triCoords[j][0] = 0;
                triCoords[j][1] = 0;
                break;
            }

        for (int k=0; k<2; k++) {
            triCoords[(j+1)%3][k] = (currentTriPoints[(j+1)%3] == fullStarVertices[i]) ?
                flatCoords[i][k] : flatCoords[(i+1)%flatCoords.size()][k];

            triCoords[(j+2)%3][k] = (currentTriPoints[(j+2)%3] == fullStarVertices[i]) ?
                flatCoords[i][k] : flatCoords[(i+1)%flatCoords.size()][k];
        }

        fullStar.checkConsistency("preMerge\n");
        par->triangles(fullStarTris[i]).checkConsistency("Tri preMerge\n");
        fullStar.mergeTriangle(fullStarTris[i], triCoords, newCenterNode, nodeStack);
        fullStar.checkConsistency("afterMerge\n");
    }

    return true;
}




//////////////////////////////////////////////////////////////////////////////////
// this routine 'flattens' a star given in 3D space into a 2D parameter domain
// using a polar map
//////////////////////////////////////////////////////////////////////////////////

void ParamToolBox::flattenStar(int center, const std::vector<int> &threeDStarVertices,
                               std::vector<StaticVector<float,2> >& twoDVertexPos, const PSurface<2,float>* par)
{
    int k,K = threeDStarVertices.size();

    twoDVertexPos.resize(K);

    std::vector<float>   theta(K+1);

    // compute the (accumulated) angles at the center point
    theta[0] = 0;

    for (k=1; k<K+1; k++){
        const Vertex<float>& pLeft  = par->vertices(threeDStarVertices[k-1]);
        const Vertex<float>& pRight = par->vertices(threeDStarVertices[k%K]);

        theta[k] = theta[k-1] + (pLeft - par->vertices(center)).angle(pRight - par->vertices(center));

    }

    const float a = 2*M_PI/theta.back();

    // compute parameter domain coordinates
    for (k=0; k<K; k++){
        const float r = (par->vertices(threeDStarVertices[k]) - par->vertices(center)).length();
        twoDVertexPos[k] = pow(r, a)*StaticVector<float,2>(cos(theta[k]*a), sin(theta[k]*a));
    }
}

//////////////////////////////////////////////////
// the same thing for a half star
//////////////////////////////////////////////////

void ParamToolBox::flattenHalfStar(int center,
                                   const std::vector<int>& threeDStarVertices,
                                   std::vector<StaticVector<float,2> >& twoDVertexPos,
                                   const PSurface<2,float>* par)
{
    int k,K = threeDStarVertices.size();

    twoDVertexPos.resize(K);

    std::vector<float>   theta(K);

    // compute the (accumulated) angles at the center point
    theta[0] = 0;

    for (k=1; k<K; k++){
        const Vertex<float>& pLeft  = par->vertices(threeDStarVertices[k-1]);
        const Vertex<float>& pRight = par->vertices(threeDStarVertices[k]);

        theta[k] = theta[k-1] + (pLeft - par->vertices(center)).angle(pRight - par->vertices(center));
    }

    float a = M_PI/theta.back();

    // compute parameter domain coordinates
    for (k=0; k<K; k++){

        const float r = (par->vertices(threeDStarVertices[k]) - par->vertices(center)).length();
        const float rPowA = pow(r, a);

        twoDVertexPos[k] = rPowA*StaticVector<float,2>(cos(theta[k]*a), sin(theta[k]*a));
    }
}



///////////////////////////////////////////////////////////////////
// this routine cuts a polygon into pizza slices given an
// appropriate interior graph node
void ParamToolBox::pizzaCutter(DomainPolygon& fullStar, NodeIdx newCenterNode,
                               int& newVertex, std::vector<int>& newTriangles)
{
    int i, j, k;

    if (!fullStar.nodes[newCenterNode].isINTERIOR_NODE())
        printf("Warning:  New centernode is not INTERIOR_NODE!\n");

    PSurface<2,float>* par = fullStar.par;


    ////////////////////////////////////////////////////
    // create new vertex and triangles, if necessary

    if (newVertex==-1)
        newVertex = par->newVertex(StaticVector<float,3>(0));

    if (newTriangles.size() != fullStar.boundaryPoints.size()){
        newTriangles.resize(fullStar.boundaryPoints.size());
        std::fill(newTriangles.begin(),newTriangles.end(),-1U);
    }

    for (i=0; i<newTriangles.size(); i++)
        if (newTriangles[i]<0){
            newTriangles[i] = par->createSpaceForTriangle(newVertex,
                                                          fullStar.boundaryPoints[i],
                                                          fullStar.boundaryPoints[(i+1) %
                                                                                  fullStar.boundaryPoints.size()]);

        }

    // /////////////////////////////////////////////////////////////////
    // cut the slices
    //par->vertices(newVertex) = fullStar.nodes[newCenterNode].getImagePos(par->iPos);
    par->vertices(newVertex) = par->iPos[fullStar.nodes[newCenterNode].getNodeNumber()];

    for (i=0; i<newTriangles.size(); i++){
        fullStar.slice(newCenterNode, newVertex, i*3);
        fullStar.checkConsistency("Slicing");
    }

    ///////////////////////////////////////////////////////////////////
    // the polygon has been cut.  Move each slice to its original triangle.

    int cN;
    std::vector<int> nodeLocs(fullStar.nodes.size());

    for (i=0; i<newTriangles.size(); i++) {

        for (cN=0; cN<fullStar.nodes.size(); cN++)
            nodeLocs[cN] = DomainPolygon::IN_POLYGON;

        DomainTriangle<float>& cT = par->triangles(newTriangles[i]);

//         int offset=0;
//         if (cT.vertices[0]==centerPoint)
//             offset = 0;
//         else if (cT.vertices[1]==centerPoint)
//             offset = 1;
//         else
//             offset = 2;

        int offset = cT.getCorner(newVertex);

        // copy edgePoint arrays
        for (j=0; j<3; j++){
            cT.edgePoints[(j+offset)%3] = fullStar.edgePoints[(3*i+1+j)%fullStar.edgePoints.size()];
            fullStar.edgePoints[(3*i+1+j)%fullStar.edgePoints.size()].clear();
        }

        // copy nodes using a graph-search algorithm
        for (j=0; j<3; j++){
            for (k=0; k<cT.edgePoints[j].size(); k++)
                if (nodeLocs[cT.edgePoints[j][k]] != DomainPolygon::IN_TRIANGLE)
                    moveSubGraph(cT.edgePoints[j][k], fullStar, nodeLocs, newCenterNode);
        }

        // make a copy of the centerNode
        fullStar.nodes.push_back(Node<float>());
        //int localCenterNode = fullStar.nodes.appendSpace(1);
        int localCenterNode = fullStar.nodes.size()-1;

        fullStar.nodes[localCenterNode].setValue(fullStar.nodes[newCenterNode].domainPos(),
                                                 fullStar.nodes[newCenterNode].getNodeNumber(),
                                                 Node<float>::CORNER_NODE);

        nodeLocs.push_back(DomainPolygon::IN_TRIANGLE);

        for (j=fullStar.nodes[newCenterNode].degree()-1; j>=0; j--)
            if (nodeLocs[fullStar.nodes[newCenterNode].neighbors(j)] == DomainPolygon::IN_TRIANGLE) {
                fullStar.nodes[localCenterNode].appendNeighbor(fullStar.nodes[newCenterNode].neighbors(j));
                fullStar.nodes[fullStar.nodes[newCenterNode].neighbors(j)].replaceReferenceTo(newCenterNode, localCenterNode);
                fullStar.nodes[newCenterNode].removeNeighbor(j);
            }

        cT.edgePoints[0+offset][0] = localCenterNode;
        cT.edgePoints[(2+offset)%3].back() = localCenterNode;

        //////////////////////////////////////
        // sort out the nodes that belong onto the triangle
        int numTriNodes = 0;
        int triNode;

        for (triNode=0; triNode<fullStar.nodes.size(); triNode++)
            if (nodeLocs[triNode] == DomainPolygon::IN_TRIANGLE)
                numTriNodes++;

        int triCount = 0;
        cT.nodes.resize(numTriNodes);
        std::vector<int> offArr(fullStar.nodes.size());

        for (triNode=0; triNode<fullStar.nodes.size(); triNode++)
            if (nodeLocs[triNode] == DomainPolygon::IN_TRIANGLE) {
                cT.nodes[triCount] = fullStar.nodes[triNode];
                fullStar.invalidate(triNode);
                offArr[triNode] = triCount;
                triCount++;
            }

        for (j=0; j<numTriNodes; j++)
            for (k=0; k<cT.nodes[j].degree(); k++)
                cT.nodes[j].neighbors(k) = offArr[cT.nodes[j].neighbors(k)];

        for (j=0; j<3; j++)
            for (k=0; k<cT.edgePoints[j].size(); k++)
                cT.edgePoints[j][k] = offArr[cT.edgePoints[j][k]];



        /////////////////////////////////////
        cT.checkConsistency("After Vertex Relaxation\n");

        cT.installBarycentricCoordinates();


        fullStar.checkConsistency("before garbage collection\n");

        // reuse offArr
        fullStar.garbageCollection(offArr);
        newCenterNode -= offArr[newCenterNode];
        for (j=0; j<offArr.size(); j++)
            nodeLocs[j-offArr[j]] = nodeLocs[j];

        fullStar.checkConsistency("after garbage collection\n");

    }

}


void ParamToolBox::moveSubGraph(int startingNode, DomainPolygon& from, std::vector<int>& nodeLocs, int centerNode)
{
    if (startingNode==centerNode)
        return;

    nodeLocs[startingNode] = DomainPolygon::IN_TRIANGLE;

    for (int i=0; i<from.nodes[startingNode].degree(); i++)
        if (nodeLocs[from.nodes[startingNode].neighbors(i)] != DomainPolygon::IN_TRIANGLE)
            moveSubGraph(from.nodes[startingNode].neighbors(i), from, nodeLocs, centerNode);
}


bool ParamToolBox::removeRegularPoint(PSurface<2,float>* par, int centerPoint, const QualityRequest &quality,
                                      MultiDimOctree<Edge, EdgeIntersectionFunctor, float, 3>* edgeOctree)
{
    std::vector<unsigned int> nodeStack;
    int i, j;

    std::vector<std::vector<int> > halfStarTris(0);
    std::vector<std::vector<int> > halfStarVertices(0);

    std::vector<int>               patches;

    std::vector<int> fullStarVertices(0);
    std::vector<int> fullStarTris(0);

    // the point is center of a star that's homeomorphic to
    // a disc.  It can be removed without any special handling

    // any two edges will do here:
    int featureEdgeA = par->vertices(centerPoint).edges[0];
    int featureEdgeB = par->vertices(centerPoint).edges[1];

    // finds the two halfstars that make up the full star
    if (!ParamToolBox::findAllHalfStars(centerPoint, featureEdgeA, featureEdgeB,
                                        halfStarVertices, halfStarTris, patches, par)){
        return false;
    }

    ParamToolBox::makeFullStarOutOfHalfStars(halfStarVertices[0], halfStarTris[0],
                                             halfStarVertices[1], halfStarTris[1],
                                             fullStarVertices, fullStarTris);

    if (singleTetrahedronTest(par, fullStarVertices)){
        return false;
    }

    // test for coplanar triangles
    if (fullStarTris.size()<3) {
        printf("coplanar triangle found!\n");
        return false;
    }

    // turn the first triangle such that its points coincide with the first
    // three points in #fullStarVertices#

    DomainTriangle<float>& firstTri = par->triangles(fullStarTris[0]);

    while (firstTri.vertices[0] != centerPoint)
        firstTri.rotate();


    // flip star, if necessary
    bool flipped = false;
    if (firstTri.vertices[1] != fullStarVertices[0]){
        flipped = true;

        std::reverse(fullStarVertices.begin(),fullStarVertices.end());
        std::rotate(fullStarVertices.begin(), fullStarVertices.end()-2,fullStarVertices.end());

        std::reverse(fullStarTris.begin(),fullStarTris.end());
        std::rotate(fullStarTris.begin(),fullStarTris.end()-1,fullStarTris.end());

    }

    // flatten star and do a retriangulation of the hole
    std::vector<StaticVector<float,2> > diskCoords;
    CircularPatch<float> fillIn(fullStarTris.size()-2, par);


    Triangulator::triangulateStar(fullStarVertices, centerPoint, fillIn, diskCoords, par);

    // test for possible topology changes
    if (fillIn.inducesTopologyChange()){
        fillIn.killAll();
        return false;
    }

    // test for small dihedral angles
    if (quality.smallDihedralAngles &&
        fillIn.hasSmallDihedralAngles(quality.dihedralAngleThreshold, par, &par->vertices(centerPoint))) {
        printf("Small Dihedral Angle found\n");
        fillIn.killAll();
        return false;
    }

    // convexify the polygon if needed
    //ParamToolBox::convexify(diskCoords);
    for (j=0; j<diskCoords.size(); j++)
        diskCoords[j].normalize();

    // merge the star

    StaticVector<float,2> triCoords[3];

    for (j=0; j<3; j++)
        if (par->triangles(fullStarTris[0]).vertices[j]==centerPoint){
            triCoords[j][0] = 0;
            triCoords[j][1] = 0;
            break;
        }

    for (int i=0; i<2; i++) {
        triCoords[(j+1)%3][i] = (par->triangles(fullStarTris[0]).vertices[(j+1)%3] == fullStarVertices[0]) ? diskCoords[0][i] : diskCoords[1][i];
        triCoords[(j+2)%3][i] = (par->triangles(fullStarTris[0]).vertices[(j+2)%3] == fullStarVertices[0]) ? diskCoords[0][i] : diskCoords[1][i];
    }

    DomainPolygon fullStar(par);
    fullStar.init(par->triangles(fullStarTris[0]), triCoords);

    for (i=1; i<fullStarTris.size(); i++){

        const std::array<int, 3>& currentTriPoints = par->triangles(fullStarTris[i]).vertices;

        for (j=0; j<3; j++)
            if (currentTriPoints[j] == centerPoint){
                triCoords[j][0] = triCoords[j][1] = 0;
                break;
            }

        for (int k=0; k<2; k++) {
            triCoords[(j+1)%3][k] = (currentTriPoints[(j+1)%3] == fullStarVertices[i]) ?
                diskCoords[i][k] : diskCoords[(i+1)%diskCoords.size()][k];
            triCoords[(j+2)%3][k] = (currentTriPoints[(j+2)%3] == fullStarVertices[i]) ?
                diskCoords[i][k] : diskCoords[(i+1)%diskCoords.size()][k];
        }

        fullStar.checkConsistency("preMerge");
        par->triangles(fullStarTris[i]).checkConsistency("Tri preMErge\n");

        int dummy;
        fullStar.mergeTriangle(fullStarTris[i], triCoords, dummy, nodeStack);
        fullStar.checkConsistency("afterMerge");

    }

    // correct possible foldovers in the parametrization
    fullStar.unflipTriangles();

    // triangulate
    std::vector<unsigned int> tempNodeStack = nodeStack;

    if (!fullStar.triangulate(fillIn, tempNodeStack)){
        // the retriangulation has failed.  This can occur due to numerical problems.
        // Up to this point we have not actually changed anything yet.
        // therefore we can simply abort the removal process!  (ain't that nice?)
        printf("rejecting point because of numerical problems while triangulating!\n");
        fillIn.killAll();
        return false;
    }

    nodeStack = tempNodeStack;

    // /////////////////////////////////////////////////////////
    // incorporate new triangle group

    // first remove old edges from the octree
    if (quality.intersections){
        for (i=0; i<par->vertices(centerPoint).degree(); i++){
            //printf("removing edge %d\n", par->vertices(centerPoint).edges[i]);
            edgeOctree->remove(&par->edges(par->vertices(centerPoint).edges[i]));
        }
    }

    // remove the old triangles from the base grid
    for (i=0; i<fullStarTris.size(); i++){
        //printf("removing tri %d\n", fullStarTris[i]);
        par->removeTriangle(fullStarTris[i]);
    }

    // delete the vertex
    par->removeVertex(centerPoint);

    // add the new triangles to the base grid
    for (i=0; i<fillIn.size(); i++){
        par->triangles(fillIn[i]).checkConsistency("newDomainTriangle");
        par->triangles(fillIn[i]).patch = patches[0];

        par->integrateTriangle(fillIn[i]);
    }

    // append new edges to octree
    if (quality.intersections){

        EdgeIntersectionFunctor edgeIntersectionFunctor(&par->vertices(0));

        for (i=0; i<fillIn.innerEdges.size(); i++) {
            //printf("inserting edge %d\n", par->findEdge(fillIn.innerEdges[i][0], fillIn.innerEdges[i][1]));
            Edge& cE = par->edges(par->findEdge(fillIn.innerEdges[i][0], fillIn.innerEdges[i][1]));

            edgeOctree->insert(&cE);
        }
    }

    return true;
}


bool ParamToolBox::removeFeatureLinePoint(PSurface<2,float>* par,
                                          int centerPoint,
                                          const QualityRequest &quality,
                                          int numHalfStars,
                                          int featureEdgeA,
                                          int featureEdgeB,
                                          MultiDimOctree<Edge, EdgeIntersectionFunctor, float, 3>* edgeOctree)
{
    int i, j;

    std::vector<std::vector<int> > halfStarTris(0);
    std::vector<std::vector<int> > halfStarVertices(0);

    std::vector<int>               patches;

    std::vector<unsigned int> nodeStack;

    // the point sits on a feature line of order #featureStatus#
    // i.e. the feature lines is the touching point of #featureStatus#
    // halfstars that have to be parametrized separately.
    if (!ParamToolBox::findAllHalfStars(centerPoint, featureEdgeA, featureEdgeB,
                                        halfStarVertices, halfStarTris, patches, par)){
        printf("Problems finding all Halfstars!\n");
        return false;
    }

    int newFeatureEdgeFrom = par->edges(featureEdgeA).theOtherVertex(centerPoint);
    int newFeatureEdgeTo   = par->edges(featureEdgeB).theOtherVertex(centerPoint);

    // test for coplanar triangles
    int numUnitaryHalfStars = 0;
    for (i=0; i<numHalfStars; i++)
        numUnitaryHalfStars += (halfStarTris[i].size()==1) ? 1 : 0;

    if (numUnitaryHalfStars>1) {
        printf("coplanar triangle found!\n");
        return false;
    }

    // test for topology change
    if (par->findEdge(newFeatureEdgeFrom, newFeatureEdgeTo)!=-1 && numUnitaryHalfStars==0){
        printf("Topology change\n");
        return false;
    }

    // test whether the removal of the vertex would lead to a triangle consisting of
    // three #nodes# that are all on the same DomainTriangle edge
    for (i=0; i<numHalfStars; i++){

        bool singleConnection = true;

        for (j=0; j<halfStarTris[i].size(); j++) {

            //              if (debugCounter>479)
            //                  printf("i %d   j %d\n", i, j);

            const DomainTriangle<float>& cT = par->triangles(halfStarTris[i][j]);

            int corner = cT.getCorner(centerPoint);
            assert(corner>=0);

            int forwardNeighbor  = cT.edgePoints[corner][1];
            int backwardNeighbor = cT.edgePoints[(corner+2)%3][cT.edgePoints[(corner+2)%3].size()-2];

            singleConnection &= cT.nodes[backwardNeighbor].isConnectedTo(forwardNeighbor);

        }

        //printf("SingleConnection = %d\n", singleConnection);
        if (singleConnection){
            //printf("Single connection\n");
            return false;
        }
    }

    // ///////////////////////////////////////////////////////////////////////////
    // check whether there is a halfstar that consists of only one triangle and
    // can't even be removed with our special handling below
    for (i=0; i<numHalfStars; i++){
        if (halfStarTris[i].size()==1){
            printf("special case ... ");

            const Edge& oppEdge = par->edges(par->triangles(halfStarTris[i][0]).getOppositeEdge(centerPoint));
            if (oppEdge.numTriangles()!=2 ||
                par->triangles(oppEdge.triangles[0]).patch !=
                par->triangles(oppEdge.triangles[1]).patch){
                printf("rejected B\n");
                return false;
            }
            printf("accepted\n");
        }
    }

    // ///////////////////////////////////////////////////

    // \bug This comes to early
    if (quality.intersections){

        // first remove old edges from the octree
        for (i=0; i<par->vertices(centerPoint).degree(); i++)
            edgeOctree->remove(&par->edges(par->vertices(centerPoint).edges[i]));
    }


    // first do the retriangulations
    std::vector<CircularPatch<float> > fillIns(numHalfStars);
    std::vector<std::vector<StaticVector<float,2> > > flatCoords(numHalfStars);

    for (i=0; i<numHalfStars; i++){
        fillIns[i].par = par;
        if (halfStarTris[i].size()!=1) {

            // merge the flattened half star into one DomainPolygon
            DomainTriangle<float>& firstTri = par->triangles(halfStarTris[i][0]);

            while (firstTri.vertices[0] != centerPoint)
                firstTri.rotate();

            // do a retriangulation of the hole
            fillIns[i].resize(halfStarTris[i].size()-1);
            Triangulator::triangulateHalfStar(halfStarVertices[i],
                                              centerPoint, fillIns[i], flatCoords[i], par);

            // test for possible topology changes
            if (fillIns[i].inducesTopologyChange()){
                for (j=0; j<fillIns.size(); j++)
                    fillIns[j].killAll();

                printf("Topology change B\n");
                return false;
            }

            // test for small dihedral angles
            if (quality.smallDihedralAngles &&
                fillIns[i].hasSmallDihedralAngles(quality.dihedralAngleThreshold, par,
                                                  &par->vertices(centerPoint))) {
                for (j=0; j<fillIns.size(); j++)
                    fillIns[j].killAll();

                printf("Small Dihedral Angles\n");
                return false;
            }

            //
            ParamToolBox::convexifyHalfStar(flatCoords[i]);

        }
    }

    if (quality.smallDihedralAngles &&
        (par->smallestDihedralAngle(featureEdgeA) < quality.dihedralAngleThreshold ||
         par->smallestDihedralAngle(featureEdgeB) < quality.dihedralAngleThreshold)) {
        for (j=0; j<fillIns.size(); j++)
            fillIns[j].killAll();

        printf("Small Dihedral Angles B\n");
        return false;
    }

    // first do all the recutting of the halfstars.  Proceed only when
    // all have been successfully cut.

    for (i=0; i<numHalfStars; i++){
        // special handling for stars that consist of only one triangle
        if (halfStarTris[i].size()==1) {
            // \bug This code is buggy!  If it performs a modification and
            // then triangulate fails on another iteration of the outer
            // loop --> the parametrization data structure is corrupted.
            // One really needs to do all calls to triangulate first, see
            // if they work, and *only* if they all do, modify the base grid

            printf("%d Halfstar with one triangle\n", i);
            printf("WARNING:  THIS CODE POTENTIALLY LEADS TO BUGS!\n");
            int tri1Idx = halfStarTris[i][0];
            const Edge& outerEdge = par->edges(par->triangles(tri1Idx).getOppositeEdge(centerPoint));
            int tri2Idx = (outerEdge.triangles[0]==halfStarTris[i][0])
                ? outerEdge.triangles[1]
                : outerEdge.triangles[0];

            StaticVector<float,2> quadCoords[4];
            DomainPolygon quadri(par);
            bool flipped;

            ParamToolBox::mergeTwoTrianglesIntoQuadrangle(tri1Idx, tri2Idx,
                                                          quadri, flipped, quadCoords, nodeStack, par);


            // make quadrangle a topological triangle
            quadri.removeVertex(centerPoint);

            // make it also a geometrical triangle
            for (j=0; j<quadri.edgePoints[0].size(); j++)
                quadri.nodes[quadri.edgePoints[0][j]].setDomainPos(quadri.nodes[quadri.edgePoints[0][0]].domainPos()
                                                                   + j*(quadri.nodes[quadri.edgePoints[0].back()].domainPos()
                                                                        - quadri.nodes[quadri.edgePoints[0][0]].domainPos()) /
                                                                   (quadri.edgePoints[0].size()-1));

            quadri.checkConsistency("small halfstar");
            quadri.applyParametrization();

            // revert back to triangle
            CircularPatch<float> cP(1, par);
            cP[0] = tri2Idx;
            //par->checkConsistency("EE.7");
            std::vector<unsigned int> tempNodeStack = nodeStack;
            if (!quadri.triangulate(cP, tempNodeStack)){
                printf("[removeFeatureLinePoint] triangulate failed A!!\n");
                for (j=0; j<fillIns.size(); j++)
                    fillIns[j].killAll();
                return false;
            }

            par->integrateTriangle(cP[0]);

            nodeStack = tempNodeStack;

            if (flipped)
                par->triangles(cP[0]).flip();

            cP[0]=-1U;
            par->removeTriangle(tri1Idx);
            par->triangles(tri2Idx).checkConsistency("tri2");

        }else {

            DomainTriangle<float>& firstTri = par->triangles(halfStarTris[i][0]);

            while (firstTri.vertices[0] != centerPoint)
                firstTri.rotate();

            // flip triangle, if necessary
            bool flipped = false;
            if (firstTri.vertices[1] != halfStarVertices[i][0]){
                flipped = true;
                firstTri.flip();
            }

            const std::array<int, 3>& firstTriPoints = par->triangles(halfStarTris[i][0]).vertices;

            StaticVector<float,2> triCoords[3];

            for (j=0; j<3; j++)
                if (firstTriPoints[j] == centerPoint){
                    triCoords[j][0] = triCoords[j][1] = 0;
                    break;
                }

            for (int k=0; k<2; k++) {
                triCoords[(j+1)%3][k] = (firstTriPoints[(j+1)%3] == halfStarVertices[i][0]) ? flatCoords[i][0][k] : flatCoords[i][1][k];
                triCoords[(j+2)%3][k] = (firstTriPoints[(j+2)%3] == halfStarVertices[i][0]) ? flatCoords[i][0][k] : flatCoords[i][1][k];
            }

            DomainPolygon halfStar(par);
            halfStar.init(par->triangles(halfStarTris[i][0]), triCoords);

            for (int k=1; k<halfStarTris[i].size(); k++){

                const std::array<int, 3>& currentTriPoints = par->triangles(halfStarTris[i][k]).vertices;

                for (j=0; j<3; j++)
                    if (currentTriPoints[j] == centerPoint){
                        triCoords[j][0] = triCoords[j][1] = 0;
                        break;
                    }

                for (int l=0; l<2; l++) {
                    triCoords[(j+1)%3][l] = (currentTriPoints[(j+1)%3] == halfStarVertices[i][k]) ?
                        flatCoords[i][k][l] : flatCoords[i][(k+1)%flatCoords[i].size()][l];
                    triCoords[(j+2)%3][l] = (currentTriPoints[(j+2)%3] == halfStarVertices[i][k]) ?
                        flatCoords[i][k][l] : flatCoords[i][(k+1)%flatCoords[i].size()][l];
                }

                int dummy;
                halfStar.mergeTriangle(halfStarTris[i][k], triCoords, dummy, nodeStack);
            }
            // remove centerPoint as a polygon vertex, it remains as a touching node
            //halfStar.checkConsistency("feature before removeVertex");

            halfStar.removeVertex(centerPoint);
            //halfStar.checkConsistency("feature after removeVertex");

            // correct possible foldovers in the parametrization
            halfStar.unflipTriangles();

            // triangulate
            halfStar.checkConsistency("feature before triangulate");

            std::vector<unsigned int> tempNodeStack = nodeStack;

            if (!halfStar.triangulate(fillIns[i], tempNodeStack)) {
                printf("[removeFeatureLinePoint] triangulate failed! B\n");
                for (j=0; j<fillIns.size(); j++)
                    fillIns[j].killAll();
                return false;
            }
            nodeStack = tempNodeStack;

            for (j=0; j<fillIns[i].size(); j++){
                if (flipped)
                    par->triangles(fillIns[i][j]).flip();

                par->triangles(fillIns[i][j]).checkConsistency("feature newDomainTriangle");
                par->triangles(fillIns[i][j]).patch = patches[i];

            }

        }
    }

    // All halfstars have been successfully cut.  Now when safely enter them
    // into the surface
    for (i=0; i<numHalfStars; i++){
        // special handling for stars that consist of only one triangle
        if (halfStarTris[i].size()==1) {

        }else {

            // incorporate new triangle group
            for (j=0; j<halfStarTris[i].size(); j++)
                par->removeTriangle(halfStarTris[i][j]);

            for (j=0; j<fillIns[i].size(); j++){
//                 if (flipped)
//                     par->triangles(fillIns[i][j]).flip();

                par->triangles(fillIns[i][j]).checkConsistency("feature newDomainTriangle");
                //par->triangles(fillIns[i][j]).patch = patches[i];

                par->integrateTriangle(fillIns[i][j]);
            }

            // update edge octree, if used
            if (quality.intersections){

                // add new edges
                for (j=0; j<fillIns[i].innerEdges.size(); j++) {

                    Edge& cE = par->edges(par->findEdge(fillIns[i].innerEdges[j][0], fillIns[i].innerEdges[j][1]));
                    edgeOctree->insert(&cE);
                }

            }
        }
    }

    // enter new feature edge into octree
    if (quality.intersections){

        int newFeatureEdge = par->findEdge(newFeatureEdgeFrom, newFeatureEdgeTo);
        assert(newFeatureEdge>=0);

        edgeOctree->insert(&par->edges(newFeatureEdge));

    }

    par->removeVertex(centerPoint);

    return true;
}


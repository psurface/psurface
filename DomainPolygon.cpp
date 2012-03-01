#ifdef _MSC_VER
    // Required to make cmath define M_PI etc.
    #define _USE_MATH_DEFINES
#endif
#include <cmath>

#include "DomainPolygon.h"
#include "CircularPatch.h"

// Check for VC9 / VS2008 without SP1, which lacks the C99 math conformance stuff.
#if defined(_MSC_VER) && _MSC_VER==1500
    #include <float.h>

    namespace std {
        inline double isnan(double x) {
            return _isnan(x);
        }
    }
#endif

using namespace psurface;

static int counter = 0;

void DomainPolygon::init(const DomainTriangle<float>& tri, const StaticVector<float,2> coords[3]){

    nodes = tri.nodes;     
    
    boundaryPoints.resize(3);
    boundaryPoints[0] = tri.vertices[0];
    boundaryPoints[1] = tri.vertices[1];
    boundaryPoints[2] = tri.vertices[2];
    
    edgePoints.resize(3);
    edgePoints[0] = tri.edgePoints[0];
    edgePoints[1] = tri.edgePoints[1];
    edgePoints[2] = tri.edgePoints[2];
    
    // turns the coordinates into world coordinates
    PlaneParam<float>::installWorldCoordinates(coords[0], coords[1], coords[2]);
    
    removeExtraEdges();
    par->hasUpToDatePointLocationStructure = false;

}

void DomainPolygon::mergeTriangle(int tri, StaticVector<float,2> coords[3], int& newCenterNode,
                                  std::vector<unsigned int>& nodeStack)
{

    // for debugging

    DomainTriangle<float>& cT = par->triangles(tri);

    cT.removeExtraEdges();
    par->hasUpToDatePointLocationStructure = false;

    int thePolygonEdge[2];
    int theTriangleEdge[2];

    cT.checkConsistency("Tri merge AAAA");
    

    bool reverse;

    thePolygonEdge[0] = -1;
    theTriangleEdge[0] = -1;
    
    // find corresponding edge
    int numMatchingEdges = 0;
    for (size_t i=0; i<boundaryPoints.size(); i++){
        for (int j=0; j<3; j++){

            if ((cT.vertices[j]==boundaryPoints[i] && 
                 cT.vertices[(j+1)%3]==boundaryPoints[(i+1)%boundaryPoints.size()])){

                // a common edge is found and the edgePoint arrays match
                thePolygonEdge[numMatchingEdges]    = i;
                theTriangleEdge[numMatchingEdges] = j;
                
                reverse = false;

                numMatchingEdges++;
            }


            if (cT.vertices[(j+1)%3]==boundaryPoints[i] && 
                cT.vertices[j]==boundaryPoints[(i+1)%boundaryPoints.size()]){

                // a common edge is found and the edgePoint arrays are reversed
                thePolygonEdge[numMatchingEdges]    = i;
                theTriangleEdge[numMatchingEdges] = j;

                reverse = true;

                numMatchingEdges++;
            }
        }
    }

    assert(numMatchingEdges<3);
    
    cT.checkConsistency("Tri merge A");

    int newNodeIdx = nodes.size();

    for (size_t k=0; k<cT.nodes.size(); k++)
        nodes.push_back(cT.nodes[k]);

    installWorldCoordinates(newNodeIdx, coords[0], coords[1], coords[2]);

    std::tr1::array<std::vector<int>, 3> tmpEdgePoints = cT.edgePoints;

    //checkConsistency("Poly merge before augmentNeighborIdx\n");

    augmentNeighborIdx(newNodeIdx, tmpEdgePoints);

    checkConsistency("Poly merge A");

    assert(tmpEdgePoints[theTriangleEdge[0]].size() == edgePoints[thePolygonEdge[0]].size());
    
    // handle all edges that cross the unifying edges
    for (int k=0; k<numMatchingEdges; k++){

        for (int i=1; i<edgePoints[thePolygonEdge[k]].size()-1; i++){
            
            if (nodes[edgePoints[thePolygonEdge[k]][i]].isINTERSECTION_NODE()){
                
                // it's an intersection point
                int intersectPointPoly = edgePoints[thePolygonEdge[k]][i];
                int intersectPointTri;

                edgePoints[thePolygonEdge[k]][i] = -1;

                if (reverse){
                    intersectPointTri  = tmpEdgePoints[theTriangleEdge[k]][tmpEdgePoints[theTriangleEdge[k]].size()-i-1];
                } else {
                    intersectPointTri  = tmpEdgePoints[theTriangleEdge[k]][i];
                }

                assert(nodes[intersectPointPoly].isINTERSECTION_NODE() &&
                       nodes[intersectPointTri].isINTERSECTION_NODE());

                assert(nodes[intersectPointPoly].getNodeNumber()==nodes[intersectPointTri].getNodeNumber());

                int innerPointPoly     = nodes[intersectPointPoly].neighbors(0);
                int innerPointTri      = nodes[intersectPointTri].neighbors(0);
                
                nodes[innerPointTri].replaceReferenceTo(intersectPointTri, innerPointPoly);
                nodes[innerPointPoly].replaceReferenceTo(intersectPointPoly, innerPointTri);

                // remove the two intersection points 
                nodeStack.push_back(nodes[intersectPointPoly].getNodeNumber());

                invalidate(intersectPointPoly);
                invalidate(intersectPointTri);

            }else{
                
                // it is a touching point
                int touchingPointPoly = edgePoints[thePolygonEdge[k]][i];
                int touchingPointTri  = -1;

                if (reverse)
                    touchingPointTri  = tmpEdgePoints[theTriangleEdge[k]][tmpEdgePoints[theTriangleEdge[k]].size()-i-1];
                else
                    touchingPointTri  = tmpEdgePoints[theTriangleEdge[k]][i];


                assert(nodes[touchingPointPoly].isTOUCHING_NODE());
                assert(nodes[touchingPointTri].isTOUCHING_NODE());

                mergeNodes(touchingPointPoly, touchingPointTri);

                nodes[touchingPointPoly].makeInteriorNode();
                    
            }     
        }
    }

    checkConsistency("Poly merge C\n");

    ////////////////////////////////////////
    // the cornerNodes get special handling
    if (numMatchingEdges==1) {

        newCenterNode = -1;

        // one corner
        int cornerNodePoly = edgePoints[thePolygonEdge[0]][0];

        int cornerNodeTri  = (reverse) 
            ? tmpEdgePoints[theTriangleEdge[0]].back() 
            : tmpEdgePoints[theTriangleEdge[0]][0];

        assert(nodes[cornerNodePoly].isCORNER_NODE());
        assert(nodes[cornerNodeTri].isCORNER_NODE());

        mergeNodes(cornerNodePoly, cornerNodeTri);

        // brute force update of the adjacent edgePoint arrays
        updateEdgePoints(tmpEdgePoints, cornerNodeTri, cornerNodePoly);

        // the other corner
        cornerNodePoly = edgePoints[thePolygonEdge[0]].back();

        cornerNodeTri  = (reverse) 
            ? tmpEdgePoints[theTriangleEdge[0]][0] 
            : tmpEdgePoints[theTriangleEdge[0]].back();

        assert(nodes[cornerNodePoly].isCORNER_NODE());
        assert(nodes[cornerNodeTri].isCORNER_NODE());

        mergeNodes(cornerNodePoly, cornerNodeTri);

        // update the references from the adjacent edgePoint arrays
        updateEdgePoints(tmpEdgePoints, cornerNodeTri, cornerNodePoly);

    } else {

        int centerNodePoly, centerNodeTri;

        // find the centerNode on the polygon side
        if (edgePoints[thePolygonEdge[0]][0] == edgePoints[thePolygonEdge[1]].back()){
            centerNodePoly = edgePoints[thePolygonEdge[0]][0];
        }else if (edgePoints[thePolygonEdge[1]][0] == edgePoints[thePolygonEdge[0]].back()){
            centerNodePoly = edgePoints[thePolygonEdge[1]][0];
        }else{
            printf("centerNode NOT found!\n");
            print(true, false, false);
            //triangle->print(true, false, false);
            assert(false);
        }

        // merge centerNode
        if (tmpEdgePoints[theTriangleEdge[0]][0] == tmpEdgePoints[theTriangleEdge[1]].back()){
            centerNodeTri = cT.edgePoints[theTriangleEdge[0]][0];
        }else if (tmpEdgePoints[theTriangleEdge[1]][0] == tmpEdgePoints[theTriangleEdge[0]].back()){
            centerNodeTri = tmpEdgePoints[theTriangleEdge[1]][0];
        }else{
            printf("centerNode Tri NOT found!\n");
            assert(false);
        }

        mergeNodes(centerNodePoly, centerNodeTri);
        nodes[centerNodePoly].makeInteriorNode();
        newCenterNode = centerNodePoly;

        updateEdgePoints(tmpEdgePoints, centerNodeTri, centerNodePoly);

        // the first cornerNode
        int cornerNodePoly, cornerNodeTri;

        if (edgePoints[thePolygonEdge[0]][0]!=centerNodePoly)
            cornerNodePoly = edgePoints[thePolygonEdge[0]][0];
        else 
            cornerNodePoly = edgePoints[thePolygonEdge[0]].back();

        if (tmpEdgePoints[theTriangleEdge[0]][0]!=centerNodePoly)
            cornerNodeTri = tmpEdgePoints[theTriangleEdge[0]][0];
        else 
            cornerNodeTri = tmpEdgePoints[theTriangleEdge[0]].back();

        mergeNodes(cornerNodePoly, cornerNodeTri);

        updateEdgePoints(tmpEdgePoints, cornerNodeTri, cornerNodePoly);

        // the second cornerNode
        if (edgePoints[thePolygonEdge[1]][0]!=centerNodePoly)
            cornerNodePoly = edgePoints[thePolygonEdge[1]][0];
        else 
            cornerNodePoly = edgePoints[thePolygonEdge[1]].back();

        if (tmpEdgePoints[theTriangleEdge[1]][0]!=centerNodePoly)
            cornerNodeTri = tmpEdgePoints[theTriangleEdge[1]][0];
        else 
            cornerNodeTri = tmpEdgePoints[theTriangleEdge[1]].back();

        mergeNodes(cornerNodePoly, cornerNodeTri);

        updateEdgePoints(tmpEdgePoints, cornerNodeTri, cornerNodePoly);

    }

    

    checkConsistency("Poly merge D\n");
        
    if (numMatchingEdges == 1){

        boundaryPoints.insert(boundaryPoints.begin() + thePolygonEdge[0]+1, cT.vertices[(theTriangleEdge[0]+2)%3]);
        
        edgePoints.erase(edgePoints.begin() + thePolygonEdge[0]);

        if (reverse){
            edgePoints.push_back(tmpEdgePoints[(theTriangleEdge[0]+1)%3]);
            edgePoints.push_back(tmpEdgePoints[(theTriangleEdge[0]+2)%3]);
        }else{
            edgePoints.push_back(tmpEdgePoints[(theTriangleEdge[0]+2)%3]);
            edgePoints.push_back(tmpEdgePoints[(theTriangleEdge[0]+1)%3]);

            std::reverse(edgePoints.back().begin(), edgePoints.back().end());
            std::reverse(edgePoints[edgePoints.size()-2].begin(), edgePoints[edgePoints.size()-2].end());
        }
        
    }else{
        
        boundaryPoints.erase(boundaryPoints.begin() + thePolygonEdge[0]);

        if (thePolygonEdge[0] < thePolygonEdge[1]){
            edgePoints.erase(edgePoints.begin() + thePolygonEdge[1]);
            edgePoints.erase(edgePoints.begin() + thePolygonEdge[0]);
        }else{
            edgePoints.erase(edgePoints.begin() + thePolygonEdge[0]);
            edgePoints.erase(edgePoints.begin() + thePolygonEdge[1]);
        }

        int thirdEdge = (theTriangleEdge[0]+1)%3;
        if (thirdEdge == theTriangleEdge[1])
            thirdEdge = (thirdEdge+1)%3;
        
        edgePoints.push_back(tmpEdgePoints[thirdEdge]);
        if (!reverse)
            std::reverse(edgePoints.back().begin(),edgePoints.back().end());
    }

    for (size_t i=0; i<edgePoints.size(); i++)
        for (int j=0; j<edgePoints[i].size(); j++)
            assert(edgePoints[i][j]>=0);
    
    garbageCollection();

    checkConsistency("Poly merge after garbage collection\n");

}

////////////////////////////////////////////////////////////////////////////////////
// from Graphics Gems III, page 199
// in our case we can safely assume that the denominator is not equal to zero
////////////////////////////////////////////////////////////////////////////////////

float DomainPolygon::computeIntersection(float &mu, const StaticVector<float,2> &p1, const StaticVector<float,2> &p2, const StaticVector<float,2> &p3, const StaticVector<float,2> &p4)
{

    StaticVector<float,2> A = p2 - p1;
    StaticVector<float,2> B = p3 - p4;
    StaticVector<float,2> C = p1 - p3;

    float det = A[1]*B[0] - A[0]*B[1];

    // if the determinant is zero we end up with a nan which is caught later
    //assert(det!=0);
    
    mu = (A[0]*C[1] - A[1]*C[0]) / det;
    return (B[1]*C[0] - B[0]*C[1]) / det;
}



bool DomainPolygon::triangulate(CircularPatch<float>& fillIn, std::vector<unsigned int>& nodeStack)
{
    int i, j, k;

    // debugging counter

    if (fillIn.size()>=2) {
        insertExtraEdges();
        createPointLocationStructure();
        removeExtraEdges();
    }

    for (i=0; i<nodes.size(); i++)
        assert(!std::isnan(nodes[i].domainPos()[0]) && !std::isnan(nodes[i].domainPos()[1]));

    for (i=0; i<fillIn.size()-1; i++){

        const int N = boundaryPoints.size();
        for (j=0; j<N; j++) {

            for (k=0; k<edgePoints[j].size()-1; k++){
                nodes[edgePoints[j][k]].setDomainEdge(j);
                nodes[edgePoints[j][k]].setDomainEdgePosition(k);
            }

        }

        DomainTriangle<float>& cT = par->triangles(fillIn[i]);

        //////////////////////////////////////////////////////////////////////////
        // cut off one triangle from the polygon and update the parametrization
        //////////////////////////////////////////////////////////////////////////

        // determine the planar world coordinates of the triangle about to be separated
        std::tr1::array<StaticVector<float,2>, 3> newTriangleCoords;
        int boundaryIdx;
        for (j=0; j<boundaryPoints.size(); j++)
            if (boundaryPoints[j] == par->triangles(fillIn[i]).vertices[0])
                boundaryIdx = j;
        
        for (j=0; j<3; j++)
            for (k=0; k<boundaryPoints.size(); k++)
                if (par->triangles(fillIn[i]).vertices[j] == boundaryPoints[k]){
                    newTriangleCoords[j] = nodes[cornerNode(k)].domainPos();
                    break;
                }


        cT.edgePoints[0] = edgePoints[boundaryIdx];
        cT.edgePoints[1] = edgePoints[(boundaryIdx+1)%boundaryPoints.size()];

        ///////////////////////////////////////////////////////////////////////////////////////
        // classify each node as either in the new triangle, on the separating segment
        // or in the remaining polygon
        ///////////////////////////////////////////////////////////////////////////////////////

        std::vector<int> nodeLocs(nodes.size());

        // edge nodes on polygon
        for (j=0; j<edgePoints.size(); j++)
            for (k=0; k<edgePoints[j].size(); k++)
                nodeLocs[edgePoints[j][k]] = IN_POLYGON;

        // edge nodes on triangle
        for (j=0; j<edgePoints[boundaryIdx].size(); j++)
            nodeLocs[edgePoints[boundaryIdx][j]] = IN_TRIANGLE;
        
        for (j=0; j<edgePoints[(boundaryIdx+1)%boundaryPoints.size()].size(); j++)
            nodeLocs[edgePoints[(boundaryIdx+1)%boundaryPoints.size()][j]] = IN_TRIANGLE;

        // two segment nodes
        nodeLocs[cornerNode(boundaryIdx)] = ON_SEGMENT;
        nodeLocs[edgePoints[(boundaryIdx+1)%boundaryPoints.size()].back()] = ON_SEGMENT;
        
        int cN;

        /////////////////////////////////////////////////////////////////
        // interior nodes.  here we have to make a geometric decision
        /////////////////////////////////////////////////////////////////
        for (cN=0; cN<nodes.size(); cN++) {

            if (nodes[cN].isINTERIOR_NODE()) {

                // This works only if the polygon is convex
                switch (orientation(newTriangleCoords[0], newTriangleCoords[2], nodes[cN].domainPos())) {
                case 1:
                    nodeLocs[cN] = IN_POLYGON;
                    break;
                case 0:
                    nodeLocs[cN] = ON_SEGMENT;
                    break;
                case -1:
                    nodeLocs[cN] = IN_TRIANGLE;
                    break;
                }

            }
        }
        
        // prepare the two new edgePoint lists
        std::vector<int> triNewEdgePoints;
        std::vector<int> polyNewEdgePoints;

        // ///////////////////////////////////////////////////////////////
        // transfer, and possibly cut, the parameterEdges
        // ///////////////////////////////////////////////////////////////
        cutParameterEdges(boundaryIdx, edgePoints[boundaryIdx][0],
                          edgePoints[(boundaryIdx+1)%boundaryPoints.size()].back(),
                          nodeLocs,
                          cT, newTriangleCoords, triNewEdgePoints, polyNewEdgePoints, nodeStack);
        
        // update boundaryPoints and edgePoints
        edgePoints[boundaryIdx] = polyNewEdgePoints;
        edgePoints.erase(edgePoints.begin() + (boundaryIdx+1)%edgePoints.size());
        
        boundaryPoints.erase(boundaryPoints.begin() + (boundaryIdx+1)%boundaryPoints.size());
        
        cT.edgePoints[2] = triNewEdgePoints;
        
        // sort out the nodes that belong onto the triangle
        int numTriNodes = 0;
        NodeIdx triNode;
        for (triNode=0; triNode<nodes.size(); triNode++) {
            if (nodeLocs[triNode] == IN_TRIANGLE)
                numTriNodes++;
            else if (nodeLocs[triNode] == ON_SEGMENT && 
                     nodeLocs[nodes[triNode].neighbors(0)] == IN_TRIANGLE) {
                    nodeLocs[triNode] = IN_TRIANGLE;
                    numTriNodes++;
            }
        }
        
        // debug
        assert(nodes.size()==nodeLocs.size());
        for (triNode=0; triNode<nodes.size(); triNode++) 
            if (nodeLocs[triNode] == IN_TRIANGLE)
                for (int triNode2=0; triNode2<nodes[triNode].degree(); triNode2++)
                    assert(nodeLocs[nodes[triNode].neighbors(triNode2)]==IN_TRIANGLE);
            
        
        int triCount = 0;
        cT.nodes.resize(numTriNodes);
        std::vector<int> offArr(nodes.size());
        
        for (triNode=0; triNode<nodes.size(); triNode++)
            if (nodeLocs[triNode] == IN_TRIANGLE) {
                cT.nodes[triCount] = nodes[triNode];
                invalidate(triNode);
                offArr[triNode] = triCount;
                triCount++;
            } else
                offArr[triNode] = -1;
        
        for (j=0; j<numTriNodes; j++)
            for (k=0; k<cT.nodes[j].degree(); k++)
                cT.nodes[j].neighbors(k) = offArr[cT.nodes[j].neighbors(k)];
        
        for (j=0; j<3; j++)
            for (k=0; k<cT.edgePoints[j].size(); k++){
                cT.edgePoints[j][k] = offArr[cT.edgePoints[j][k]];
            }


        ////////////////////////////////////////////////////////////////////////
        // Check whether two nonadjacent nodes in the new edgepoints array
        // are connected.  If yes abort
        ////////////////////////////////////////////////////////////////////////
        /** \todo what do I do with this? */
        for (j=0; j<cT.edgePoints[2].size()-2; j++) {

            const Node<float>& nA = cT.nodes[cT.edgePoints[2][j]];
            const Node<float>& nB = cT.nodes[cT.edgePoints[2][j+2]];

            if (!nA.isInvalid() && !nB.isInvalid() &&
                (nA.isConnectedTo(cT.edgePoints[2][j+2]) || 
                 nB.isConnectedTo(cT.edgePoints[2][j]))) {

                // This case may happen when smoothing horizontally.
                // It doesn't seem to be a programming error.
                assert(false);
            }
        }

        cT.installBarycentricCoordinates();

        // ward off more degeneracies in the input
        for (j=0; j<cT.nodes.size(); j++)
            if (std::isnan(cT.nodes[j].domainPos()[0]) || std::isnan(cT.nodes[j].domainPos()[1])) {
                printf("[DomainPolygon::triangulate]  Rejecting this vertex because\n");
                printf("  the domain position of one of its nodes contains NaNs.\n");
                printf("  This can be due to a degeneracy in the input data or to a\n");
                printf("  bug in the program.  Who knows?  But fear not:\n");
                printf("  Getting this message doesn't mean getting bad output.\n");
                return false;
            }

        cT.adjustTouchingNodes();

        checkConsistency("poly cut before GC\n");

        garbageCollection(offArr);
        for (j=0; j<offArr.size(); j++)
            nodeLocs[j-offArr[j]] = nodeLocs[j];


        cT.checkConsistency("recently cut\n");
        checkConsistency("poly recentlycut\n");
    }
    
    // there's only on triangle left
    DomainTriangle<float>& cT = par->triangles(fillIn.last());
    
    cT.nodes = nodes;

    cT.vertices[0] = boundaryPoints[0];
    cT.vertices[1] = boundaryPoints[1];
    cT.vertices[2] = boundaryPoints[2];
    
    cT.edgePoints[0] = edgePoints[0];
    cT.edgePoints[1] = edgePoints[1];
    cT.edgePoints[2] = edgePoints[2];



    cT.installBarycentricCoordinates();

    // ward off more degeneracies in the input
    for (j=0; j<cT.nodes.size(); j++)
        if (std::isnan(cT.nodes[j].domainPos()[0]) || std::isnan(cT.nodes[j].domainPos()[1])) {
            printf("[DomainPolygon::triangulate]  Rejecting this vertex because\n");
            printf("  the domain position of one of its nodes contains NaNs.\n");
            printf("  This can be due to a degeneracy in the input data or to a\n");
            printf("  bug in the program.  Who knows?  But fear not:\n");
            printf("  Getting this message doesn't mean getting bad output.\n");
            return false;
        }

    cT.adjustTouchingNodes();
    
    counter++;
    return true;
}

//#define PATHLOG

void DomainPolygon::cutParameterEdges(int boundaryIdx, NodeIdx startNode, NodeIdx lastNode,
                                      std::vector<int>& nodeLocs,
                                      DomainTriangle<float>& cT,
                                      const std::tr1::array<StaticVector<float,2>, 3>& newTriangleCoords,
                                      std::vector<int>& triNewEdgePoints,
                                      std::vector<int>& polyNewEdgePoints,
                                      std::vector<unsigned int>& nodeStack)
{
    int i;

    NodeIdx currentPolyNode = startNode;
    NodeIdx currentTriNode  = startNode;

    NodeIdx nextPolyNode = -1;
    NodeIdx nextTriNode  = -1;

    NodeIdx preNextPolyNode = -1;
    NodeIdx preNextTriNode  = -1;
    
    triNewEdgePoints.clear();
    polyNewEdgePoints.clear();

#ifdef PATHLOG
    FILE* fp = fopen("pathlog", "w");
    fprintf(fp, "# LineSet 0.1\n\nnDataVals 0\n\n{ ");
#endif

    while (currentPolyNode != lastNode) {

#ifdef PATHLOG
        fprintf(fp, "{\n");
        fprintf(fp, "%f %f 0\n", nodes[currentPolyNode].domainPos[0], nodes[currentPolyNode].domainPos[1]);
        fprintf(fp, "%f %f 0\n", nodes[currentTriNode].domainPos[0], nodes[currentTriNode].domainPos[1]);
        fprintf(fp, "} ");
#endif

        if (currentTriNode != currentPolyNode) {


            assert(nodes[currentPolyNode].isConnectedTo(currentTriNode));

            // //////////////////////////////////////
            // Find next pair of current###Nodes
            // //////////////////////////////////////

            // preNextPolyNode
            switch (nodes[currentPolyNode].type) {
            case Node<float>::INTERIOR_NODE: {
                DirectedEdgeIterator cPE = getDirectedEdgeIterator(currentPolyNode, currentTriNode);
                assert(cPE.isValid());

                DirectedEdgeIterator Onext = cPE.getONext();
                preNextPolyNode = Onext.to();
                break;
            }
            case Node<float>::TOUCHING_NODE:
            case Node<float>::CORNER_NODE:
                if (nodes[currentPolyNode].isLastNeighbor(currentTriNode)) {
                    preNextPolyNode = getPreviousEdgeNode(currentPolyNode);
                } else {

                    DirectedEdgeIterator cPE = getDirectedEdgeIterator(currentPolyNode, currentTriNode);
                    assert(cPE.isValid());

                    DirectedEdgeIterator Onext = cPE.getONext();
                    preNextPolyNode = Onext.to();
                }
                break;

            case Node<float>::INTERSECTION_NODE:
                preNextPolyNode = getPreviousEdgeNode(currentPolyNode);
                break;
            }

            
            // preNextTriNode
            switch (nodes[currentTriNode].type) {
            case Node<float>::INTERIOR_NODE: {
                DirectedEdgeIterator cPE = getDirectedEdgeIterator(currentPolyNode, currentTriNode);
                assert(cPE.isValid());

                DirectedEdgeIterator Dprev = cPE.getDPrev();
                preNextTriNode = Dprev.from();
                break;
            }
            case Node<float>::TOUCHING_NODE:
            case Node<float>::CORNER_NODE:
                if (nodes[currentTriNode].isFirstNeighbor(currentPolyNode)) {
                    
                    preNextTriNode = getNextEdgeNode(currentTriNode);
                } else {
                    DirectedEdgeIterator cPE = getDirectedEdgeIterator(currentPolyNode, currentTriNode);
                    assert(cPE.isValid());

                    DirectedEdgeIterator Dprev = cPE.getDPrev();
                    preNextTriNode = Dprev.from();
                }
                break;

            case Node<float>::INTERSECTION_NODE:
                preNextTriNode = getNextEdgeNode(currentTriNode);
                break;
            }

            nextPolyNode = currentPolyNode;
            nextTriNode = currentTriNode;

            if (nodeLocs[preNextPolyNode] == ON_SEGMENT) {

                nextPolyNode = nextTriNode = preNextPolyNode;

            } else if (nodeLocs[preNextTriNode] == ON_SEGMENT) {

                nextPolyNode = nextTriNode = preNextTriNode;

            } else if (nodeLocs[preNextPolyNode] == IN_POLYGON &&
                       nodeLocs[preNextTriNode] == IN_POLYGON) {

                    nextPolyNode = preNextTriNode;

            } else if (nodeLocs[preNextPolyNode] == IN_TRIANGLE &&
                       nodeLocs[preNextTriNode] == IN_POLYGON) {

                printf("shouldn't happen\n"); assert(false);

            } else if (nodeLocs[preNextPolyNode] == IN_POLYGON &&
                       nodeLocs[preNextTriNode] == IN_TRIANGLE) {

                if (nodes[preNextPolyNode].isConnectedTo(preNextTriNode)){
                    nextPolyNode = preNextPolyNode;
                    nextTriNode  = preNextTriNode;
                } else {

                    assert(nodes[preNextPolyNode].isINTERSECTION_NODE());
                    assert(nodes[preNextTriNode].isINTERSECTION_NODE());

                    if (nodes[preNextPolyNode].neighbors(0) == nodes[preNextTriNode].neighbors(0)) {

                        NodeIdx prePreNode = nodes[preNextPolyNode].neighbors(0);
                        
                        switch (nodeLocs[prePreNode]) {
                        case IN_POLYGON:
                            nextPolyNode = prePreNode;
                            nextTriNode  = preNextTriNode;
                            break;
                        case IN_TRIANGLE:
                            nextPolyNode = preNextPolyNode;
                            nextTriNode  = prePreNode;
                            break;
                        case ON_SEGMENT:
                            nextPolyNode = nextTriNode = prePreNode;
                            break;
                        }

                    } else if (nodes[getPreviousEdgeNode(preNextPolyNode)].isConnectedTo(preNextTriNode)) {
                        nextPolyNode = getPreviousEdgeNode(preNextPolyNode);
                        nextTriNode  = preNextTriNode;
                    } else if (nodes[getNextEdgeNode(preNextTriNode)].isConnectedTo(preNextPolyNode)) {
                        nextPolyNode = preNextPolyNode;
                        nextTriNode  = getNextEdgeNode(preNextTriNode);
                    } else {
                        NodeIdx prePrePolyNode = nodes[preNextPolyNode].neighbors(0);
                        NodeIdx prePreTriNode  = nodes[preNextTriNode].neighbors(0);

                        if (nodeLocs[prePrePolyNode] == IN_TRIANGLE) {
                            nextPolyNode = preNextPolyNode;
                            nextTriNode  = prePrePolyNode;
                        } else if (nodeLocs[prePrePolyNode] == ON_SEGMENT) {
                            nextPolyNode = nextTriNode = prePrePolyNode;
                        } else if (nodeLocs[prePreTriNode] == IN_POLYGON) {
                            nextPolyNode = prePreTriNode;
                            nextTriNode  = preNextTriNode;
                        } else if (nodeLocs[prePreTriNode] == ON_SEGMENT) {
                            nextPolyNode = nextTriNode = prePreTriNode;
                        } else {
                            printf("I don't know what else to do...\n");
                            assert(false);
                        }
                    }
                    

                }

                

            } else {

                nextTriNode = preNextPolyNode;

            }

            // /////////////////////////////////////////////////////////////////////
            // cut that edge and place the parts where they are supposed to be
            // /////////////////////////////////////////////////////////////////////
            float mu=0;
                    
            float lambda = computeIntersection(mu, nodes[currentTriNode].domainPos(),
                                         nodes[currentPolyNode].domainPos(),
                                         newTriangleCoords[2], newTriangleCoords[0]);
                    
            StaticVector<float,2> newDomainPos = linearInterpol(lambda, 
                                                  nodes[currentTriNode].domainPos(),
                                                  nodes[currentPolyNode].domainPos());  

            StaticVector<float,3> newImagePos  = linearInterpol(lambda,
                                                  par->iPos[nodes[currentTriNode].getNodeNumber()],
                                                  par->iPos[nodes[currentPolyNode].getNodeNumber()]);

            nodes.resize(nodes.size()+2);
            int newTriNode  = nodes.size()-2;
            int newPolyNode = nodes.size()-1;

            nodeLocs.resize(nodeLocs.size()+2);
            
            int newNodeNumber = createNodePosition(par->iPos, nodeStack, newImagePos);
            
            nodes[newTriNode].setValue(newDomainPos, newNodeNumber, Node<float>::INTERSECTION_NODE);
            nodes[newPolyNode].setValue(newDomainPos, newNodeNumber, Node<float>::INTERSECTION_NODE);
            
            nodeLocs[newTriNode]  = IN_TRIANGLE;
            nodeLocs[newPolyNode] = IN_POLYGON;
            
            // split edge           
            nodes[currentTriNode].replaceReferenceTo(currentPolyNode, newTriNode);
            nodes[currentPolyNode].replaceReferenceTo(currentTriNode, newPolyNode);
            
            nodes[newPolyNode].setNeighbor(currentPolyNode);
            nodes[newTriNode].setNeighbor(currentTriNode);
            
            // insert this node-pair into the edgepoints arrays
            polyNewEdgePoints.push_back(newPolyNode);
            triNewEdgePoints.push_back(newTriNode);


            
        } else {


            // //////////////////////////////////////
            // Find next pair of current###Nodes
            // //////////////////////////////////////

            const Node<float>& cN = nodes[currentPolyNode];
            assert(nodeLocs[currentPolyNode] == ON_SEGMENT);

            nextPolyNode = -1;

            for (i=0; i<cN.degree(); i++) {

                if (nodeLocs[cN.neighbors(i)] == ON_SEGMENT) {
                    nextPolyNode = nextTriNode = cN.neighbors(i);
                    break;
                }

                if (nodeLocs[cN.neighbors(i)] == IN_TRIANGLE &&
                    nodeLocs[cN.neighbors((i+1)%cN.degree())] == IN_POLYGON) {
                    nextTriNode = cN.neighbors(i);
                    nextPolyNode  = cN.neighbors((i+1)%cN.degree());
                    break;
                }
            }

            if (nextPolyNode == -1) {

                // all outgoing edges are either IN_TRIANGLE or IN_POLYGON
                assert(cN.isCORNER_NODE());
                if (cN.degree()==0) {

                    nextPolyNode = getPreviousEdgeNode(currentPolyNode);
                    nextTriNode  = getNextEdgeNode(currentTriNode);

                } else if (nodeLocs[cN.neighbors(0)] == IN_POLYGON) {

                    nextPolyNode = cN.neighbors(0);
                    if (nodes[nextPolyNode].isINTERSECTION_NODE()) {
                        nextPolyNode = getPreviousEdgeNode(nextPolyNode);
                    }
                    nextTriNode  = getNextEdgeNode(currentTriNode);

                } else {

                    nextPolyNode = getPreviousEdgeNode(currentPolyNode);
                    nextTriNode = cN.nbs.back();
                    if (nodes[nextTriNode].isINTERSECTION_NODE()) {
                        nextTriNode = getNextEdgeNode(nextTriNode);
                    }

                }

            } else {

                if (nodes[nextTriNode].isINTERSECTION_NODE()) {
                    nextTriNode = getNextEdgeNode(nextTriNode);
                }
                
                if (nodes[nextPolyNode].isINTERSECTION_NODE()) {
                    nextPolyNode = getPreviousEdgeNode(nextPolyNode);
                }

            }
                
            // ////////////////////////////////////////
            // split this node
            // ////////////////////////////////////////
        
            currentTriNode = splitNode(currentPolyNode, nodeLocs);
                
            /** \todo The second 'if' should be useless */
            if (cT.edgePoints[0][0] == currentPolyNode){
                cT.edgePoints[0][0] = currentTriNode;
            }else if (cT.edgePoints[1].back() == currentPolyNode){
                cT.edgePoints[1].back() = currentTriNode;
            }
            
            
            polyNewEdgePoints.push_back(currentPolyNode);
            triNewEdgePoints.push_back(currentTriNode);
        }
        
        
        currentPolyNode = nextPolyNode;
        currentTriNode  = nextTriNode;


    }

    // ////////////////////////////////////////
    // split the last corner node
    // ////////////////////////////////////////
    
    currentTriNode = splitNode(currentPolyNode, nodeLocs);
    
    /** \todo The first if should be useless */
    if (cT.edgePoints[0][0] == currentPolyNode){
        cT.edgePoints[0][0] = currentTriNode;
    }else if (cT.edgePoints[1].back() == currentPolyNode){ 
        cT.edgePoints[1].back() = currentTriNode;
    }

    polyNewEdgePoints.push_back(currentPolyNode);
    triNewEdgePoints.push_back(currentTriNode);

    std::reverse(triNewEdgePoints.begin(), triNewEdgePoints.end());

#ifdef PATHLOG
    fclose(fp);
#endif

}

NodeIdx DomainPolygon::splitNode(NodeIdx cN, std::vector<int>& nodeLocs)
{
    int i;

    assert(nodeLocs[cN]==ON_SEGMENT);

    nodes.resize(nodes.size()+1);
    NodeIdx newNode = nodes.size()-1;

    nodeLocs.resize(nodeLocs.size()+1);
          
    if (nodes[cN].isCORNER_NODE()){
        nodes[newNode].setValue(nodes[cN].domainPos(),
                                nodes[cN].getNodeNumber(), Node<float>::CORNER_NODE);
    }else{
        nodes[newNode].setValue(nodes[cN].domainPos(),
                                nodes[cN].getNodeNumber(), Node<float>::TOUCHING_NODE);
        nodes[cN].makeTouchingNode();
    }

    nodeLocs[newNode] = IN_TRIANGLE;
    nodeLocs[cN]      = IN_POLYGON;
            
    // ///////////////////////////////////////////////////////////
    // Readjust the neighborhood relations of the split nodes
    // ///////////////////////////////////////////////////////////

    int lastNonePolygonNeighbor;
    lastNonePolygonNeighbor = -1;
    for (i=0; i<nodes[cN].degree(); i++)
        if (nodeLocs[nodes[cN].neighbors(i)] != IN_POLYGON)
            lastNonePolygonNeighbor = i;

    if (lastNonePolygonNeighbor==-1)
        return newNode;

    
    DomainTriangle<float>::rotate(nodes[cN].nbs, -(lastNonePolygonNeighbor+1));

    
    for (i=nodes[cN].degree()-1; i>=0; i--){
        
        switch (nodeLocs[nodes[cN].neighbors(i)]){
        case  IN_POLYGON:
            // do nothing
            
            break;
            
        case  ON_SEGMENT: {
            // add edge from newNode to nodes[cN].neighbors(i)
            // The new edge has to appear topologically 'next to' the old one
            
            nodes[newNode].appendNeighbor(Node<float>::NeighborReference(nodes[cN].neighbors(i)));

            Node<float>& thisNeighbor = nodes[nodes[cN].neighbors(i)];
            for (int j=0; j<thisNeighbor.degree(); j++)
                if (thisNeighbor.neighbors(j)==cN){
                    thisNeighbor.nbs.insert(thisNeighbor.nbs.begin() + (j+1)%thisNeighbor.degree(), 
                                            Node<float>::NeighborReference(newNode));
                    break;
                }
            
            
            break;
        }
        case  IN_TRIANGLE:
            
            nodes[nodes[cN].neighbors(i)].replaceReferenceTo(cN, newNode);
            nodes[newNode].appendNeighbor(nodes[cN].neighbors(i));
            nodes[cN].removeNeighbor(i);
            
            break;
        }
        
    }
    
    return newNode;
}

unsigned int DomainPolygon::createNodePosition(std::vector<StaticVector<float,3> >& nodePositions, 
                                               std::vector<unsigned int>& nodeStack,
                                               const StaticVector<float,3>& newImagePos)
{

    if (nodeStack.size()!=0) {

        unsigned int newNodeNumber;
        newNodeNumber = nodeStack.back();
        nodeStack.pop_back();

        nodePositions[newNodeNumber] = newImagePos;

        return newNodeNumber;

    } else {

        nodePositions.push_back(newImagePos);
        return nodePositions.size()-1;

    }

}



void DomainPolygon::removeVertex(int point)
{
    int idx;
    int size = boundaryPoints.size();

    for (idx=0; idx<size; idx++)
        if (boundaryPoints[idx] == point)
            break;

    boundaryPoints.erase(boundaryPoints.begin() + idx);

    nodes[edgePoints[idx][0]].makeTouchingNode();

    edgePoints[(idx+size-1)%size].pop_back();

    for (int i=0; i<edgePoints[idx].size(); i++)
        edgePoints[(idx+size-1)%size].push_back(edgePoints[idx][i]);

    edgePoints.erase(edgePoints.begin() + idx);

}

void DomainPolygon::slice(int centerNode, int centerVertex, int bVertex)
{             
    std::vector<unsigned char> nodeLocs(nodes.size());

    nodes[centerNode].makeCornerNode();

    int i;

    int boundaryCutNode = cornerNode(bVertex);

    StaticVector<float,2> segmentFrom = nodes[centerNode].domainPos();
    StaticVector<float,2> segmentTo   = nodes[boundaryCutNode].domainPos();
        

    boundaryPoints.insert(boundaryPoints.begin() + bVertex, centerVertex);

    // die Null am Ende dieser Zeile ist höchstwahrscheinlich falsch!!
    boundaryPoints.insert(boundaryPoints.begin() + bVertex, boundaryPoints[0]);

    // prepare the two new edgePoint lists
    //edgePoints.insert(edgePoints.begin() + bVertex, 2, 0);
    edgePoints.insert(edgePoints.begin() + bVertex, std::vector<int>());
    edgePoints.insert(edgePoints.begin() + bVertex, std::vector<int>());

    std::vector<int>& newEdgePoints1 = edgePoints[bVertex];
    std::vector<int>& newEdgePoints2 = edgePoints[bVertex+1];

    newEdgePoints1.resize(2);
    newEdgePoints2.resize(2);
    
    std::vector<float> lambda1(2);
    
    newEdgePoints1[0] = boundaryCutNode;
    newEdgePoints1[1] = centerNode;
    
    lambda1[0] = 0;
    lambda1[1] = 1;
    
    newEdgePoints2[0] = centerNode;
    newEdgePoints2[1] = boundaryCutNode;
    
    
    // split nodes that lie on the cutting segment
    int newNode;

    // classify nodes
    for (size_t cN=0; cN<nodes.size(); cN++) {
        if (nodes[cN].isOnSegment(segmentFrom, segmentTo, 0.000001))
            
            if (nodes[cN].isINTERSECTION_NODE() || nodes[cN].isTOUCHING_NODE()){
                // rare type of numerical inconsistency
                printf("intersection node on segment!\n");
                assert(false);
            } else
                nodeLocs[cN] = ON_SEGMENT;
        else 
            switch(PlaneParam<float>::orientation(segmentFrom, segmentTo, nodes[cN].domainPos())){
            case  1:
            case  0:
                nodeLocs[cN] = IN_POLYGON;
                break;
            case  -1:
                nodeLocs[cN] = IN_TRIANGLE;
                break;
            }
    }

    // the two corner nodes that limit the cutting segment.  They have to get classified
    // ON_SEGMENT.  But sometimes they aren't due to numerical problems
    nodeLocs[centerNode]      = ON_SEGMENT;
    nodeLocs[boundaryCutNode] = ON_SEGMENT;


    for (int cN=nodes.size()-1; cN>=0; cN--) {

        if (nodeLocs[cN]==ON_SEGMENT && cN!=centerNode) {

            nodes.resize(nodes.size()+1);
            newNode = nodes.size()-1;
            nodeLocs.resize(nodeLocs.size()+1);

            nodes[newNode].setValue(nodes[cN].domainPos(), nodes[cN].getNodeNumber(), Node<float>::TOUCHING_NODE);

            if (nodes[cN].isCORNER_NODE())
                nodes[newNode].makeCornerNode();
            else {
                nodes[cN].makeTouchingNode();
                nodes[newNode].makeTouchingNode();
            }
            nodeLocs[cN]      = IN_POLYGON;
            nodeLocs[newNode] = IN_TRIANGLE;
            
            for (i=nodes[cN].degree()-1; i>=0; i--){
                
                switch(nodeLocs[nodes[cN].neighbors(i)]) {

                case  IN_POLYGON:
                    // do nothing
                    
                    break;
                    
                case ON_SEGMENT:                   

                    addEdge(newNode, nodes[cN].neighbors(i));

                    break;
                    
                case IN_TRIANGLE:

                    nodes[nodes[cN].neighbors(i)].replaceReferenceTo(cN, newNode);                  
                    nodes[newNode].appendNeighbor(nodes[cN].neighbors(i));
                    nodes[cN].removeNeighbor(i);
                    
                    break;
                }
                
            }
            
            // enter the split nodes in the edgePoints arrays
            
            if (cN == boundaryCutNode) {
                
                edgePoints[(bVertex+edgePoints.size()-1)%edgePoints.size()].back() = newNode;
                newEdgePoints1[0] = newNode;
                
            }else{
                float lambda = (nodes[newNode].domainPos() - segmentTo).length() /
                    (segmentFrom - segmentTo).length();

                int j=0;
                while (lambda1[j]<lambda) j++;
                
                lambda1.insert(lambda1.begin() + j, lambda);
                newEdgePoints1.insert(newEdgePoints1.begin()+j, newNode);
                
                newEdgePoints2.insert(newEdgePoints2.begin()+newEdgePoints2.size()-j, cN);
                
            }
            
        }
        
    }

    for (i=0; i<newEdgePoints1.size(); i++){
        nodeLocs[newEdgePoints1[i]] = ON_SEGMENT;
        nodeLocs[newEdgePoints2[i]] = ON_SEGMENT;
    }

    // bisect edges that cross the cut
    
    int polyNode;

    for (size_t triNode=0; triNode<nodes.size(); triNode++) {

        if (nodeLocs[triNode] != IN_TRIANGLE)
            continue;

        for (i=0; i<nodes[triNode].degree(); i++) {

            if (nodeLocs[nodes[triNode].neighbors(i)] == IN_POLYGON) {
                    
                polyNode = nodes[triNode].neighbors(i);
    
                float mu=0;
                float lambda;

        
                lambda = computeIntersection(mu, nodes[triNode].domainPos(), nodes[polyNode].domainPos(), 
                                             segmentTo, segmentFrom);
        
                if (lambda<=0.0f || lambda>=1.0f || mu<=0.0f || mu>=1.0f)
                    continue;
            

                StaticVector<float,2> newDomainPos = nodes[triNode].domainPos() + 
                    lambda*(nodes[polyNode].domainPos() - nodes[triNode].domainPos());
                StaticVector<float,3> newImagePos  = par->iPos[nodes[triNode].getNodeNumber()]  + 
                    lambda*(par->iPos[nodes[polyNode].getNodeNumber()]  - 
                             par->iPos[nodes[triNode].getNodeNumber()]);

                nodes.resize(nodes.size()+2);
                int newTriNode  = nodes.size()-2;
                int newPolyNode = nodes.size()-1;

                nodeLocs.resize(nodeLocs.size()+2);

                par->iPos.push_back(newImagePos);
                int newNodeNumber = par->iPos.size()-1;

                nodes[newTriNode].setValue(newDomainPos, newNodeNumber, Node<float>::INTERSECTION_NODE);
                nodes[newPolyNode].setValue(newDomainPos, newNodeNumber, Node<float>::INTERSECTION_NODE);
                
                nodeLocs[newTriNode]  = ON_SEGMENT;
                nodeLocs[newPolyNode] = ON_SEGMENT;
        
                nodes[triNode].replaceReferenceTo(polyNode, newTriNode);
                nodes[polyNode].replaceReferenceTo(triNode, newPolyNode);
        
                nodes[newTriNode].appendNeighbor(triNode);
                nodes[newPolyNode].appendNeighbor(polyNode);

                // insert new nodes in edgePoint arrays
                int j=0;
                while (lambda1[j]<mu) j++;
        
                lambda1.insert(lambda1.begin() + j, mu);
                newEdgePoints1.insert(newEdgePoints1.begin()+j, newTriNode);
        
                newEdgePoints2.insert(newEdgePoints2.begin()+newEdgePoints2.size()-j, newPolyNode);
        
        
            }

        }

    }

}


void DomainPolygon::garbageCollection(std::vector<int>& offArr)
{
    int offset = 0;

    offArr.resize(nodes.size());

    for (size_t i=0; i<nodes.size(); i++){
        offArr[i] = offset;

        if (nodes[i].isInvalid()) 
            offset++;
    }

    ////////////////////
    for (size_t i=0; i<offArr.size(); i++)
        nodes[i-offArr[i]] = nodes[i];

    nodes.resize(nodes.size()-offset);

    ///////////////////
    for (size_t i=0; i<nodes.size(); i++)
        for (int j=0; j<nodes[i].degree(); j++)
            nodes[i].neighbors(j) -= offArr[nodes[i].neighbors(j)];

    //////////////////
    for (size_t i=0; i<edgePoints.size(); i++)
        for (int j=0; j<edgePoints[i].size(); j++)
            edgePoints[i][j] -= offArr[edgePoints[i][j]];

}


void DomainPolygon::createPointLocationStructure() 
{
    checkConsistency("BeforeCreate");

    for (size_t i=0; i<nodes.size(); i++){
        
        if (nodes[i].isINTERIOR_NODE()){
            makeCyclicInteriorNode(nodes[i]);
        }
    }

    checkConsistency("AfterInterior");

    const int N = boundaryPoints.size();

    for (int i=0; i<N; i++) {
        
        checkConsistency("Edge");

        makeCyclicBoundaryNode(nodes[edgePoints[i][0]], 
                               edgePoints[i][1], 
                               edgePoints[(i+N-1)%N][edgePoints[(i+N-1)%N].size()-2]);

        checkConsistency("AfterCorners");

        for (int j=1; j<edgePoints[i].size()-1; j++){
            makeCyclicBoundaryNode(nodes[edgePoints[i][j]], edgePoints[i][j+1], edgePoints[i][j-1]);

            if (nodes[edgePoints[i][j]].isINTERSECTION_NODE()){
                nodes[edgePoints[i][j]].setDomainEdge(i);
                nodes[edgePoints[i][j]].setDomainEdgePosition(j);
            }
        }

        checkConsistency("AfterEdges");

    }

}


void DomainPolygon::insertExtraEdges()
{
    int i,j;

    const int N = boundaryPoints.size();

    for (i=0; i<N; i++){
        
        for (j=1; j<edgePoints[i].size(); j++){
            
            if (nodes[edgePoints[i][j]].isINTERSECTION_NODE() ||
                nodes[edgePoints[i][j-1]].isINTERSECTION_NODE()){

                
                addEdge(edgePoints[i][j-1], edgePoints[i][j], true);
            }
        }

    }

    for (i=0; i<N; i++){
    
        for (j=1; j<edgePoints[i].size()-1; j++)
            
            if (nodes[edgePoints[i][j]].isINTERSECTION_NODE()){

                int interiorPoint = nodes[edgePoints[i][j]].neighbors(0);


                if (!nodes[interiorPoint].isConnectedTo(edgePoints[i][j-1])){

                    addEdge(edgePoints[i][j-1], interiorPoint, true);
                }
            }   
    }
}


void DomainPolygon::print(bool showEdgePoints, bool showParamEdges, bool showNodes) const 
{
    printf("--------------------------------------------------------\n");
    printf("--  Print Polygon  -------------------------------------\n");

//     printf("points:  ");
//     for (i=0; i<boundaryPoints.size(); i++)
//      printf("%d  ", (int)boundaryPoints[i]);

    printf("\n");


    if (showEdgePoints){
        
        for (size_t i=0; i<edgePoints.size(); i++){
            printf("edgePoints %d:\n", i);
            for (size_t j=0; j<edgePoints[i].size(); j++){
                printf("  %d) -- ", edgePoints[i][j]);
                nodes[edgePoints[i][j]].print();
            }
        }
        
        printf("\n");
    }

    if (showNodes){
        for (int cN=0; cN<nodes.size(); cN++){
            printf("%d  ", cN);
            nodes[cN].print(showParamEdges);
        }
    }

    printf("---------------------------------------------------------\n\n");
}




void DomainPolygon::checkConsistency(const char* where)
{
#ifndef NDEBUG
    
    int i, j;
    PlaneParam<float>::checkConsistency(where);

    // check whether all corner nodes are of type CORNER_NODE
    for (i=0; i<edgePoints.size(); i++){
//      if (edgePoints[i][0]->type != Node::CORNER_NODE || 
//          edgePoints[i].last()->type != Node::CORNER_NODE){
//          printf(where);
//          printf("***** corner node is not CORNER_NODE *****\n");
//          assert(FALSE);
//      }
        // check if two subsequent TOUCHING_NODES are connected by an edge
        for (j=0; j<edgePoints[i].size()-1; j++){

            //assert(edgePoints[i][j]>=0 && edgePoints[i][j+1]>=0);
            if (edgePoints[i][j]<0 || edgePoints[i][j+1]<0) {
//              printf(where);
//              printf("***** negative edgePoints ********\n");
//              assert(false);
                continue;
            }

            Node<float>& nA = nodes[edgePoints[i][j]];
            Node<float>& nB = nodes[edgePoints[i][j+1]];

            if (!nA.isInvalid() && !nB.isInvalid() && 
                nA.isTOUCHING_NODE() && nB.isTOUCHING_NODE() &&
                !nA.isConnectedTo(edgePoints[i][j+1])){
                printf(where);
                printf("***** two adjacent TOUCHING NODES are not connected! *****\n");
                assert(false);
            }
        }

        // check whether nodes that are not neighbors in the edgePoint array are connected
        for (j=0; j<edgePoints[i].size()-2; j++) {
            if (edgePoints[i][j]<0 || edgePoints[i][j+2]<0)
                continue;

            const Node<float>& nA = nodes[edgePoints[i][j]];
            const Node<float>& nB = nodes[edgePoints[i][j+2]];
            
            if (!nA.isInvalid() && !nB.isInvalid() &&
                (nA.isConnectedTo(edgePoints[i][j+2]) || nB.isConnectedTo(edgePoints[i][j]))) {

                printf(where);
                printf("%d  and  %d\n", edgePoints[i][j], edgePoints[i][j+2]);
                printf("****** two nonadjacent nodes are connected!! *******\n");
                assert(false);
            }
            
        }
    }
#endif
}

/////////////////////////////////////////////////////////////////
//
// $Id: Parametrization.cpp,v 1.2 2007/12/17 10:13:39 sander Exp $
// 
// $Log: Parametrization.cpp,v $
//
/////////////////////////////////////////////////////////////////

#include <tr1/array>
#include <vector>

#if defined HAVE_AMIRAMESH || !defined PSURFACE_STANDALONE
#include <amiramesh/AmiraMesh.h>
#endif

#include <psurface/PSurface.h>
#include <psurface/GlobalNodeIdx.h>

#if defined HAVE_AMIRAMESH || !defined PSURFACE_STANDALONE
template <int dim, class ctype>
PSurface<dim,ctype>::PSurface(HxParamBundle* bundle)
{
    if (bundle) {
        params = bundle;
        hasOwnParamBundle = false;
    } else {
        params = new HxParamBundle;
        hasOwnParamBundle = true;
    }
}
#else
template <int dim, class ctype>
PSurface<dim,ctype>::PSurface()
{}
#endif

template <int dim, class ctype>
PSurface<dim,ctype>::~PSurface() 
{ 
#if defined HAVE_AMIRAMESH || !defined PSURFACE_STANDALONE
    if (hasOwnParamBundle)
        delete params;
#endif
}

template <int dim, class ctype>
void PSurface<dim,ctype>::clear()
{
    surface = NULL;
    patches.clear();

#if defined HAVE_AMIRAMESH || !defined PSURFACE_STANDALONE
    if (hasOwnParamBundle)
        delete params;
#endif

    iPos.clear();
    paths.clear();
    McSurfaceBase<McVertex<float>, McEdge, DomainTriangle>::clear();
}

template <int dim, class ctype>
void PSurface<dim,ctype>::getBoundingBox(Box<float,3>& bbox) const
{
    if (getNumVertices()==0)
        return;

    bbox.set(vertices(0), vertices(0));

    for (int i=1; i<getNumVertices(); i++)
        bbox.extendBy(vertices(i));
}

template <int dim, class ctype>
void PSurface<dim,ctype>::init(const PSurface* other)
{
    // copy domain surface.  
    *this = *other;

    surface = new Surface();
    *surface = *other->surface;
    
}

template <int dim, class ctype>
StaticVector<float,2> PSurface<dim,ctype>::getLocalTargetCoords(const GlobalNodeIdx& n, int targetTri) const
{
    const Node<float>& cN = triangles(n.tri).nodes[n.idx];
    
    switch (cN.type) {
    case Node<float>::GHOST_NODE:
    case Node<float>::INTERSECTION_NODE: {

        StaticVector<float,3> iPos = imagePos(n.tri, n.idx);
        return triangles(n.tri).computeBarycentricCoords(iPos, *(StaticVector<float,3>*)&surface->points[surface->triangles[targetTri].points[0]][0], 
                                                         *(StaticVector<float,3>*)&surface->points[surface->triangles[targetTri].points[1]][0], 
                                                         *(StaticVector<float,3>*)&surface->points[surface->triangles[targetTri].points[2]][0]);
    }
    default:
        if (cN.getNodeNumber()==surface->triangles[targetTri].points[0])
            return StaticVector<float,2>(1, 0);
        else if (cN.getNodeNumber()==surface->triangles[targetTri].points[1])
            return StaticVector<float,2>(0, 1);
        else if (cN.getNodeNumber()==surface->triangles[targetTri].points[2])
            return StaticVector<float,2>(0, 0);
        else {
            printf("The node is not related to the targetTri!\n");
            throw ParamError();
        }
    }
}


template <int dim, class ctype>
GlobalNodeIdx PSurface<dim,ctype>::getOtherEndNode(int triIdx, NodeIdx cN) const
{
    int i;
    
    while (triangles(triIdx).nodes[cN].isINTERSECTION_NODE()) {
        
        const DomainTriangle& cT = triangles(triIdx);
        
        // get the edge the node is on and its position in the edgePoints array
        int edge  = cT.nodes[cN].getDomainEdge();
        int edgePos = cT.nodes[cN].getDomainEdgePosition();
        
        // get adjacent triangle
        const int cE = cT.getOppositeEdge(cT.vertices[(edge+2)%3]);

#ifndef NDEBUG
        if (edges(cE).numTriangles()!=2) {
            printf("Edge:  %d --> %d\n", edges(cE).from, edges(cE).to);
            for (i=0; i<edges(cE).numTriangles(); i++)
                triangles(edges(cE).triangles[i]).print(true, true, true);
        }
#endif

        assert(edges(cE).numTriangles()==2);
        
        const int oppT = (edges(cE).triangles[0]==triIdx) 
            ? edges(cE).triangles[1] 
            : edges(cE).triangles[0];
        
        // get the opposite edgePoint array
        int oppEdge = -1;
        bool reverse = false;
        for (i=0; i<3; i++) {
            if (triangles(oppT).vertices[i] == cT.vertices[edge] && 
                triangles(oppT).vertices[(i+1)%3]==cT.vertices[(edge+1)%3]) {
                
                oppEdge = i;
                reverse = false;
                break;
            } else if (triangles(oppT).vertices[i] == cT.vertices[(edge+1)%3] && 
                       triangles(oppT).vertices[(i+1)%3] == cT.vertices[edge]) {
                
                oppEdge = i;
                reverse = true;
                break;
            }
        }

        assert(oppEdge!=-1);
        
        int oppEdgePos = (reverse) ? cT.edgePoints[edge].size()-edgePos-1 : edgePos;

        if (triangles(oppT).nodes[triangles(oppT).edgePoints[oppEdge][oppEdgePos]].getNodeNumber()
            != cT.nodes[cT.edgePoints[edge][edgePos]].getNodeNumber()) {

            printf("Condition triangles(oppT).nodes[triangles(oppT).edgePoints[oppEdge][oppEdgePos]].getNodeNumber() != cT.nodes[cT.edgePoints[edge][edgePos]].getNodeNumber()  failed!\n");
            
            throw ParamError();

        }
        
        int newPoint = triangles(oppT).nodes[triangles(oppT).edgePoints[oppEdge][oppEdgePos]].theInteriorNode();
        
        // get the new triangle and node
        triIdx = oppT;
        cN = newPoint;
        
    }
    
    return GlobalNodeIdx(triIdx, cN);
}

template <int dim, class ctype>
int PSurface<dim,ctype>::getNumNodes() const
{
    int n = 0;
    for (int i=0; i<getNumTriangles(); i++)
        n += triangles(i).nodes.size();
    return n;
}

template <int dim, class ctype>
int PSurface<dim,ctype>::getNumTrueNodes() 
{
    int highestTrueNodeNumber = -1;

    for (int j(0); j<getNumTriangles(); j++) {

        const DomainTriangle& cT = triangles(j);

        for (int i=0; i<cT.nodes.size(); i++){
            int a = cT.nodes[i].getNodeNumber();
            if (!cT.nodes[i].isINTERSECTION_NODE() && (a-highestTrueNodeNumber)>0){
                highestTrueNodeNumber = a;
            }
        }
    }
    
    return highestTrueNodeNumber+1;
}

template <int dim, class ctype>
void PSurface<dim,ctype>::removeExtraEdges()
{

    for (int i(0); i<getNumTriangles(); i++) {
        triangles(i).removeExtraEdges();
    }

    hasUpToDatePointLocationStructure = false;
}

template <int dim, class ctype>
void PSurface<dim,ctype>::createPointLocationStructure()
{
    for (int i(0); i<getNumTriangles(); i++){
        //        printf("######## Triangle: %d\n", i);
        triangles(i).checkConsistency("Before Insert");
        triangles(i).insertExtraEdges();
//         if (i==95)
//             triangles(i).print(false, true, true);
        triangles(i).createPointLocationStructure();
    }

    hasUpToDatePointLocationStructure = true;
}

    /// \todo The copying could be sped up considerably...
template <int dim, class ctype>
void PSurface<dim,ctype>::garbageCollection()
{
    int i, j;
#ifndef NDEBUG
    printf("This is the Parametrization garbage collection...\n");
    std::cout << freeVertexStack.size() << " vertices, "
              << freeEdgeStack.size()   << " edges, "
              << freeTriangleStack.size() << " triangles removed" << std::endl;
#endif
    
    std::vector<bool> isInvalid;
    
    // clean up vertices
    if (freeVertexStack.size()) {
        
        int offset = 0;
        
        std::vector<int> vertexOffsets(vertexArray.size());
        isInvalid.resize(vertexArray.size());
        for (i=0; i<isInvalid.size(); i++)
            isInvalid[i] = false;
        
        for (i=0; i<freeVertexStack.size(); i++)
            isInvalid[freeVertexStack[i]] = true;
        
        for (i=0; i<vertexArray.size(); i++){
            vertexOffsets[i] = offset;
            
            if (isInvalid[i]) 
                offset++;
        }
        
        ////////////////////
        for (i=0; i<vertexOffsets.size(); i++)
            vertexArray[i-vertexOffsets[i]] = vertexArray[i];
        
        vertexArray.resize(vertexArray.size()-offset);
        
        // Adjust edges
        for (i=0; i<edgeArray.size(); i++) {
            edgeArray[i].from -= vertexOffsets[edgeArray[i].from];
            edgeArray[i].to   -= vertexOffsets[edgeArray[i].to];
        }                
        
        // Adjust triangles
        for (i=0; i<triangleArray.size(); i++) 
            for (j=0; j<3; j++) 
                triangleArray[i].vertices[j] -= vertexOffsets[triangleArray[i].vertices[j]];
        
        // Adjust paths
        for (i=0; i<paths.size(); i++)
            for (j=0; j<paths[i].points.size(); j++)
                paths[i].points[j] -= vertexOffsets[paths[i].points[j]];

        freeVertexStack.clear();                

    }

    McSurfaceBase<McVertex<float>, McEdge, DomainTriangle>::garbageCollection();
    
#ifndef NDEBUG
        printf("   ...Garbage collection finished!\n");
#endif
    }
    

#if defined HAVE_AMIRAMESH || !defined PSURFACE_STANDALONE
template <int dim, class ctype>
int PSurface<dim,ctype>::writeAmiraMesh(PSurface<dim,ctype>* par, const char* filename)
{
    AmiraMesh am;
    am.parameters = *(par->params);

    int i, j, k;
    int numVertices  = par->getNumVertices();
    int numTriangles = par->getNumTriangles();

    am.parameters.set("ContentType", "Parametrization");

    par->savePaths(am.parameters);
    

    ////////////////////////////////////////////
    // the patches array
    AmiraMesh::Location* patchesLoc = new AmiraMesh::Location("Patches", par->numPatches());
    am.insert(patchesLoc);

    AmiraMesh::Data* patchesData = new AmiraMesh::Data("Patches", patchesLoc, 
                                                       McPrimType::mc_int32, 3, 
                                                       (void*)&par->patches[0]);
    am.insert(patchesData);

    //////////////////////////////////////
    // store base grid vertices
    AmiraMesh::Location* vertices = new AmiraMesh::Location("BaseGridVertexCoords", numVertices);
    am.insert(vertices);

    std::vector<StaticVector<float,3> > baseGridVertexCoordsArray(numVertices);
    for (i=0; i<par->getNumVertices(); i++)
        baseGridVertexCoordsArray[i] = par->vertices(i);

    AmiraMesh::Data* vertexCoords = new AmiraMesh::Data("BaseGridVertexCoords", vertices, 
                                                        McPrimType::mc_float, 3, 
                                                        (void*)&baseGridVertexCoordsArray[0]);
    am.insert(vertexCoords);

    /////////////////////////////////////
    // store base grid triangles
    AmiraMesh::Location* triangles = new AmiraMesh::Location("BaseGridTriangles", numTriangles);
    am.insert(triangles);

    std::vector<std::tr1::array<int, 3> > baseGridTriArray(numTriangles);

    for (int i(0); i<par->getNumTriangles(); i++)
        baseGridTriArray[i] = par->triangles(i).vertices;
    

    AmiraMesh::Data* triangleCoords = new AmiraMesh::Data("BaseGridTriangles", triangles,
                                                          McPrimType::mc_int32, 3, 
                                                          (void*)&baseGridTriArray);
    am.insert(triangleCoords);

    ////////////////////////////////////////
    // the image positions of the nodes
    AmiraMesh::Location* nodePosLoc = new AmiraMesh::Location("NodePositions", par->iPos.size());
    am.insert(nodePosLoc);

    AmiraMesh::Data* nodePosData = new AmiraMesh::Data("NodePositions", nodePosLoc,
                                                       McPrimType::mc_float, 3, 
                                                       (void*)&par->iPos[0]);
    am.insert(nodePosData);

    ////////////////////////////////////////////////////////////
    // the number of nodes and parameterEdges per triangle and the
    // patch it belongs to.  Plus the NodeNumber of the three corner nodes
    //
    AmiraMesh::Location* numNodesAndEdges = new AmiraMesh::Location("NumNodesAndParameterEdgesPerTriangle", numTriangles);
    am.insert(numNodesAndEdges);

    int numNodes      = 0;
    int numParamEdges = 0;
    int numEdgePoints = 0;

    std::vector<int> numNodesAndEdgesArray(11*numTriangles);
    
    for (i=0; i<numTriangles; i++) {

        const DomainTriangle& cT = par->triangles(i);

        int numIntersectionNodes;
        int numTouchingNodes;
        int numInteriorNodes;

        cT.countNodes(numIntersectionNodes, numTouchingNodes, numInteriorNodes);
        int numEdges = cT.getNumRegularEdges();

        numNodesAndEdgesArray[11*i+0] = numIntersectionNodes;
        numNodesAndEdgesArray[11*i+1] = numTouchingNodes;
        numNodesAndEdgesArray[11*i+2] = numInteriorNodes;
        numNodesAndEdgesArray[11*i+3] = numEdges;
        numNodesAndEdgesArray[11*i+4] = cT.patch;

        numNodesAndEdgesArray[11*i+5] = cT.edgePoints[0].size()-2;
        numNodesAndEdgesArray[11*i+6] = cT.edgePoints[1].size()-2;
        numNodesAndEdgesArray[11*i+7] = cT.edgePoints[2].size()-2;

        numNodesAndEdgesArray[11*i+8] = cT.nodes[cT.cornerNode(0)].getNodeNumber();
        numNodesAndEdgesArray[11*i+9] = cT.nodes[cT.cornerNode(1)].getNodeNumber();
        numNodesAndEdgesArray[11*i+10] = cT.nodes[cT.cornerNode(2)].getNodeNumber();

        numNodes += numIntersectionNodes;
        numNodes += numTouchingNodes;
        numNodes += numInteriorNodes;

        numEdgePoints += cT.edgePoints[0].size() + cT.edgePoints[1].size() + cT.edgePoints[2].size() - 6;

        numParamEdges += numEdges;

    }

    AmiraMesh::Data*  numNodesAndEdgesData = new AmiraMesh::Data("NumNodesAndParameterEdgesPerTriangle", 
                                                                 numNodesAndEdges,
                                                                 McPrimType::mc_int32, 11, 
                                                                 (void*)&numNodesAndEdgesArray[0]);
    am.insert(numNodesAndEdgesData);
    
    /////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////

    std::vector<StaticVector<float,2> > domainPositions(numNodes);
    std::vector<int>     nodeNumbers(numNodes);
    std::vector<std::tr1::array<int,2> > parameterEdgeArray(numParamEdges);
    std::vector<int>     edgePointsArray(numEdgePoints);

    int cN;    

    int arrayIdx           = 0;
    int edgeArrayIdx       = 0;
    int edgePointsArrayIdx = 0;

    for (i=0; i<numTriangles; i++) {

        const DomainTriangle& cT = par->triangles(i);

        std::vector<int> newIdx(cT.nodes.size());
        int localArrayIdx = 3;
        // the cornerNode are not saved, because everything about them
        // can be deduced from the base grid
        newIdx[cT.cornerNode(0)] = 0;
        newIdx[cT.cornerNode(1)] = 1;
        newIdx[cT.cornerNode(2)] = 2;

        // the three remaining types are saved separately, in order to avoid
        // explicitly saving the type for each node.
        for (cN=0; cN<cT.nodes.size(); cN++) {
            if (cT.nodes[cN].isINTERSECTION_NODE()){

                domainPositions[arrayIdx] = cT.nodes[cN].domainPos();
                nodeNumbers[arrayIdx]     = cT.nodes[cN].getNodeNumber();

                newIdx[cN] = localArrayIdx;
                
                arrayIdx++;
                localArrayIdx++;
            }
        }

        for (cN=0; cN<cT.nodes.size(); cN++) {

            if (cT.nodes[cN].isTOUCHING_NODE()){

                domainPositions[arrayIdx] = cT.nodes[cN].domainPos();
                nodeNumbers[arrayIdx]     = cT.nodes[cN].getNodeNumber();

                newIdx[cN] = localArrayIdx;

                arrayIdx++;
                localArrayIdx++;
            }
        }

        for (cN=0; cN<cT.nodes.size(); cN++) {

            if (cT.nodes[cN].isINTERIOR_NODE()){

                domainPositions[arrayIdx] = cT.nodes[cN].domainPos();
                nodeNumbers[arrayIdx]     = cT.nodes[cN].getNodeNumber();

                newIdx[cN] = localArrayIdx;
                
                arrayIdx++;
                localArrayIdx++;
            }
        }
        
        /////////////////////////////////////////
        // the parameterEdges for this triangle

        PlaneParam<float>::UndirectedEdgeIterator cE;
        for (cE = cT.firstUndirectedEdge(); cE.isValid(); ++cE){
            if (cE.isRegularEdge()) {
                parameterEdgeArray[edgeArrayIdx][0] = newIdx[cE.from()];
                parameterEdgeArray[edgeArrayIdx][1] = newIdx[cE.to()];
                edgeArrayIdx++;
            }
        }


        ///////////////////////////////////////
        // the edgePoints for this triangle
        for (j=0; j<3; j++){

            for (k=1; k<cT.edgePoints[j].size()-1; k++)
                edgePointsArray[edgePointsArrayIdx++] = newIdx[cT.edgePoints[j][k]];

        }
    }

    AmiraMesh::Location* nodeLoc = new AmiraMesh::Location("Nodes", numNodes);
    am.insert(nodeLoc);

    AmiraMesh::Data*  nodeData   = new AmiraMesh::Data("Nodes", nodeLoc,
                                                       McPrimType::mc_float, 2, (void*)&domainPositions[0]);
    am.insert(nodeData);

    AmiraMesh::Location* nNLoc = new AmiraMesh::Location("NodeNumbers", numNodes);
    am.insert(nNLoc);

    AmiraMesh::Data* nNData    = new AmiraMesh::Data("NodeNumbers", nNLoc,
                                                     McPrimType::mc_int32, 1, (void*)&nodeNumbers[0]);
    am.insert(nNData);

    
    AmiraMesh::Location* parameterEdges = new AmiraMesh::Location("ParameterEdges", numParamEdges);
    am.insert(parameterEdges);

    AmiraMesh::Data*  parameterEdgeData = new AmiraMesh::Data("ParameterEdges", parameterEdges,
                                                              McPrimType::mc_int32, 2, (void*)&parameterEdgeArray[0]);
    am.insert(parameterEdgeData);

    /////////////////////////////////////////
    // The edgePoints arrays

    AmiraMesh::Location* edgePoints = new AmiraMesh::Location("EdgePoints", edgePointsArray.size());
    am.insert(edgePoints);

    AmiraMesh::Data* edgePointsData = new AmiraMesh::Data("EdgePoints", edgePoints,
                                                          McPrimType::mc_int32, 1, 
                                                          (void*)&edgePointsArray[0]);
    am.insert(edgePointsData);

    //////////////////////////////
    // actually write the file
    if (!am.write(filename) ) {
        printf("An error has occured writing file %s.\n", filename);
        return 0;
    }

    return 1;

}

template <int dim, class ctype>
void* PSurface<dim,ctype>::readAmiraMesh(AmiraMesh* am, const char* filename)
{
    PSurface<dim,ctype>* par = new PSurface<dim,ctype>;

    Surface* surf = new Surface;

    if (!par->initFromAmiraMesh(am, filename, surf)) {
        delete par;
        delete surf;

        return NULL;
    }

    *par->params = am->parameters;

    return par;
}

template <int dim, class ctype>
bool PSurface<dim,ctype>::initFromAmiraMesh(AmiraMesh* am, const char* filename, Surface* surf)
{

    surface = surf;
    
    getPaths(am->parameters);

    ///////////////////////////////////////////////
    // test for file format
    AmiraMesh::Data* AMnodePos = am->findData("NodePositions", HxFLOAT, 3, "NodePositions");

    if (!AMnodePos){
        printf("You're trying to read a file in the old-style AmiraMesh format.\n");
        printf("But the format is not supported anymore :-(\n");
        return false;
    }

    int i, j;

    //////////////////////////////////////////////
    // the innerRegions and outerRegions arrays
    AmiraMesh::Data* AMpatches = am->findData("Patches", HxINT32, 3, "Patches");
    if (!AMpatches){
        printf("AmiraMesh: Field 'Patches' not found!\n");
        return false;
    }
    
    patches.resize(AMpatches->location()->dims()[0]);

    for (i=0; i<patches.size(); i++){
        patches[i] = ((Patch*)AMpatches->dataPtr())[i];
    }
    
    ///////////////////////////////
    // load the base grid vertices
    AmiraMesh::Data* AMvertices = am->findData("BaseGridVertexCoords", HxFLOAT, 3, "BaseGridVertexCoords");
    if (!AMvertices){
        printf("AmiraMesh: Field 'BaseGridVertexCoords' not found!\n");
        return false;
    }
    
    int numPoints = AMvertices->location()->dims()[0];

    const StaticVector<float,3>* vertexCoords = (StaticVector<float,3>*)AMvertices->dataPtr();

    // copy points
    for (i=0; i<numPoints; i++)
        newVertex(vertexCoords[i]);

    ///////////////////////////////////////////////////////
    // copy node positions
    iPos.resize(AMnodePos->location()->dims()[0]);
    
    for (i=0; i<iPos.size(); i++)
        iPos[i] = ((StaticVector<float,3>*)AMnodePos->dataPtr())[i];
    
    ////////////////////////////////////////////////////
    // copy triangles.  This takes care of the edges, too
    AmiraMesh::Data* AMtriangles = am->findData("BaseGridTriangles", HxINT32, 3, "BaseGridTriangles");
    if (!AMtriangles){
        printf("AmiraMesh: Field 'BaseGridTriangles' not found!\n");
        return false;
    }
    
    const std::tr1::array<int,3>* triIdx = (std::tr1::array<int,3>*)AMtriangles->dataPtr();
    const int numTriangles = AMtriangles->location()->dims()[0];

    AmiraMesh::Data* numNodesAndEdges = am->findData("NumNodesAndParameterEdgesPerTriangle", HxINT32, 11, 
                                                     "NumNodesAndParameterEdgesPerTriangle");
    if (!numNodesAndEdges){
        printf("AmiraMesh: Field 'NumNodesAndParameterEdgesPerTriangle' not found!\n");
        return false;
    }
    AmiraMesh::Data* nodes            = am->findData("Nodes", HxFLOAT, 2, "Nodes");
    if (!nodes){
        printf("AmiraMesh: Field 'Nodes' not found!\n");
        return false;
    }
    
    AmiraMesh::Data* AMnodeNumbers    = am->findData("NodeNumbers", HxINT32, 1, "NodeNumbers");
    if (!AMnodeNumbers){
        printf("AmiraMesh: Field 'NodeNumbers' not found!\n");
        return false;
    }
    
    AmiraMesh::Data* AMedges            = am->findData("ParameterEdges", HxINT32, 2, "ParameterEdges");
    if (!AMedges){
        printf("AmiraMesh: Field 'ParameterEdges' not found!\n");
        return false;
    }
    AmiraMesh::Data* AMedgePoints      = am->findData("EdgePoints", HxINT32, 1, "EdgePoints");
    if (!AMedgePoints){
        printf("AmiraMesh: Field 'EdgePoints' not found!\n");
        return false;
    }

    const int* numNodesAndEdgesData = (int*)numNodesAndEdges->dataPtr();
    const StaticVector<float,2>*   nodeData       = (StaticVector<float,2>*)nodes->dataPtr();
    const int* nodeNumbers          = (int*)AMnodeNumbers->dataPtr();
    const std::tr1::array<int,2>* edgeData  = (std::tr1::array<int,2>*)AMedges->dataPtr();
    const int*     edgePointData    = (int*)AMedgePoints->dataPtr();


    int edgeCounter=0, edgePointCounter=0;
    int nodeArrayIdx = 0;
    
    for (i=0; i<numTriangles; i++){
        
        int newTriIdx = createSpaceForTriangle(triIdx[i][0], triIdx[i][1], triIdx[i][2]);

        triangles(newTriIdx).patch = numNodesAndEdgesData[11*i+4];

        ///////////////////////////////////////////////
        // get the parametrization on this triangle

        int numIntersectionNodes = numNodesAndEdgesData[11*i+0];
        int numTouchingNodes     = numNodesAndEdgesData[11*i+1];
        int numInteriorNodes     = numNodesAndEdgesData[11*i+2];

        int numParamEdges        = numNodesAndEdgesData[11*i+3];

        ////////////////////////////////////
        // first the nodes

        int cornerNodes[3];

        triangles(newTriIdx).nodes.resize(numIntersectionNodes + numTouchingNodes + numInteriorNodes + 3);

        // three corner nodes
        StaticVector<float,2> domainPos(1, 0);
        int nodeNumber = numNodesAndEdgesData[11*i+8];
        triangles(newTriIdx).nodes[0].setValue(domainPos, nodeNumber, Node<float>::CORNER_NODE);
        cornerNodes[0] = 0;

        domainPos = StaticVector<float,2>(0, 1);
        nodeNumber = numNodesAndEdgesData[11*i+9];
        triangles(newTriIdx).nodes[1].setValue(domainPos, nodeNumber, Node<float>::CORNER_NODE);
        cornerNodes[1] = 1;

        domainPos = StaticVector<float,2>(0, 0);
        nodeNumber = numNodesAndEdgesData[11*i+10];
        triangles(newTriIdx).nodes[2].setValue(domainPos, nodeNumber, Node<float>::CORNER_NODE);
        cornerNodes[2] = 2;

        int nodeCounter = 3;

        // the intersection nodes
        for (j=0; j<numIntersectionNodes; j++, nodeCounter++, nodeArrayIdx++){
            StaticVector<float,2> domainPos = nodeData[nodeArrayIdx];
            int nodeNumber    = nodeNumbers[nodeArrayIdx];

            triangles(newTriIdx).nodes[nodeCounter].setValue(domainPos, nodeNumber, Node<float>::INTERSECTION_NODE);
        }

        // the touching nodes
        for (j=0; j<numTouchingNodes; j++, nodeCounter++, nodeArrayIdx++){
            StaticVector<float,2> domainPos = nodeData[nodeArrayIdx];
            int nodeNumber    = nodeNumbers[nodeArrayIdx];
            
            triangles(newTriIdx).nodes[nodeCounter].setValue(domainPos, nodeNumber, Node<float>::TOUCHING_NODE);
        }

        // the interior nodes
        for (j=0; j<numInteriorNodes; j++, nodeCounter++, nodeArrayIdx++){
            StaticVector<float,2> domainPos = nodeData[nodeArrayIdx];
            int nodeNumber    = nodeNumbers[nodeArrayIdx];

            triangles(newTriIdx).nodes[nodeCounter].setValue(domainPos, nodeNumber, Node<float>::INTERIOR_NODE);
        }

        ///////////////////////////////
        // the parameterEdges
        for (j=0; j<numParamEdges; j++, edgeCounter++){
            triangles(newTriIdx).addEdge(edgeData[edgeCounter][0], edgeData[edgeCounter][1]);
        }

        //////////////////////////////////
        // the edgePoints arrays
        for (j=0; j<3; j++){

            triangles(newTriIdx).edgePoints[j].resize(numNodesAndEdgesData[11*i+5+j] + 2);

            triangles(newTriIdx).edgePoints[j][0]     = j;
            triangles(newTriIdx).edgePoints[j].back() = (j+1)%3;

            for (int k=0; k<numNodesAndEdgesData[11*i+5+j]; k++){
                triangles(newTriIdx).edgePoints[j][k+1] = edgePointData[edgePointCounter];
                edgePointCounter++;
            }

        }

        integrateTriangle(newTriIdx);
    }
    
    // sad but true ...
    hasUpToDatePointLocationStructure = false;

    setupOriginalSurface();

    return true;
}
#endif


template <int dim, class ctype>
void PSurface<dim,ctype>::setupOriginalSurface()
{
    int i, j, k;
    
    if (!hasUpToDatePointLocationStructure)
        createPointLocationStructure();
    
    // transfer materials

    

#ifndef PSURFACE_STANDALONE
    // create patches.  Only the Amira Surface class needs this
    surface->patches.resize(numPatches());

    for (i=0; i<surface->patches.size(); i++){
        surface->patches[i] = new Surface::Patch();
        surface->patches[i]->innerRegion = patches[i].innerRegion;
        surface->patches[i]->outerRegion = patches[i].outerRegion;
        surface->patches[i]->boundaryId  = patches[i].boundaryId;
    }
#endif
    ////////////////////////////////////////////
    //
    surface->points.resize(getNumTrueNodes());
    for (i=0; i<surface->points.size(); i++)
        for (int j=0; j<3; j++)
            surface->points[i][j] = iPos[i][j];

    ////////////////////////////////////////////
    //

    for (k=0; k<getNumTriangles(); k++) {
        //printf("i = %d von %d\n", i, getNumTriangles());
        DomainTriangle& cT = triangles(k);

        ////////////////////////////////
        int numNodes = cT.nodes.size();

        for (i=0; i<numNodes; i++) {

            Node<float>& cN = cT.nodes[i];
            std::tr1::array<int,3> v;

            v[0] = cN.nodeNumber;

            switch (cN.type) {

            case Node<float>::INTERSECTION_NODE:
                continue;

            case Node<float>::INTERIOR_NODE:
                
                for (j=0; j<cN.degree(); j++) {
                    
                    if (!cN.neighbors(j).isRegular())
                        continue;
                    
//                     printf("A\n");
//                     nodes(GlobalNodeIdx(k, cN.neighbors(j))).print();
                    //v[1] = getOtherEndNodeNumber(k, cN.neighbors(j));
                    v[1] = nodes(getOtherEndNode(k, cN.neighbors(j))).getNodeNumber();

                    int nN = (j+1)%cN.degree();
                    
                    if (!cN.neighbors(nN).isRegular())
                        nN = (nN+1)%cN.degree();
                    
//                     printf("B\n");
//                     nodes(GlobalNodeIdx(k, cN.neighbors(nN))).print();
                    //v[2] = getOtherEndNodeNumber(k, cN.neighbors(nN));
                    v[2] = nodes(getOtherEndNode(k, cN.neighbors(nN))).getNodeNumber();

                    // insert triangle
                    if (v[0] < v[1] && v[0] < v[2]) {
                        appendTriangleToOriginalSurface(v, cT.patch);
                    }
                }
                
                break;

            case Node<float>::TOUCHING_NODE:
            case Node<float>::CORNER_NODE:
                

                int firstRegNeighbor = -1;

                for (j=0; j<cN.degree(); j++)
                    if (cN.neighbors(j).isRegular()) {
                        firstRegNeighbor = j;
                        break;
                    }

                // corner node has no regular neighbors  --> do nothing
                if (firstRegNeighbor == -1)
                    break;

                if (firstRegNeighbor != 0) {

//                     printf("C\n");
//                     nodes(GlobalNodeIdx(k, cN.neighbors(0))).print();
                    v[1] = nodes(getOtherEndNode(k, cN.neighbors(0))).getNodeNumber();
//                     printf("D\n");
//                     nodes(GlobalNodeIdx(k, cN.neighbors(firstRegNeighbor))).print();
                    v[2] = nodes(getOtherEndNode(k, cN.neighbors(firstRegNeighbor))).getNodeNumber();

                    // insert triangle
                    if (v[0] < v[1] && v[0] < v[2]) {
                        appendTriangleToOriginalSurface(v, cT.patch);
                    }
                }

                // /////////////////////////////////////////
                int nextRegNeighbor;
                do {
                    for (nextRegNeighbor=firstRegNeighbor+1; 
                         nextRegNeighbor<cN.degree(); nextRegNeighbor++)
                        if (cN.neighbors(nextRegNeighbor).isRegular()) 
                            break;
                                        
                    if (nextRegNeighbor <cN.degree()) {
                    
                        v[1] = nodes(getOtherEndNode(k, cN.neighbors(firstRegNeighbor))).getNodeNumber();
                        v[2] = nodes(getOtherEndNode(k, cN.neighbors(nextRegNeighbor))).getNodeNumber();

                        // insert triangle
                        if (v[0] < v[1] && v[0] < v[2]) {
                            appendTriangleToOriginalSurface(v, cT.patch);
                        }
                    }

                    firstRegNeighbor = nextRegNeighbor;

                } while (nextRegNeighbor<cN.degree());

                // ////////////////////////
                if (!cN.neighbors(cN.degree()-1).isRegular()) {

#if 0
                    int otherTri;
                    v[1] = v[2];
                    v[2] = getOtherEndNodeNumber(k, cN.neighbors(cN.degree()-1), &otherTri);
                    
                    // insert triangle
                    if (v[0] < v[1] && v[0] < v[2] && k<otherTri) {
                        
                        appendTriangleToOriginalSurface(v, cT.patch);
                    }
#endif
                }

                break;
            }
        }
        
    }
}

template <int dim, class ctype>
void PSurface<dim,ctype>::appendTriangleToOriginalSurface(const std::tr1::array<int,3>& v, int patch)
{
    surface->triangles.push_back(Surface::Triangle());
                        
    surface->triangles.back().points[0] = v[0];
    surface->triangles.back().points[1] = v[1];
    surface->triangles.back().points[2] = v[2];

#ifndef PSURFACE_STANDALONE    
    // The Amira Surface class needs a consistent 'patch' structure
    surface->triangles.back().patch = patch;
    surface->patches[patch]->triangles.push_back(surface->triangles.size()-1);
#endif
}


#if defined HAVE_AMIRAMESH || !defined PSURFACE_STANDALONE
template <int dim, class ctype>
void PSurface<dim,ctype>::getPaths(const HxParamBundle& parameters)
{
    int i;
    paths.resize(0);

    for (i=0 ; i<parameters.size() ; i++) {
        if (strncmp("SurfacePath", parameters[i]->name(),11)==0) {
            HxParameter* p = (HxParameter*) parameters[i];
            if (!p->isBundle() && p->primType()== McPrimType::mc_int32) {
                paths.resize(paths.size()+1);
                for (int j=0 ; j<p->dim() ; j+=2) {
                    paths.back().points.push_back(p->getNum(j));
                    paths.back().isFix.push_back(p->getNum(j+1));
                }
            }
        }
    }
}

template <int dim, class ctype>
void PSurface<dim,ctype>::savePaths(HxParamBundle& parameters)
{
    int i;

    for (i=0 ; i<parameters.size() ; i++) 
        if (strncmp("SurfacePath", parameters[i]->name(), 11)==0)
            parameters.remove(parameters[i--]);

    for ( i=0 ; i<paths.size() ; i++) {
        char buf[64];
        sprintf(buf,"SurfacePath%03d",i);
        parameters.remove(buf);
        std::vector<int> tmp;
        for (int j=0 ; j<paths[i].points.size() ; j++) {
            tmp.push_back(paths[i].points[j]);
            tmp.push_back(paths[i].isFix[j]);
        }
        HxParameter* p = new HxParameter(buf, tmp.size(), &tmp[0]);
        parameters.insert(p);
    }
}
#endif


template <int dim, class ctype>
int PSurface<dim,ctype>::map(int triIdx, StaticVector<float,2>& p, std::tr1::array<int,3>& vertices, 
                         StaticVector<float,2>& coords, int seed) const
{
    int i;
    const DomainTriangle& tri = triangles(triIdx);
    const std::vector<StaticVector<float,3> >& nP = iPos;

    // this is boundary handling
    if (p[0] < 0.001){
        for (i=0; i<tri.edgePoints[1].size()-1; i++){

            const StaticVector<float,2>& a = tri.nodes[tri.edgePoints[1][i]].domainPos();
            const StaticVector<float,2>& b = tri.nodes[tri.edgePoints[1][i+1]].domainPos();

            if (a[1]+1e-5 > p[1] && p[1] > b[1]-1e-5) {

                std::tr1::array<GlobalNodeIdx, 3> targetNodes;
                handleMapOnEdge(triIdx, p, a, b, 1, i, targetNodes, coords);
                vertices[0] = nodes(targetNodes[0]).getNodeNumber();
                vertices[1] = nodes(targetNodes[1]).getNodeNumber();
                vertices[2] = nodes(targetNodes[2]).getNodeNumber();
                return true;
            }
            
        }

    }else if (p[1] < 0.001){

        for (i=0; i<tri.edgePoints[2].size()-1; i++){
            const StaticVector<float,2>& a = tri.nodes[tri.edgePoints[2][i]].domainPos();
            const StaticVector<float,2>& b = tri.nodes[tri.edgePoints[2][i+1]].domainPos();

            if (b[0]+1e-5 > p[0] && p[0] > a[0]-1e-5) {

                std::tr1::array<GlobalNodeIdx, 3> targetNodes;
                handleMapOnEdge(triIdx, p, a, b, 2, i, targetNodes, coords);
                vertices[0] = nodes(targetNodes[0]).getNodeNumber();
                vertices[1] = nodes(targetNodes[1]).getNodeNumber();
                vertices[2] = nodes(targetNodes[2]).getNodeNumber();
                return true;
            }
        }

    }else if (p[0]+p[1] > 0.999) {

        for (i=0; i<tri.edgePoints[0].size()-1; i++){
            const StaticVector<float,2>& a = tri.nodes[tri.edgePoints[0][i]].domainPos();
            const StaticVector<float,2>& b = tri.nodes[tri.edgePoints[0][i+1]].domainPos();

            if (a[0]+1e-5>p[0] && p[0]>b[0]-1e-5) {

                std::tr1::array<GlobalNodeIdx, 3> targetNodes;
                handleMapOnEdge(triIdx, p, a, b, 0, i, targetNodes, coords);
                vertices[0] = nodes(targetNodes[0]).getNodeNumber();
                vertices[1] = nodes(targetNodes[1]).getNodeNumber();
                vertices[2] = nodes(targetNodes[2]).getNodeNumber();
                return true;

            }
        }

        return false;

    }

    std::tr1::array<NodeIdx, 3> v;
    int status = tri.map(p, v, coords, seed);

    if (!status)
        return 0;

    StaticVector<float,3> imagePos = PlaneParam<float>::linearInterpol<StaticVector<float,3> >(coords, nP[tri.nodes[v[0]].getNodeNumber()],
                                                           nP[tri.nodes[v[1]].getNodeNumber()],
                                                           nP[tri.nodes[v[2]].getNodeNumber()]);

    // ///////////////////////////////////////////////////////
    // make sure we don't return intersection nodes
    std::tr1::array<GlobalNodeIdx, 3> resultNodes;
    getActualVertices(triIdx, v, resultNodes);
    vertices[0] = nodes(resultNodes[0]).getNodeNumber();
    vertices[1] = nodes(resultNodes[1]).getNodeNumber();
    vertices[2] = nodes(resultNodes[2]).getNodeNumber();

    coords = tri.computeBarycentricCoords(imagePos, nP[vertices[0]], nP[vertices[1]], nP[vertices[2]]);

    return true;
}

template <int dim, class ctype>
void PSurface<dim,ctype>::getActualVertices(int tri, const std::tr1::array<NodeIdx, 3>& nds,
                                        std::tr1::array<GlobalNodeIdx, 3>& vertices) const
{
    const DomainTriangle& cT = triangles(tri);
    //cT.print(true, true, true);
    int mode = cT.nodes[nds[0]].isINTERSECTION_NODE() +
        2*cT.nodes[nds[1]].isINTERSECTION_NODE() +
        4*cT.nodes[nds[2]].isINTERSECTION_NODE();
    //printf("MODE %d \n", mode);
    for (int i=0; i<3; i++)
        vertices[i] = getOtherEndNode(tri, nds[i]);
    //printf("***************MODE %d \n", mode);
    if (mode==6) {

        if (triangles(vertices[1].tri).nodes[vertices[1].idx].getNodeNumber() == 
            triangles(vertices[2].tri).nodes[vertices[2].idx].getNodeNumber()) {

            int int1 = cT.nodes[nds[1]].theInteriorNode();
            int int2 = cT.nodes[nds[2]].theInteriorNode();
            assert(int1!=int2);
            assert(int1==nds[0] || int2==nds[0]);
            
            if (int1 == nds[0])
                vertices[2] = getOtherEndNode(tri, int2);
            else
                vertices[1] = getOtherEndNode(tri, int1);

        }

    } else if (mode==5) {

        if (nodes(vertices[0]).getNodeNumber() == nodes(vertices[2]).getNodeNumber()) {

            int int0 = cT.nodes[nds[0]].theInteriorNode();
            int int2 = cT.nodes[nds[2]].theInteriorNode();
            assert(int0!=int2);
            assert(int0==nds[1] || int2==nds[1]);
            
            if (int0 == nds[1])
                vertices[2] = getOtherEndNode(tri, int2);
            else
                vertices[0] = getOtherEndNode(tri, int0);

        }

    } else if (mode==3) {

        if (triangles(vertices[1].tri).nodes[vertices[1].idx].getNodeNumber() == 
            triangles(vertices[0].tri).nodes[vertices[0].idx].getNodeNumber()) {

            int int1 = cT.nodes[nds[1]].theInteriorNode();
            int int0 = cT.nodes[nds[0]].theInteriorNode();
            assert(int1!=int0);
            assert(int1==nds[2] || int0==nds[2]);
            
            if (int1 == nds[2])
                vertices[0] = getOtherEndNode(tri, int0);
            else
                vertices[1] = getOtherEndNode(tri, int1);

        }

    } else {
        // They are all three intersection nodes...

        if (nodes(vertices[1]).getNodeNumber() == nodes(vertices[0]).getNodeNumber() &&
            cT.nodes[nds[1]].isINTERSECTION_NODE() &&
            cT.nodes[nds[0]].isINTERSECTION_NODE()) {

            NodeIdx int1 = cT.nodes[nds[1]].theInteriorNode();
            NodeIdx int0 = cT.nodes[nds[0]].theInteriorNode();
            assert(int1!=int0);
            assert(int1==nds[2] || int0==nds[2]);
            
            if (int1 == nds[2])
                vertices[0] = getOtherEndNode(tri, int0);
            else
                vertices[1] = getOtherEndNode(tri, int1);

        } else if (nodes(vertices[1]).getNodeNumber() == nodes(vertices[2]).getNodeNumber() &&
                   cT.nodes[nds[1]].isINTERSECTION_NODE() &&
                   cT.nodes[nds[2]].isINTERSECTION_NODE()) {

            int int1 = cT.nodes[nds[1]].theInteriorNode();
            int int2 = cT.nodes[nds[2]].theInteriorNode();
            assert(int1!=int2);
            assert(int1==nds[0] || int2==nds[0]);
            
            if (int1 == nds[0])
                vertices[2] = getOtherEndNode(tri, int2);
            else
                vertices[1] = getOtherEndNode(tri, int1);

        } else if (nodes(vertices[2]).getNodeNumber() == nodes(vertices[0]).getNodeNumber() &&
                   cT.nodes[nds[2]].isINTERSECTION_NODE() &&
                   cT.nodes[nds[0]].isINTERSECTION_NODE()) {

            int int0 = cT.nodes[nds[0]].theInteriorNode();
            int int2 = cT.nodes[nds[2]].theInteriorNode();
            assert(int0!=int2);
            assert(int0==nds[1] || int2==nds[1]);
            
            if (int0 == nds[1])
                vertices[2] = getOtherEndNode(tri, int2);
            else
                vertices[0] = getOtherEndNode(tri, int0);

        }

    }
}

template <int dim, class ctype>
int PSurface<dim,ctype>::getImageSurfaceTriangle(int tri,
                                             const std::tr1::array<NodeIdx, 3>& nds
                                             ) const
{
    int i;
    
    std::tr1::array<GlobalNodeIdx, 3> actualVertices;
    std::tr1::array<std::vector<int>, 3> trianglesPerNode;

    getActualVertices(tri, nds, actualVertices);
    
    assert(surface->trianglesPerPoint.size());
    
    for (i=0; i<3; i++)
        trianglesPerNode[i] = getTargetTrianglesPerNode(actualVertices[i]);
    
    for (i=0; i<trianglesPerNode[0].size(); i++) {

#ifdef PSURFACE_STANDALONE
        if (std::find(trianglesPerNode[1].begin(), 
                      trianglesPerNode[1].end(), 
                      trianglesPerNode[0][i]) != trianglesPerNode[1].end() &&
            std::find(trianglesPerNode[2].begin(),
                      trianglesPerNode[2].end(),
                      trianglesPerNode[0][i]) != trianglesPerNode[2].end())
            return trianglesPerNode[0][i];
#else
        if (mcSmallArray::index(trianglesPerNode[1], trianglesPerNode[0][i])!=-1 &&
            mcSmallArray::index(trianglesPerNode[2], trianglesPerNode[0][i])!=-1)
            return trianglesPerNode[0][i];
#endif

    }
            
    return -1;
}

template <int dim, class ctype>
std::vector<int> PSurface<dim,ctype>::getTargetTrianglesPerNode(const GlobalNodeIdx& n) const
{
    assert(surface->trianglesPerPoint.size());
    const Node<float>& cN = triangles(n.tri).nodes[n.idx];
    const float eps = 1e-6;

    switch (cN.type) {
    case Node<float>::GHOST_NODE: {
        
        std::vector<int> result(1);
        result[0] = cN.getNodeNumber();
        if (cN.dP[0] + cN.dP[1] > 1-eps) {
            // append the triangles bordering on edge 0
            int p = surface->triangles[result[0]].points[0];
            int q = surface->triangles[result[0]].points[1];

            getTrianglesPerEdge(p, q, result, result[0]);

        } else if (cN.dP[0] < eps) {
            // append the triangles bordering on edge 1
            int p = surface->triangles[result[0]].points[1];
            int q = surface->triangles[result[0]].points[2];

            getTrianglesPerEdge(p, q, result, result[0]);

        } else if (cN.dP[1] < eps) {
            // append the triangles bordering on edge 2
            int p = surface->triangles[result[0]].points[2];
            int q = surface->triangles[result[0]].points[0];

            getTrianglesPerEdge(p, q, result, result[0]);
            
        }
        
        return result;
        
        
    }
    case Node<float>::INTERSECTION_NODE:
        assert(false);
    }

    // Copying from a McSmallVector to a std::vector
    std::vector<int> result(surface->trianglesPerPoint[cN.getNodeNumber()].size());
    for (int i=0; i<result.size(); i++)
        result[i] = surface->trianglesPerPoint[cN.getNodeNumber()][i];

    return result;
    
}

/// This is a service routine only for getTargetTrianglesPerNode
template <int dim, class ctype>
void PSurface<dim,ctype>::getTrianglesPerEdge(int from, int to, std::vector<int>& tris, int exception) const
{
    for (int i=0; i<surface->trianglesPerPoint[from].size(); i++) {

#ifdef PSURFACE_STANDALONE
        if (std::find(surface->trianglesPerPoint[to].begin(),
                      surface->trianglesPerPoint[to].end(),
                      surface->trianglesPerPoint[from][i]) != surface->trianglesPerPoint[to].end() &&
            surface->trianglesPerPoint[from][i] != exception)

#else
        if (mcSmallArray::index(surface->trianglesPerPoint[to], surface->trianglesPerPoint[from][i]) != -1 &&
            surface->trianglesPerPoint[from][i] != exception)
#endif
            tris.push_back(surface->trianglesPerPoint[from][i]);

    }

}

template <int dim, class ctype>
void PSurface<dim,ctype>::handleMapOnEdge(int triIdx, const StaticVector<float,2>& p, const StaticVector<float,2>& a, const StaticVector<float,2>& b,
                                      int edge, int edgePos, std::tr1::array<GlobalNodeIdx, 3>& vertices, StaticVector<float,2>& coords) const
{
    const DomainTriangle& tri = triangles(triIdx);
    float lambda = (p-a).length() / (a-b).length();

    StaticVector<float,3> targetPos = PlaneParam<float>::linearInterpol<StaticVector<float,3> >(lambda, 
                                                            imagePos(triIdx, tri.edgePoints[edge][edgePos]),
                                                            imagePos(triIdx, tri.edgePoints[edge][edgePos+1]));

    int n1 = tri.edgePoints[edge][edgePos];
    int n2 = tri.edgePoints[edge][edgePos+1];

    vertices[0] = getOtherEndNode(triIdx, n1);
    vertices[1] = getOtherEndNode(triIdx, n2);

    ///////////////////////////////////////////////
    if (tri.nodes[n1].isINTERSECTION_NODE() && tri.nodes[n2].isINTERSECTION_NODE()) {
        
        int intNode1 = tri.nodes[n1].theInteriorNode();
        int intNode2 = tri.nodes[n2].theInteriorNode();
        
        int intNodeNumber1 = nodes(getOtherEndNode(triIdx, intNode1)).getNodeNumber();
        int intNodeNumber2 = nodes(getOtherEndNode(triIdx, intNode2)).getNodeNumber();

        assert(nodes(vertices[0]).getNodeNumber() != nodes(vertices[1]).getNodeNumber() || 
               intNodeNumber1!=intNodeNumber2);
        assert(nodes(vertices[0]).getNodeNumber() == nodes(vertices[1]).getNodeNumber() || 
               intNodeNumber1==intNodeNumber2);
        
        if (intNodeNumber1==intNodeNumber2){
            vertices[2] = getOtherEndNode(triIdx, intNode1);
        } else {
            vertices[1] = getOtherEndNode(triIdx, intNode2);
            vertices[2] = getOtherEndNode(triIdx, intNode1);
        }

    } else if (!tri.nodes[n1].isINTERSECTION_NODE() && tri.nodes[n2].isINTERSECTION_NODE()) {
        
        vertices[2] = getOtherEndNode(triIdx, tri.nodes[n2].theInteriorNode());

    } else if  (tri.nodes[n1].isINTERSECTION_NODE() && !tri.nodes[n2].isINTERSECTION_NODE()) {
        
        vertices[2] = getOtherEndNode(triIdx, tri.nodes[n1].theInteriorNode());

    } else {
        
        PlaneParam<float>::DirectedEdgeIterator edge = tri.getDirectedEdgeIterator(n1, n2);
        vertices[2] = getOtherEndNode(triIdx, edge.getONext().to());

    }

    assert(nodes(vertices[0]).getNodeNumber() != nodes(vertices[1]).getNodeNumber() && 
           nodes(vertices[1]).getNodeNumber() != nodes(vertices[2]).getNodeNumber() && 
           nodes(vertices[2]).getNodeNumber() != nodes(vertices[0]).getNodeNumber());
    
    coords = tri.computeBarycentricCoords(targetPos, imagePos(vertices[0]), 
                                          imagePos(vertices[1]), 
                                          imagePos(vertices[2]));

}


template <int dim, class ctype>
int PSurface<dim,ctype>::positionMap(int triIdx, StaticVector<float,2>& p, StaticVector<float,3>& result) const
{
    StaticVector<float,2> localCoords;
    std::tr1::array<int,3> tri;

    int status = map(triIdx, p, tri, localCoords);
    
    if (!status) {
        printf("p: (%f %f)\n", p[0], p[1]);
        triangles(triIdx).print(true, true, false);
        assert(false);
        return false;
    }

    result = PlaneParam<float>::linearInterpol<StaticVector<float,3> >(localCoords, iPos[tri[0]], iPos[tri[1]], iPos[tri[2]]);

    return true;
}



template <int dim, class ctype>
int PSurface<dim,ctype>::directNormalMap(int triIdx, StaticVector<float,2>& p, StaticVector<float,3>& result) const
{
    StaticVector<float,2> localCoords;
    std::tr1::array<int,3> tri;

    int status = map(triIdx, p, tri, localCoords);
    
    if (!status)
        return false;

    const StaticVector<float,3> a = iPos[tri[1]] - iPos[tri[0]];
    const StaticVector<float,3> b = iPos[tri[2]] - iPos[tri[0]];
    result = a.cross(b);
    result.normalize();

    assert(!isnan(result[0]) && !isnan(result[1]) && !isnan(result[2]));

    return true;
}

template <int dim, class ctype>
int PSurface<dim,ctype>::invertTriangles(int patch)
{
    
    int i;
    int count=0;
    
    for (i=int(0); i<int(getNumTriangles()); i++) 
        if (patch==-1 || triangles(i).patch==patch){
            
            triangles(i).flip();
            count++;
            
            if (hasUpToDatePointLocationStructure) {
                
                for (int j=0; j<triangles(i).nodes.size(); j++)
                    triangles(i).nodes[j].reverseNeighbors();
            }
            
        }
    
    return count;
}

template <int dim, class ctype>
NodeIdx PSurface<dim,ctype>::addInteriorNode(int tri, const StaticVector<float,2>& dom, int nodeNumber)
{
    triangles(tri).nodes.push_back(Node<float>(dom, nodeNumber, Node<float>::INTERIOR_NODE));
    return triangles(tri).nodes.size()-1;
}

template <int dim, class ctype>
NodeIdx PSurface<dim,ctype>::addGhostNode(int tri, int corner, int targetTri, const StaticVector<float,2>& localTargetCoords)
{
    triangles(tri).nodes.push_back(Node<float>());
    triangles(tri).nodes.back().makeGhostNode(corner, targetTri, localTargetCoords);
    return triangles(tri).nodes.size()-1;
}

template <int dim, class ctype>
NodeIdx PSurface<dim,ctype>::addCornerNode(int tri, int corner, int nodeNumber)
{
    DomainTriangle& cT = triangles(tri);

    cT.nodes.push_back(Node<float>());
    cT.nodes.back().makeCornerNode(corner, nodeNumber);
    return cT.nodes.size()-1;
}

// BUG: The node needs to be entered in the edgepoint arrays
template <int dim, class ctype>
NodeIdx PSurface<dim,ctype>::addIntersectionNodePair(int tri1, int tri2,
                                                const StaticVector<float,2>& dP1, const StaticVector<float,2>& dP2, 
                                                int edge1, int edge2, const StaticVector<float,3>& range)
{
    DomainTriangle& cT1 = triangles(tri1);
    DomainTriangle& cT2 = triangles(tri2);

    iPos.push_back(range);
    int nodeNumber = iPos.size()-1;

    cT1.nodes.push_back(Node<float>());
    int newNode1 = cT1.nodes.size()-1;
    cT2.nodes.push_back(Node<float>());
    
    cT1.nodes.back().setValue(dP1, nodeNumber, Node<float>::INTERSECTION_NODE);
    cT2.nodes.back().setValue(dP2, nodeNumber, Node<float>::INTERSECTION_NODE);

    cT1.nodes.back().setDomainEdge(edge1);
    cT2.nodes.back().setDomainEdge(edge2);
    return newNode1;
}

// BUG: The node needs to be entered in the edgepoint arrays
template <int dim, class ctype>
NodeIdx PSurface<dim,ctype>::addTouchingNode(int tri, const StaticVector<float,2>& dP, int edge, int nodeNumber)
{
    DomainTriangle& cT = triangles(tri);

    cT.nodes.push_back(Node<float>());
    
    cT.nodes.back().setValue(dP, nodeNumber, Node<float>::TOUCHING_NODE);
    cT.nodes.back().setDomainEdge(edge);
    return cT.nodes.size()-1;
}

// BUG: The node needs to be entered in the edgepoint arrays
template <int dim, class ctype>
NodeIdx PSurface<dim,ctype>::addTouchingNodePair(int tri1, int tri2,
                                            const StaticVector<float,2>& dP1, const StaticVector<float,2>& dP2, 
                                            int edge1, int edge2, int nodeNumber)
{
    DomainTriangle& cT1 = triangles(tri1);
    DomainTriangle& cT2 = triangles(tri2);

    cT1.nodes.push_back(Node<float>());
    cT2.nodes.push_back(Node<float>());
    
    cT1.nodes.back().setValue(dP1, nodeNumber, Node<float>::TOUCHING_NODE);
    cT2.nodes.back().setValue(dP2, nodeNumber, Node<float>::TOUCHING_NODE);

    cT1.nodes.back().setDomainEdge(edge1);
    cT2.nodes.back().setDomainEdge(edge2);

    return cT1.nodes.size()-1;
}

template <int dim, class ctype>
void PSurface<dim,ctype>::addParTriangle(int tri, const std::tr1::array<int,3>& p)
{
    DomainTriangle& cT = triangles(tri);

    assert(p[0]>=0 && p[0]<cT.nodes.size());
    assert(p[1]>=0 && p[1]<cT.nodes.size());
    assert(p[2]>=0 && p[2]<cT.nodes.size());

    if (!cT.nodes[p[0]].isConnectedTo(p[1]))
        cT.addEdge(p[0], p[1]);
    if (!cT.nodes[p[1]].isConnectedTo(p[2]))
        cT.addEdge(p[1], p[2]);
    if (!cT.nodes[p[2]].isConnectedTo(p[0]))
        cT.addEdge(p[2], p[0]);

}

template <int dim, class ctype>
NodeBundle PSurface<dim,ctype>::getNodeBundleAtVertex(int v) const
{
    NodeBundle result;
    std::vector<int> neighbors = getTrianglesPerVertex(v);

    result.resize(neighbors.size());

    for (int i=0; i<neighbors.size(); i++) {

        result[i].tri = neighbors[i];
        const DomainTriangle& cT = triangles(neighbors[i]);
        int corner = cT.getCorner(v);

        for (int j=0; j<cT.nodes.size(); j++)
            if ((cT.nodes[j].isCORNER_NODE() || cT.nodes[j].isGHOST_NODE()) &&
                cT.nodes[j].getCorner()==corner) {
                result[i].idx = j;
                break;
            }

    }

    return result;

}


// ////////////////////////////////////////////////////////
//   Explicit template instantiations.
//   If you need more, you can add them here.
// ////////////////////////////////////////////////////////

template class PSurface<2,float>;
template class PSurface<2,double>;

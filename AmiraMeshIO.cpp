#include <vector>

#ifdef PSURFACE_STANDALONE
#include "TargetSurface.h"
#else
#include "hxsurface/Surface.h"
#endif

#include "StaticVector.h"
#include "AmiraMeshIO.h"
#include "Domains.h"
#include "PSurface.h"
#include "PSurfaceFactory.h"

#if defined HAVE_AMIRAMESH || !defined PSURFACE_STANDALONE
#include <amiramesh/AmiraMesh.h>
#endif


/** \brief Tiny array class with controlled length

The AmiraMesh library cannot handle vectors of vectors.  Therefore the list
of triangle indices is returned as int* rather than array<int,3>*.  For easier
access I would like to cast to array<int,3>, but that is dangerous.  In various
gcc implementations of array, its size is more than sizeof(type)*N.
What I need is a replacement for array with guaranteed size, so the cast works.
*/

template <class T, int N>
struct MiniArray
{

    T& operator[](int index) {
        return data_[index];
    }

    const T& operator[](int index) const {
        return data_[index];
    }

    T data_[N];
};



#if defined HAVE_AMIRAMESH || !defined PSURFACE_STANDALONE
template <class ctype>
int AmiraMeshIO<ctype>::writeAmiraMesh(PSurface<2,ctype>* par, const char* filename)
{
    AmiraMesh am;

    int i, j;
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

    std::vector<StaticVector<ctype,3> > baseGridVertexCoordsArray(numVertices);
    for (size_t i=0; i<par->getNumVertices(); i++)
        baseGridVertexCoordsArray[i] = par->vertices(i);

    AmiraMesh::Data* vertexCoords = new AmiraMesh::Data("BaseGridVertexCoords", vertices, 
                                                        McPrimType::mc_float, 3, 
                                                        (void*)&baseGridVertexCoordsArray[0]);
    am.insert(vertexCoords);

    /////////////////////////////////////
    // store base grid triangles
    AmiraMesh::Location* triangles = new AmiraMesh::Location("BaseGridTriangles", numTriangles);
    am.insert(triangles);

    std::vector<MiniArray<int, 3> > baseGridTriArray(numTriangles);

    for (size_t i=0; i<par->getNumTriangles(); i++)
        for (int j=0; j<3; j++)
            baseGridTriArray[i][j] = par->triangles(i).vertices[j];
    

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

        const DomainTriangle<ctype>& cT = par->triangles(i);

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

    std::vector<StaticVector<ctype,2> > domainPositions(numNodes);
    std::vector<int>     nodeNumbers(numNodes);
    std::vector<MiniArray<int,2> > parameterEdgeArray(numParamEdges);
    std::vector<int>     edgePointsArray(numEdgePoints);

    int arrayIdx           = 0;
    int edgeArrayIdx       = 0;
    int edgePointsArrayIdx = 0;

    for (i=0; i<numTriangles; i++) {

        const DomainTriangle<ctype>& cT = par->triangles(i);

        std::vector<int> newIdx(cT.nodes.size());
        int localArrayIdx = 3;
        // the cornerNode are not saved, because everything about them
        // can be deduced from the base grid
        newIdx[cT.cornerNode(0)] = 0;
        newIdx[cT.cornerNode(1)] = 1;
        newIdx[cT.cornerNode(2)] = 2;

        // the three remaining types are saved separately, in order to avoid
        // explicitly saving the type for each node.
        for (size_t cN=0; cN<cT.nodes.size(); cN++) {
            if (cT.nodes[cN].isINTERSECTION_NODE()){

                domainPositions[arrayIdx] = cT.nodes[cN].domainPos();
                nodeNumbers[arrayIdx]     = cT.nodes[cN].getNodeNumber();

                newIdx[cN] = localArrayIdx;
                
                arrayIdx++;
                localArrayIdx++;
            }
        }

        for (size_t cN=0; cN<cT.nodes.size(); cN++) {

            if (cT.nodes[cN].isTOUCHING_NODE()){

                domainPositions[arrayIdx] = cT.nodes[cN].domainPos();
                nodeNumbers[arrayIdx]     = cT.nodes[cN].getNodeNumber();

                newIdx[cN] = localArrayIdx;

                arrayIdx++;
                localArrayIdx++;
            }
        }

        for (size_t cN=0; cN<cT.nodes.size(); cN++) {

            if (cT.nodes[cN].isINTERIOR_NODE()){

                domainPositions[arrayIdx] = cT.nodes[cN].domainPos();
                nodeNumbers[arrayIdx]     = cT.nodes[cN].getNodeNumber();

                newIdx[cN] = localArrayIdx;
                
                arrayIdx++;
                localArrayIdx++;
            }
        }
        
        // ///////////////////////////////////////
        // the parameterEdges for this triangle
        // ///////////////////////////////////////

        typename PlaneParam<ctype>::UndirectedEdgeIterator cE;
        for (cE = cT.firstUndirectedEdge(); cE.isValid(); ++cE){
            if (cE.isRegularEdge()) {
                parameterEdgeArray[edgeArrayIdx][0] = newIdx[cE.from()];
                parameterEdgeArray[edgeArrayIdx][1] = newIdx[cE.to()];
                edgeArrayIdx++;
            }
        }


        ///////////////////////////////////////
        // the edgePoints for this triangle
        ///////////////////////////////////////
        for (j=0; j<3; j++){

            for (size_t k=1; k<cT.edgePoints[j].size()-1; k++)
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
    //////////////////////////////
    if (!am.write(filename) ) {
        printf("An error has occured writing file %s.\n", filename);
        return 0;
    }

    return 1;

}

template <class ctype>
void* AmiraMeshIO<ctype>::readAmiraMesh(AmiraMesh* am, const char* filename)
{
    PSurface<2,ctype>* par = new PSurface<2,ctype>;

    Surface* surf = new Surface;

    if (!initFromAmiraMesh(par, am, filename, surf)) {
        delete par;
        delete surf;

        return NULL;
    }

    return par;
}

template <class ctype>
bool AmiraMeshIO<ctype>::initFromAmiraMesh(PSurface<2,ctype>* psurface, AmiraMesh* am, const char* filename, Surface* surf)
{
    // //////////////////////////////////
    //   Create PSurface factory
    // //////////////////////////////////

    PSurfaceFactory<2,ctype> factory(psurface);
    
    // Target surface already exists
    factory.setTargetSurface(surf);
    
    psurface->getPaths(am->parameters);

    ///////////////////////////////////////////////
    // test for file format
    ///////////////////////////////////////////////
    AmiraMesh::Data* AMnodePos = am->findData("NodePositions", HxFLOAT, 3, "NodePositions");

    if (!AMnodePos){
        printf("You're trying to read a file in the old-style AmiraMesh format.\n");
        printf("But the format is not supported anymore :-(\n");
        return false;
    }

    //////////////////////////////////////////////
    // the innerRegions and outerRegions arrays
    //////////////////////////////////////////////
    AmiraMesh::Data* AMpatches = am->findData("Patches", HxINT32, 3, "Patches");
    if (!AMpatches){
        printf("AmiraMesh: Field 'Patches' not found!\n");
        return false;
    }
    
    psurface->patches.resize(AMpatches->location()->dims()[0]);

    for (size_t i=0; i<psurface->patches.size(); i++)
        psurface->patches[i] = ((typename PSurface<2,ctype>::Patch*)AMpatches->dataPtr())[i];
    
    //////////////////////////////////
    // load the base grid vertices
    //////////////////////////////////
    AmiraMesh::Data* AMvertices = am->findData("BaseGridVertexCoords", HxFLOAT, 3, "BaseGridVertexCoords");
    if (!AMvertices){
        printf("AmiraMesh: Field 'BaseGridVertexCoords' not found!\n");
        return false;
    }
    
    int numPoints = AMvertices->location()->dims()[0];

    // copy points
    StaticVector<ctype,3> newVertex;
    for (int i=0; i<numPoints; i++) {
        for (int j=0; j<3; j++)
            newVertex[j] = ((float*)AMvertices->dataPtr())[i*3+j];
        factory.insertVertex(newVertex);
    }

    // /////////////////////////////////////////////////////
    //  copy node positions
    // /////////////////////////////////////////////////////
    psurface->iPos.resize(AMnodePos->location()->dims()[0]);
    
    for (size_t i=0; i<psurface->iPos.size(); i++)
        for (int j=0; j<3; j++)
            psurface->iPos[i][j] = ((float(*)[3])AMnodePos->dataPtr())[i][j];

    ////////////////////////////////////////////////////////
    // copy triangles.  This takes care of the edges, too
    ////////////////////////////////////////////////////////
    AmiraMesh::Data* AMtriangles = am->findData("BaseGridTriangles", HxINT32, 3, "BaseGridTriangles");
    if (!AMtriangles){
        printf("AmiraMesh: Field 'BaseGridTriangles' not found!\n");
        return false;
    }
    
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
    const MiniArray<float,2>*   nodeData       = (MiniArray<float,2>*)nodes->dataPtr();
    const int* nodeNumbers          = (int*)AMnodeNumbers->dataPtr();
    const MiniArray<int,2>* edgeData  = (MiniArray<int,2>*)AMedges->dataPtr();
    const int*     edgePointData    = (int*)AMedgePoints->dataPtr();


    int edgeCounter=0, edgePointCounter=0;
    int nodeArrayIdx = 0;
    
    for (int i=0; i<numTriangles; i++){
        
        //int newTriIdx = psurface->createSpaceForTriangle(triIdx[i][0], triIdx[i][1], triIdx[i][2]);
        std::tr1::array<unsigned int, 3> triangleVertices = {((int*)AMtriangles->dataPtr())[3*i+0], 
                                                             ((int*)AMtriangles->dataPtr())[3*i+1], 
                                                             ((int*)AMtriangles->dataPtr())[3*i+2]};
        int newTriIdx = factory.insertSimplex(triangleVertices);

        psurface->triangles(newTriIdx).patch = numNodesAndEdgesData[11*i+4];

        ///////////////////////////////////////////////
        // get the parametrization on this triangle
        ///////////////////////////////////////////////

        int numIntersectionNodes = numNodesAndEdgesData[11*i+0];
        int numTouchingNodes     = numNodesAndEdgesData[11*i+1];
        int numInteriorNodes     = numNodesAndEdgesData[11*i+2];

        int numParamEdges        = numNodesAndEdgesData[11*i+3];

        ////////////////////////////////////
        // first the nodes
        ////////////////////////////////////

        psurface->triangles(newTriIdx).nodes.resize(numIntersectionNodes + numTouchingNodes + numInteriorNodes + 3);

        // three corner nodes
        StaticVector<ctype,2> domainPos(1, 0);
        int nodeNumber = numNodesAndEdgesData[11*i+8];
        psurface->triangles(newTriIdx).nodes[0].setValue(domainPos, nodeNumber, Node<ctype>::CORNER_NODE);
        psurface->triangles(newTriIdx).nodes[0].makeCornerNode(0, nodeNumber);

        domainPos = StaticVector<ctype,2>(0, 1);
        nodeNumber = numNodesAndEdgesData[11*i+9];
        psurface->triangles(newTriIdx).nodes[1].setValue(domainPos, nodeNumber, Node<ctype>::CORNER_NODE);
        psurface->triangles(newTriIdx).nodes[1].makeCornerNode(1, nodeNumber);

        domainPos = StaticVector<ctype,2>(0, 0);
        nodeNumber = numNodesAndEdgesData[11*i+10];
        psurface->triangles(newTriIdx).nodes[2].setValue(domainPos, nodeNumber, Node<ctype>::CORNER_NODE);
        psurface->triangles(newTriIdx).nodes[2].makeCornerNode(2, nodeNumber);

        int nodeCounter = 3;

        // the intersection nodes
        for (int j=0; j<numIntersectionNodes; j++, nodeCounter++, nodeArrayIdx++){
            
            // float --> double
            StaticVector<ctype,2> domainPos;
            domainPos[0] = nodeData[nodeArrayIdx][0];
            domainPos[1] = nodeData[nodeArrayIdx][1];
            
            int nodeNumber    = nodeNumbers[nodeArrayIdx];

            psurface->triangles(newTriIdx).nodes[nodeCounter].setValue(domainPos, nodeNumber, Node<ctype>::INTERSECTION_NODE);
        }

        // the touching nodes
        for (int j=0; j<numTouchingNodes; j++, nodeCounter++, nodeArrayIdx++){

            // float --> double
            StaticVector<ctype,2> domainPos;
            domainPos[0] = nodeData[nodeArrayIdx][0];
            domainPos[1] = nodeData[nodeArrayIdx][1];

            int nodeNumber    = nodeNumbers[nodeArrayIdx];
            
            psurface->triangles(newTriIdx).nodes[nodeCounter].setValue(domainPos, nodeNumber, Node<ctype>::TOUCHING_NODE);
        }

        // the interior nodes
        for (int j=0; j<numInteriorNodes; j++, nodeCounter++, nodeArrayIdx++){

            // float --> double
            StaticVector<ctype,2> domainPos;
            domainPos[0] = nodeData[nodeArrayIdx][0];
            domainPos[1] = nodeData[nodeArrayIdx][1];

            int nodeNumber    = nodeNumbers[nodeArrayIdx];

            psurface->triangles(newTriIdx).nodes[nodeCounter].setValue(domainPos, nodeNumber, Node<ctype>::INTERIOR_NODE);
        }

        ///////////////////////////////
        // the parameterEdges
        ///////////////////////////////
        for (int j=0; j<numParamEdges; j++, edgeCounter++)
            psurface->triangles(newTriIdx).addEdge(edgeData[edgeCounter][0], edgeData[edgeCounter][1]);

        //////////////////////////////////
        // the edgePoints arrays
        //////////////////////////////////
        for (int j=0; j<3; j++){

            psurface->triangles(newTriIdx).edgePoints[j].resize(numNodesAndEdgesData[11*i+5+j] + 2);

            psurface->triangles(newTriIdx).edgePoints[j][0]     = j;
            psurface->triangles(newTriIdx).edgePoints[j].back() = (j+1)%3;

            for (int k=0; k<numNodesAndEdgesData[11*i+5+j]; k++){
                psurface->triangles(newTriIdx).edgePoints[j][k+1] = edgePointData[edgePointCounter];
                edgePointCounter++;
            }

        }

        //psurface->integrateTriangle(newTriIdx);
    }

    // sad but true ...
    psurface->hasUpToDatePointLocationStructure = false;

    psurface->setupOriginalSurface();

    return true;
}


// ////////////////////////////////////////////////////////
//   Explicit template instantiations.
//   If you need more, you can add them here.
// ////////////////////////////////////////////////////////

template class PSURFACE_EXPORT AmiraMeshIO<float>;
template class PSURFACE_EXPORT AmiraMeshIO<double>;

#endif


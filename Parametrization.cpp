/////////////////////////////////////////////////////////////////
//
// $Id: Parametrization.cpp,v 1.2 2007/12/17 10:13:39 sander Exp $
// 
// $Log: Parametrization.cpp,v $
// Revision 1.2  2007/12/17 10:13:39  sander
// adapted include path names from old package name 'parametrization' to new one 'psurface'
//
// Revision 1.1  2007/10/17 13:16:55  sander
// moved here from the ZIB server
//
// Revision 1.23  2006/10/20 10:18:54  bzfsande
// stupid bug fixed: the method std::vector<bool>::clear() does not do the same as McBitfield::unsetAll()
// mailtoauthor: sander@zib.de
//
// Revision 1.22  2006/03/14 16:29:43  bzfsande
// use std::vector<bool> instead of McBitfield in order to be more independent from Amira
// mailtoauthor: sander@zib.de
//
// Revision 1.21  2005/04/26 15:33:51  bzfsande
// include McDArray.h
// mailtoauthor: sander@zib.de
//
// Revision 1.20  2004/10/01 13:10:53  bzfsande
// object duplication implemented, assert replaced by exception-throwing
// mailtoauthor: sander@zib.de
//
// Revision 1.19  2003/10/09 12:07:37  bzfsande
// different fixes
//
// Revision 1.18  2003/08/27 12:22:19  bzfsande
// fixes
//
// Revision 1.17  2003/07/25 08:07:49  bzfsande
// old file reader deleted
//
// Revision 1.16  2003/06/30 12:44:44  bzfsande
// generalizations needed for the contact library
//
// Revision 1.15  2003/06/05 13:01:33  bzfsande
// introduces ghost nodes as a fifth node type
// also, the access to the node domain positions is procedural now
//
// Revision 1.14  2003/05/09 08:58:20  bzfsande
// bugfixes
//
// Revision 1.13  2003/04/04 14:59:18  bzfsande
// The base grid is now array-based
//
// Revision 1.12  2003/03/24 13:26:53  bzfsande
// separated the Parametrization object into a base grid object 
// connected to a standard Surface
//
// Revision 1.11  2002/10/08 13:00:29  bzfsande
// mapping function with an explicit seed
//
// Revision 1.10  2002/10/02 15:22:38  bzfsande
// Introduced the global #imagePos# array
// - any data defined on the vertices of the original surface
//   can now be queried
// - The base grid vertices don't have to be at the same
//   position as their homologues on the original surface anymore
//
// Revision 1.9  2001/12/02 18:13:24  bzfsande
// fixes
//
// Revision 1.8  2001/11/02 15:28:27  bzflamec
// compile on windows
//
// Revision 1.7  2001/10/16 14:30:31  bzfsande
// fixes galore
//
// Revision 1.6  2001/10/15 09:31:07  bzfsande
// fixes and cosmetics
//
// Revision 1.5  2001/09/27 08:53:59  bzfsande
// fixes
//
// Revision 1.4  2001/09/26 09:33:00  bzfsande
// edgeIterators added
//
// Revision 1.3  2001/09/25 13:42:47  bzfsande
// the plane graph data structure is now index-based
//
// Revision 1.2  2001/09/21 08:06:34  bzfsande
// .
//
// Revision 1.1  2001/09/13 12:39:52  bzfsande
// initial version
//
// Revision 1.16  2001/09/12 09:11:59  bzfsande
// .
//
// Revision 1.12  2001/08/17 15:37:58  bzfsande
// more fixes
//
// Revision 1.11  2001/08/13 15:34:47  bzfsande
// more fixes
//
// Revision 1.10  2001/08/08 14:41:11  bzfsande
// old bugs gone - new bugs added
//
// Revision 1.9  2001/06/25 08:26:52  bzfsande
// bugfixes
//
// Revision 1.8  2001/06/15 10:13:55  bzfsande
// New changes in the surface parametrization stuff:
// - The HxParametrizer modules does only the conversion from
//   HxSurface to HxParametrization now ...
// - ... plus it can now parametrize an HxSurface over an
//   explicitly given base grid!
// - The simplification stuff has been moved to an editor
// - The different smoothing algorithms have been moved to
//   another editor
// Enjoy!
//
// Revision 1.7  2001/05/09 15:30:27  bzfsande
// Floater's Parametrization now usable
//
// Revision 1.6  2001/04/27 13:04:37  bzfsande
// new stuff:
// - reader & writer for AmiraMesh
// - A histogram plotting function for evaluation of
//   the parametrization quality
// - first attempts on constructing a GLOBAL parametrization
//
// Revision 1.5  2001/04/06 12:24:49  bzfsande
// a test for self-intersections and topology changes, simplification using a 
// Hausdorff-type distance function and a brand-new internal structure
//
// Revision 1.4  2001/02/16 11:53:57  bzfsande
// added:
// - the Brown/Faigle point location algorithm
// - guarantee for foldover-free parametrizations
// - barycentric parametrization
// - tons of bugfixes
//
// Revision 1.3  2001/01/31 17:24:33  bzflamec
// compile on IRIX
//
// Revision 1.2  2001/01/31 16:20:29  bzfsande
// a complete, but buggy, version of MAPS
//
// Revision 1.1  2000/12/05 13:38:47  bzfsande
// MAPS
//
//
/////////////////////////////////////////////////////////////////

#include <tr1/array>
#include <vector>

#include <amiramesh/AmiraMesh.h>

#include <psurface/Parametrization.h>
#include <psurface/GlobalNodeIdx.h>

Parametrization::Parametrization(HxParamBundle* bundle)
{
    if (bundle) {
        params = bundle;
        hasOwnParamBundle = false;
    } else {
        params = new HxParamBundle;
        hasOwnParamBundle = true;
    }
}

Parametrization::~Parametrization() 
{ 
    if (hasOwnParamBundle)
        delete params;
}

void Parametrization::clear()
{
    surface = NULL;
    patches.clear();
    if (hasOwnParamBundle)
        delete params;

    iPos.clear();
    paths.clear();
    McSurfaceBase<DomainVertex, DomainEdge, DomainTriangle>::clear();
}

void Parametrization::getBoundingBox(Box<float,3>& bbox) const
{
    if (getNumVertices()==0)
        return;

    bbox.set(vertices(0), vertices(0));

    for (int i=1; i<getNumVertices(); i++)
        bbox.extendBy(vertices(i));
}

void Parametrization::init(const Parametrization* other)
{
    // copy domain surface.  
    *this = *other;

    surface = new Surface();
    *surface = *other->surface;
    
}

McVec2f Parametrization::getLocalTargetCoords(const GlobalNodeIdx& n, int targetTri) const
{
    const Node& cN = triangles(n.tri).nodes[n.idx];
    
    switch (cN.type) {
    case Node::GHOST_NODE:
    case Node::INTERSECTION_NODE: {

        StaticVector<float,3> iPos = imagePos(n.tri, n.idx);
        return triangles(n.tri).computeBarycentricCoords(iPos, surface->points[surface->triangles[targetTri].points[0]], 
                                            surface->points[surface->triangles[targetTri].points[1]], 
                                            surface->points[surface->triangles[targetTri].points[2]]);
    }
    default:
        if (cN.getNodeNumber()==surface->triangles[targetTri].points[0])
            return McVec2f(1, 0);
        else if (cN.getNodeNumber()==surface->triangles[targetTri].points[1])
            return McVec2f(0, 1);
        else if (cN.getNodeNumber()==surface->triangles[targetTri].points[2])
            return McVec2f(0, 0);
        else {
            printf("The node is not related to the targetTri!\n");
            throw ParamError();
        }
    }
}


GlobalNodeIdx Parametrization::getOtherEndNode(int triIdx, NodeIdx cN) const
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

int Parametrization::getNumNodes() const
{
    int n = 0;
    for (int i=0; i<getNumTriangles(); i++)
        n += triangles(i).nodes.size();
    return n;
}

int Parametrization::getNumTrueNodes() 
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

void Parametrization::removeExtraEdges()
{

    for (int i(0); i<getNumTriangles(); i++) {
        triangles(i).removeExtraEdges();
    }

    hasUpToDatePointLocationStructure = false;
}

#if 0
void Parametrization::insertExtraEdges()
{
    DomainTriangle *cT = triangles.first();
    while (cT){
        cT->insertExtraEdges();
        cT = triangles.succ(cT);
    }
}
#endif

void Parametrization::createPointLocationStructure()
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
void Parametrization::garbageCollection()
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

    McSurfaceBase<DomainVertex, DomainEdge, DomainTriangle>::garbageCollection();
    
#ifndef NDEBUG
        printf("   ...Garbage collection finished!\n");
#endif
    }
    

int Parametrization::writeAmiraMesh(Parametrization* par, const char* filename)
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

    std::vector<McSArray<int, 3> > baseGridTriArray(numTriangles);

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

    std::vector<McVec2f> domainPositions(numNodes);
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

        PlaneParam::UndirectedEdgeIterator cE;
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

Parametrization* Parametrization::readAmiraMesh(AmiraMesh* am, const char* filename)
{
    Parametrization* par = new Parametrization;

    Surface* surf = new Surface;

    if (!par->initFromAmiraMesh(am, filename, surf)) {
        delete par;
        delete surf;

        return NULL;
    }

    *par->params = am->parameters;

    return par;
}


bool Parametrization::initFromAmiraMesh(AmiraMesh* am, const char* filename, Surface* surf)
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
        patches[i] = ((Parametrization::Patch*)AMpatches->dataPtr())[i];
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
    
    const McVec3i* triIdx = (McVec3i*)AMtriangles->dataPtr();
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
    const McVec2f*   nodeData       = (McVec2f*)nodes->dataPtr();
    const int* nodeNumbers          = (int*)AMnodeNumbers->dataPtr();
    const std::tr1::array<int,2>* edgeData  = (std::tr1::array<int,2>*)AMedges->dataPtr();
    const int*     edgePointData    = (int*)AMedgePoints->dataPtr();


    int edgeCounter=0, edgePointCounter=0;
    int nodeArrayIdx = 0;
    
    for (i=0; i<numTriangles; i++){
        
        int newTriIdx = createSpaceForTriangle(triIdx[i].i, triIdx[i].j, triIdx[i].k);

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
        McVec2f domainPos(1, 0);
        int nodeNumber = numNodesAndEdgesData[11*i+8];
        triangles(newTriIdx).nodes[0].setValue(domainPos, nodeNumber, Node::CORNER_NODE);
        cornerNodes[0] = 0;

        domainPos = McVec2f(0, 1);
        nodeNumber = numNodesAndEdgesData[11*i+9];
        triangles(newTriIdx).nodes[1].setValue(domainPos, nodeNumber, Node::CORNER_NODE);
        cornerNodes[1] = 1;

        domainPos = McVec2f(0, 0);
        nodeNumber = numNodesAndEdgesData[11*i+10];
        triangles(newTriIdx).nodes[2].setValue(domainPos, nodeNumber, Node::CORNER_NODE);
        cornerNodes[2] = 2;

        int nodeCounter = 3;

        // the intersection nodes
        for (j=0; j<numIntersectionNodes; j++, nodeCounter++, nodeArrayIdx++){
            McVec2f domainPos = nodeData[nodeArrayIdx];
            int nodeNumber    = nodeNumbers[nodeArrayIdx];

            triangles(newTriIdx).nodes[nodeCounter].setValue(domainPos, nodeNumber, Node::INTERSECTION_NODE);
        }

        // the touching nodes
        for (j=0; j<numTouchingNodes; j++, nodeCounter++, nodeArrayIdx++){
            McVec2f domainPos = nodeData[nodeArrayIdx];
            int nodeNumber    = nodeNumbers[nodeArrayIdx];
            
            triangles(newTriIdx).nodes[nodeCounter].setValue(domainPos, nodeNumber, Node::TOUCHING_NODE);
        }

        // the interior nodes
        for (j=0; j<numInteriorNodes; j++, nodeCounter++, nodeArrayIdx++){
            McVec2f domainPos = nodeData[nodeArrayIdx];
            int nodeNumber    = nodeNumbers[nodeArrayIdx];

            triangles(newTriIdx).nodes[nodeCounter].setValue(domainPos, nodeNumber, Node::INTERIOR_NODE);
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
            triangles(newTriIdx).edgePoints[j].last() = (j+1)%3;

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

void Parametrization::setupOriginalSurface()
{
    int i, j, k;
    
    if (!hasUpToDatePointLocationStructure)
        createPointLocationStructure();
    
    // transfer materials

    

    // create patches

    surface->patches.resize(numPatches());

    for (i=0; i<surface->patches.size(); i++){
        surface->patches[i] = new Surface::Patch();
        surface->patches[i]->innerRegion = patches[i].innerRegion;
        surface->patches[i]->outerRegion = patches[i].outerRegion;
        surface->patches[i]->boundaryId  = patches[i].boundaryId;
    }
    
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

            Node& cN = cT.nodes[i];
            McVec3i v;

            v[0] = cN.nodeNumber;

            switch (cN.type) {

            case Node::INTERSECTION_NODE:
                continue;

            case Node::INTERIOR_NODE:
                
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

            case Node::TOUCHING_NODE:
            case Node::CORNER_NODE:
                

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

void Parametrization::appendTriangleToOriginalSurface(const McVec3i& v, int patch)
{
    surface->triangles.appendSpace(1);
                        
    surface->triangles.last().points[0] = v[0];
    surface->triangles.last().points[1] = v[1];
    surface->triangles.last().points[2] = v[2];
    
    surface->triangles.last().patch = patch;
    surface->patches[patch]->triangles.append(surface->triangles.size()-1);
}


void Parametrization::getPaths(const HxParamBundle& parameters)
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


void Parametrization::savePaths(HxParamBundle& parameters)
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



int Parametrization::map(int triIdx, McVec2f& p, McVec3i& vertices, 
                         McVec2f& coords, int seed) const
{
    int i;
    const DomainTriangle& tri = triangles(triIdx);
    const std::vector<StaticVector<float,3> >& nP = iPos;

    // this is boundary handling
    if (p.x < 0.001){
        for (i=0; i<tri.edgePoints[1].size()-1; i++){

            const McVec2f& a = tri.nodes[tri.edgePoints[1][i]].domainPos();
            const McVec2f& b = tri.nodes[tri.edgePoints[1][i+1]].domainPos();

            if (a.y+1e-5 > p.y && p.y > b.y-1e-5) {

                McSArray<GlobalNodeIdx, 3> targetNodes;
                handleMapOnEdge(triIdx, p, a, b, 1, i, targetNodes, coords);
                vertices[0] = nodes(targetNodes[0]).getNodeNumber();
                vertices[1] = nodes(targetNodes[1]).getNodeNumber();
                vertices[2] = nodes(targetNodes[2]).getNodeNumber();
                return true;
            }
            
        }

    }else if (p.y < 0.001){

        for (i=0; i<tri.edgePoints[2].size()-1; i++){
            const McVec2f& a = tri.nodes[tri.edgePoints[2][i]].domainPos();
            const McVec2f& b = tri.nodes[tri.edgePoints[2][i+1]].domainPos();

            if (b.x+1e-5 > p.x && p.x > a.x-1e-5) {

                McSArray<GlobalNodeIdx, 3> targetNodes;
                handleMapOnEdge(triIdx, p, a, b, 2, i, targetNodes, coords);
                vertices[0] = nodes(targetNodes[0]).getNodeNumber();
                vertices[1] = nodes(targetNodes[1]).getNodeNumber();
                vertices[2] = nodes(targetNodes[2]).getNodeNumber();
                return true;
            }
        }

    }else if (p.x+p.y > 0.999) {

        for (i=0; i<tri.edgePoints[0].size()-1; i++){
            const McVec2f& a = tri.nodes[tri.edgePoints[0][i]].domainPos();
            const McVec2f& b = tri.nodes[tri.edgePoints[0][i+1]].domainPos();

            if (a.x+1e-5>p.x && p.x>b.x-1e-5) {

                McSArray<GlobalNodeIdx, 3> targetNodes;
                handleMapOnEdge(triIdx, p, a, b, 0, i, targetNodes, coords);
                vertices[0] = nodes(targetNodes[0]).getNodeNumber();
                vertices[1] = nodes(targetNodes[1]).getNodeNumber();
                vertices[2] = nodes(targetNodes[2]).getNodeNumber();
                return true;

            }
        }

        return false;

    }

    McSArray<NodeIdx, 3> v;
    int status = tri.map(p, v, coords, seed);

    if (!status)
        return 0;

    StaticVector<float,3> imagePos = PlaneParam::linearInterpol<StaticVector<float,3> >(coords, nP[tri.nodes[v[0]].getNodeNumber()],
                                                           nP[tri.nodes[v[1]].getNodeNumber()],
                                                           nP[tri.nodes[v[2]].getNodeNumber()]);

    // ///////////////////////////////////////////////////////
    // make sure we don't return intersection nodes
    McSArray<GlobalNodeIdx, 3> resultNodes;
    getActualVertices(triIdx, v, resultNodes);
    vertices[0] = nodes(resultNodes[0]).getNodeNumber();
    vertices[1] = nodes(resultNodes[1]).getNodeNumber();
    vertices[2] = nodes(resultNodes[2]).getNodeNumber();

    coords = tri.computeBarycentricCoords(imagePos, nP[vertices[0]], nP[vertices[1]], nP[vertices[2]]);

    return true;
}

void Parametrization::getActualVertices(int tri, const McSArray<NodeIdx, 3>& nds,
                                        McSArray<GlobalNodeIdx, 3>& vertices) const
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

int Parametrization::getImageSurfaceTriangle(int tri,
                                             const McSArray<NodeIdx, 3>& nds
                                             ) const
{
    int i;
    
    McSArray<GlobalNodeIdx, 3> actualVertices;
    McSArray<McSmallArray<int, 6>, 3> trianglesPerNode;

    getActualVertices(tri, nds, actualVertices);
    
    assert(surface->trianglesPerPoint.size());
    
    for (i=0; i<3; i++)
        trianglesPerNode[i] = getTargetTrianglesPerNode(actualVertices[i]);
    
    for (i=0; i<trianglesPerNode[0].size(); i++) {

        if (mcSmallArray::index(trianglesPerNode[1], trianglesPerNode[0][i])!=-1 &&
            mcSmallArray::index(trianglesPerNode[2], trianglesPerNode[0][i])!=-1)
            return trianglesPerNode[0][i];

    }
            
    return -1;
}

McSmallArray<int, 6> Parametrization::getTargetTrianglesPerNode(const GlobalNodeIdx& n) const
{
    assert(surface->trianglesPerPoint.size());
    const Node& cN = triangles(n.tri).nodes[n.idx];
    const float eps = 1e-6;

    switch (cN.type) {
    case Node::GHOST_NODE: {
//         printf("targetTri: %d (%d %d %d)      targetLocalCoords: (%f %f)\n", cN.getNodeNumber(),
//                surface->triangles[cN.getNodeNumber()].points[0],
//                surface->triangles[cN.getNodeNumber()].points[1],
//                surface->triangles[cN.getNodeNumber()].points[2],
//                cN.dP.x, cN.dP.y);
        
        McSmallArray<int, 6> result(1);
        result[0] = cN.getNodeNumber();
        if (cN.dP.x + cN.dP.y > 1-eps) {
            // append the triangles bordering on edge 0
            int p = surface->triangles[result[0]].points[0];
            int q = surface->triangles[result[0]].points[1];
            //printf("1) using:  p -> %d   q-> %d \n", p, q);
            getTrianglesPerEdge(p, q, result, result[0]);

        } else if (cN.dP.x < eps) {
            // append the triangles bordering on edge 1
            int p = surface->triangles[result[0]].points[1];
            int q = surface->triangles[result[0]].points[2];
            //printf("2) using:  p -> %d   q-> %d \n", p, q);
            getTrianglesPerEdge(p, q, result, result[0]);

        } else if (cN.dP.y < eps) {
            // append the triangles bordering on edge 2
            int p = surface->triangles[result[0]].points[2];
            int q = surface->triangles[result[0]].points[0];
            //printf("3) using:  p -> %d   q-> %d \n", p, q);
            getTrianglesPerEdge(p, q, result, result[0]);
            
        }
        
        return result;
        
        
    }
    case Node::INTERSECTION_NODE:
        assert(false);
    }

    return surface->trianglesPerPoint[cN.getNodeNumber()];
    
}

/// This is a service routine only for getTargetTrianglesPerNode
void Parametrization::getTrianglesPerEdge(int from, int to, McSmallArray<int, 6>& tris, int exception) const
{
    for (int i=0; i<surface->trianglesPerPoint[from].size(); i++) {

        if (mcSmallArray::index(surface->trianglesPerPoint[to], surface->trianglesPerPoint[from][i]) != -1 &&
            surface->trianglesPerPoint[from][i] != exception)
            tris.append(surface->trianglesPerPoint[from][i]);

    }

}

void Parametrization::handleMapOnEdge(int triIdx, const McVec2f& p, const McVec2f& a, const McVec2f& b,
                                      int edge, int edgePos, McSArray<GlobalNodeIdx, 3>& vertices, McVec2f& coords) const
{
    const DomainTriangle& tri = triangles(triIdx);
    float lambda = (p-a).length() / (a-b).length();

    StaticVector<float,3> targetPos = PlaneParam::linearInterpol<StaticVector<float,3> >(lambda, 
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
        
        PlaneParam::DirectedEdgeIterator edge = tri.getDirectedEdgeIterator(n1, n2);
        vertices[2] = getOtherEndNode(triIdx, edge.getONext().to());

    }

    assert(nodes(vertices[0]).getNodeNumber() != nodes(vertices[1]).getNodeNumber() && 
           nodes(vertices[1]).getNodeNumber() != nodes(vertices[2]).getNodeNumber() && 
           nodes(vertices[2]).getNodeNumber() != nodes(vertices[0]).getNodeNumber());
    
    coords = tri.computeBarycentricCoords(targetPos, imagePos(vertices[0]), 
                                          imagePos(vertices[1]), 
                                          imagePos(vertices[2]));

}


int Parametrization::positionMap(int triIdx, McVec2f& p, StaticVector<float,3>& result) const
{
    McVec2f localCoords;
    McVec3i tri;

    int status = map(triIdx, p, tri, localCoords);
    
    if (!status) {
        printf("p: (%f %f)\n", p.x, p.y);
        triangles(triIdx).print(true, true, false);
        assert(false);
        return false;
    }

    result = PlaneParam::linearInterpol<StaticVector<float,3> >(localCoords, iPos[tri[0]], iPos[tri[1]], iPos[tri[2]]);

    return true;
}



int Parametrization::directNormalMap(int triIdx, McVec2f& p, StaticVector<float,3>& result) const
{
    McVec2f localCoords;
    McVec3i tri;

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

int Parametrization::invertTriangles(int patch)
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

// NodeIdx Parametrization::addNode(int tri, const StaticVector<float,3>& p)
// {
//     DomainTriangle& cT = pars[side].triangles(tri);

//     pars[side].iPos.append(p);

//     int nodeNumber = pars[side].iPos.size()-1;

//     return cT.nodes.append(Node(McVec2f(0,0), nodeNumber, Node::INTERIOR_NODE));  

// }

NodeIdx Parametrization::addInteriorNode(int tri, const McVec2f& dom, int nodeNumber)
{
    triangles(tri).nodes.push_back(Node(dom, nodeNumber, Node::INTERIOR_NODE));
    return triangles(tri).nodes.size()-1;
}

NodeIdx Parametrization::addGhostNode(int tri, int corner, int targetTri, const McVec2f& localTargetCoords)
{
    triangles(tri).nodes.push_back(Node());
    triangles(tri).nodes.back().makeGhostNode(corner, targetTri, localTargetCoords);
    return triangles(tri).nodes.size()-1;
}

NodeIdx Parametrization::addCornerNode(int tri, int corner, int nodeNumber)
{
    DomainTriangle& cT = triangles(tri);

    cT.nodes.push_back(Node());
    cT.nodes.back().makeCornerNode(corner, nodeNumber);
    return cT.nodes.size()-1;
}

// BUG: The node needs to be entered in the edgepoint arrays
NodeIdx Parametrization::addIntersectionNodePair(int tri1, int tri2,
                                                const McVec2f& dP1, const McVec2f& dP2, 
                                                int edge1, int edge2, const StaticVector<float,3>& range)
{
    DomainTriangle& cT1 = triangles(tri1);
    DomainTriangle& cT2 = triangles(tri2);

    iPos.push_back(range);
    int nodeNumber = iPos.size()-1;

    cT1.nodes.push_back(Node());
    int newNode1 = cT1.nodes.size()-1;
    cT2.nodes.push_back(Node());
    
    cT1.nodes.back().setValue(dP1, nodeNumber, Node::INTERSECTION_NODE);
    cT2.nodes.back().setValue(dP2, nodeNumber, Node::INTERSECTION_NODE);

    cT1.nodes.back().setDomainEdge(edge1);
    cT2.nodes.back().setDomainEdge(edge2);
    return newNode1;
}

// BUG: The node needs to be entered in the edgepoint arrays
NodeIdx Parametrization::addTouchingNode(int tri, const McVec2f& dP, int edge, int nodeNumber)
{
    DomainTriangle& cT = triangles(tri);

    cT.nodes.push_back(Node());
    
    cT.nodes.back().setValue(dP, nodeNumber, Node::TOUCHING_NODE);
    cT.nodes.back().setDomainEdge(edge);
    return cT.nodes.size()-1;
}

// BUG: The node needs to be entered in the edgepoint arrays
NodeIdx Parametrization::addTouchingNodePair(int tri1, int tri2,
                                            const McVec2f& dP1, const McVec2f& dP2, 
                                            int edge1, int edge2, int nodeNumber)
{
    DomainTriangle& cT1 = triangles(tri1);
    DomainTriangle& cT2 = triangles(tri2);

    cT1.nodes.push_back(Node());
    cT2.nodes.push_back(Node());
    
    cT1.nodes.back().setValue(dP1, nodeNumber, Node::TOUCHING_NODE);
    cT2.nodes.back().setValue(dP2, nodeNumber, Node::TOUCHING_NODE);

    cT1.nodes.back().setDomainEdge(edge1);
    cT2.nodes.back().setDomainEdge(edge2);

    return cT1.nodes.size()-1;
}

void Parametrization::addParTriangle(int tri, const McVec3i& p)
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

NodeBundle Parametrization::getNodeBundleAtVertex(int v) const
{
    NodeBundle result;
    McSmallArray<int, 12> neighbors = getTrianglesPerVertex(v);

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

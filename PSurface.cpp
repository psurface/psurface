#include <vector>

// Check for VC9 / VS2008 with installed feature pack.
#if defined(_MSC_VER) && (_MSC_VER>=1500)
    #if defined(_CPPLIB_VER) && _CPPLIB_VER>=505
        #include <array>
    #else
        #error Please install the Visual Studio 2008 SP1 for TR1 support.
    #endif
#else
    #include <tr1/array>
#endif

#ifdef _MSC_VER
    // Required to make cmath define M_PI etc.
    #define _USE_MATH_DEFINES
#endif
#include <cmath>

#ifdef PSURFACE_STANDALONE
#include "TargetSurface.h"
#else
#include "hxsurface/Surface.h"
#endif

#if defined HAVE_AMIRAMESH || !defined PSURFACE_STANDALONE
#include <amiramesh/AmiraMesh.h>
#endif

#include "PSurface.h"
#include "GlobalNodeIdx.h"
#include "Box.h"

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

template <int dim, class ctype>
PSurface<dim,ctype>::~PSurface() 
{}

template <int dim, class ctype>
void PSurface<dim,ctype>::clear()
{
    surface = NULL;
    patches.clear();

    iPos.clear();
    paths.clear();
    SurfaceBase<Vertex<ctype>, Edge, DomainTriangle<ctype> >::clear();
}

template <int dim, class ctype>
void PSurface<dim,ctype>::getBoundingBox(Box<ctype,3>& bbox) const
{
    if (this->getNumVertices()==0)
        return;

    bbox.set(this->vertices(0), this->vertices(0));

    for (int i=1; i<this->getNumVertices(); i++)
        bbox.extendBy(this->vertices(i));
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
StaticVector<ctype,2> PSurface<dim,ctype>::getLocalTargetCoords(const GlobalNodeIdx& n, int targetTri) const
{
    const Node<ctype>& cN = this->triangles(n.tri).nodes[n.idx];
    
    switch (cN.type) {
    case Node<ctype>::GHOST_NODE:
    case Node<ctype>::INTERSECTION_NODE: {

        StaticVector<ctype,3> iPos = imagePos(n.tri, n.idx);

        // Convert from McVec3f to StaticVector
        std::tr1::array<StaticVector<ctype,3>, 3> p;

        for (int i=0; i<3; i++)
            for (int j=0; j<3; j++)
                p[i][j] = surface->points[surface->triangles[targetTri].points[i]][j];

        return this->triangles(n.tri).computeBarycentricCoords(iPos, p[0], p[1], p[2]);

    }
    default:
        if (cN.getNodeNumber()==surface->triangles[targetTri].points[0])
            return StaticVector<ctype,2>(1, 0);
        else if (cN.getNodeNumber()==surface->triangles[targetTri].points[1])
            return StaticVector<ctype,2>(0, 1);
        else if (cN.getNodeNumber()==surface->triangles[targetTri].points[2])
            return StaticVector<ctype,2>(0, 0);
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
    
    while (this->triangles(triIdx).nodes[cN].isINTERSECTION_NODE()) {
        
        const DomainTriangle<ctype>& cT = this->triangles(triIdx);
        
        // get the edge the node is on and its position in the edgePoints array
        int edge  = cT.nodes[cN].getDomainEdge();
        int edgePos = cT.nodes[cN].getDomainEdgePosition();
        
        // get adjacent triangle
        const int cE = cT.getOppositeEdge(cT.vertices[(edge+2)%3]);

#ifndef NDEBUG
        if (this->edges(cE).numTriangles()!=2) {
            printf("Edge:  %d --> %d\n", this->edges(cE).from, this->edges(cE).to);
            for (i=0; i<this->edges(cE).numTriangles(); i++)
                this->triangles(this->edges(cE).triangles[i]).print(true, true, true);
        }
#endif

        assert(this->edges(cE).numTriangles()==2);
        
        const int oppT = (this->edges(cE).triangles[0]==triIdx) 
            ? this->edges(cE).triangles[1] 
            : this->edges(cE).triangles[0];
        
        // get the opposite edgePoint array
        int oppEdge = -1;
        bool reverse = false;
        for (i=0; i<3; i++) {
            if (this->triangles(oppT).vertices[i] == cT.vertices[edge] && 
                this->triangles(oppT).vertices[(i+1)%3]==cT.vertices[(edge+1)%3]) {
                
                oppEdge = i;
                reverse = false;
                break;
            } else if (this->triangles(oppT).vertices[i] == cT.vertices[(edge+1)%3] && 
                       this->triangles(oppT).vertices[(i+1)%3] == cT.vertices[edge]) {
                
                oppEdge = i;
                reverse = true;
                break;
            }
        }

        assert(oppEdge!=-1);
        
        int oppEdgePos = (reverse) ? cT.edgePoints[edge].size()-edgePos-1 : edgePos;

        if (this->triangles(oppT).nodes[this->triangles(oppT).edgePoints[oppEdge][oppEdgePos]].getNodeNumber()
            != cT.nodes[cT.edgePoints[edge][edgePos]].getNodeNumber()) {

            printf("Condition triangles(oppT).nodes[triangles(oppT).edgePoints[oppEdge][oppEdgePos]].getNodeNumber() != cT.nodes[cT.edgePoints[edge][edgePos]].getNodeNumber()  failed!\n");
            
            throw ParamError();

        }
        
        int newPoint = this->triangles(oppT).nodes[this->triangles(oppT).edgePoints[oppEdge][oppEdgePos]].theInteriorNode();
        
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
    for (int i=0; i<this->getNumTriangles(); i++)
        n += this->triangles(i).nodes.size();
    return n;
}

template <int dim, class ctype>
int PSurface<dim,ctype>::getNumTrueNodes() 
{
    int highestTrueNodeNumber = -1;

    for (int j(0); j<this->getNumTriangles(); j++) {

        const DomainTriangle<ctype>& cT = this->triangles(j);

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

    for (int i(0); i<this->getNumTriangles(); i++) {
        this->triangles(i).removeExtraEdges();
    }

    hasUpToDatePointLocationStructure = false;
}

template <int dim, class ctype>
StaticVector<ctype,3> PSurface<dim,ctype>::imagePos(int tri, NodeIdx node) const 
{
    const Node<ctype>& cN = this->triangles(tri).nodes[node];
    
    switch (cN.type) {
    case Node<ctype>::GHOST_NODE: {
        const Surface::Triangle& cT = surface->triangles[cN.getNodeNumber()];
        
        StaticVector<ctype,3> p0(surface->points[cT.points[0]][0],surface->points[cT.points[0]][1],surface->points[cT.points[0]][2]);
        StaticVector<ctype,3> p1(surface->points[cT.points[1]][0],surface->points[cT.points[1]][1],surface->points[cT.points[1]][2]);
        StaticVector<ctype,3> p2(surface->points[cT.points[2]][0],surface->points[cT.points[2]][1],surface->points[cT.points[2]][2]);
        
        return PlaneParam<ctype>::linearInterpol(cN.dP, p0, p1, p2);
    }
    case Node<ctype>::INTERSECTION_NODE:
        return iPos[cN.getNodeNumber()];
        
    default:
        // Do a componentwise copy to get from McVec3f to StaticVector<ctype>
        StaticVector<ctype,3> result;
        return StaticVector<ctype,3>(surface->points[cN.getNodeNumber()][0],
                                     surface->points[cN.getNodeNumber()][1],
                                     surface->points[cN.getNodeNumber()][2]);
    }
    
}


template <int dim, class ctype>
void PSurface<dim,ctype>::createPointLocationStructure()
{
    for (int i(0); i<this->getNumTriangles(); i++){
        this->triangles(i).checkConsistency("Before Insert");
        this->triangles(i).insertExtraEdges();
        this->triangles(i).createPointLocationStructure();
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
    std::cout << this->freeVertexStack.size() << " vertices, "
              << this->freeEdgeStack.size()   << " edges, "
              << this->freeTriangleStack.size() << " triangles removed" << std::endl;
#endif
    
    std::vector<bool> isInvalid;
    
    // clean up vertices
    if (this->freeVertexStack.size()) {
        
        int offset = 0;
        
        std::vector<int> vertexOffsets(this->vertexArray.size());
        isInvalid.resize(this->vertexArray.size());
        for (i=0; i<isInvalid.size(); i++)
            isInvalid[i] = false;
        
        for (i=0; i<this->freeVertexStack.size(); i++)
            isInvalid[this->freeVertexStack[i]] = true;
        
        for (i=0; i<this->vertexArray.size(); i++){
            vertexOffsets[i] = offset;
            
            if (isInvalid[i]) 
                offset++;
        }
        
        ////////////////////
        for (i=0; i<vertexOffsets.size(); i++)
            this->vertexArray[i-vertexOffsets[i]] = this->vertexArray[i];
        
        this->vertexArray.resize(this->vertexArray.size()-offset);
        
        // Adjust edges
        for (i=0; i<this->edgeArray.size(); i++) {
            this->edgeArray[i].from -= vertexOffsets[this->edgeArray[i].from];
            this->edgeArray[i].to   -= vertexOffsets[this->edgeArray[i].to];
        }                
        
        // Adjust triangles
        for (i=0; i<this->triangleArray.size(); i++) 
            for (j=0; j<3; j++) 
                this->triangleArray[i].vertices[j] -= vertexOffsets[this->triangleArray[i].vertices[j]];
        
        // Adjust paths
        for (i=0; i<paths.size(); i++)
            for (j=0; j<paths[i].points.size(); j++)
                paths[i].points[j] -= vertexOffsets[paths[i].points[j]];

        this->freeVertexStack.clear();                

    }

    SurfaceBase<Vertex<ctype>, Edge, DomainTriangle<ctype> >::garbageCollection();
    
#ifndef NDEBUG
        printf("   ...Garbage collection finished!\n");
#endif
}
    


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

    for (k=0; k<this->getNumTriangles(); k++) {

        DomainTriangle<ctype>& cT = this->triangles(k);

        ////////////////////////////////
        int numNodes = cT.nodes.size();

        for (i=0; i<numNodes; i++) {

            Node<ctype>& cN = cT.nodes[i];
            std::tr1::array<int,3> v;

            v[0] = cN.nodeNumber;

            switch (cN.type) {

            case Node<ctype>::INTERSECTION_NODE:
                continue;

            case Node<ctype>::INTERIOR_NODE:
                
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

            case Node<ctype>::TOUCHING_NODE:
            case Node<ctype>::CORNER_NODE:
                

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
bool PSurface<dim,ctype>::map(int triIdx, const StaticVector<ctype,2>& p, std::tr1::array<int,3>& vertices, 
                         StaticVector<ctype,2>& coords, int seed) const
{
    int i;
    const DomainTriangle<ctype>& tri = this->triangles(triIdx);
    const std::vector<StaticVector<ctype,3> >& nP = iPos;

    // this is boundary handling
    if (p[0] < 0.001){
        for (i=0; i<tri.edgePoints[1].size()-1; i++){

            const StaticVector<ctype,2>& a = tri.nodes[tri.edgePoints[1][i]].domainPos();
            const StaticVector<ctype,2>& b = tri.nodes[tri.edgePoints[1][i+1]].domainPos();

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
            const StaticVector<ctype,2>& a = tri.nodes[tri.edgePoints[2][i]].domainPos();
            const StaticVector<ctype,2>& b = tri.nodes[tri.edgePoints[2][i+1]].domainPos();

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
            const StaticVector<ctype,2>& a = tri.nodes[tri.edgePoints[0][i]].domainPos();
            const StaticVector<ctype,2>& b = tri.nodes[tri.edgePoints[0][i+1]].domainPos();

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
        return false;

    StaticVector<ctype,3> imagePos = PlaneParam<ctype>::template linearInterpol<StaticVector<ctype,3> >(coords, nP[tri.nodes[v[0]].getNodeNumber()],
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
    const DomainTriangle<ctype>& cT = this->triangles(tri);
    //cT.print(true, true, true);
    int mode = cT.nodes[nds[0]].isINTERSECTION_NODE() +
        2*cT.nodes[nds[1]].isINTERSECTION_NODE() +
        4*cT.nodes[nds[2]].isINTERSECTION_NODE();
    //printf("MODE %d \n", mode);
    for (int i=0; i<3; i++)
        vertices[i] = getOtherEndNode(tri, nds[i]);
    //printf("***************MODE %d \n", mode);
    if (mode==6) {

        if (this->triangles(vertices[1].tri).nodes[vertices[1].idx].getNodeNumber() == 
            this->triangles(vertices[2].tri).nodes[vertices[2].idx].getNodeNumber()) {

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

        if (this->triangles(vertices[1].tri).nodes[vertices[1].idx].getNodeNumber() == 
            this->triangles(vertices[0].tri).nodes[vertices[0].idx].getNodeNumber()) {

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
    const Node<ctype>& cN = this->triangles(n.tri).nodes[n.idx];
    const ctype eps = 1e-6;

    switch (cN.type) {
    case Node<ctype>::GHOST_NODE: {
        
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
    case Node<ctype>::INTERSECTION_NODE:
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
void PSurface<dim,ctype>::handleMapOnEdge(int triIdx, const StaticVector<ctype,2>& p, const StaticVector<ctype,2>& a, const StaticVector<ctype,2>& b,
                                      int edge, int edgePos, std::tr1::array<GlobalNodeIdx, 3>& vertices, StaticVector<ctype,2>& coords) const
{
    const DomainTriangle<ctype>& tri = this->triangles(triIdx);
    ctype lambda = (p-a).length() / (a-b).length();

    StaticVector<ctype,3> targetPos = PlaneParam<ctype>::template linearInterpol<StaticVector<ctype,3> >(lambda, 
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
        
        typename PlaneParam<ctype>::DirectedEdgeIterator edge = tri.getDirectedEdgeIterator(n1, n2);
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
bool PSurface<dim,ctype>::positionMap(int triIdx, const StaticVector<ctype,2>& p, StaticVector<ctype,3>& result) const
{
    StaticVector<ctype,2> localCoords;
    std::tr1::array<int,3> tri;

    int status = map(triIdx, p, tri, localCoords);
    
    if (!status) {
        printf("p: (%f %f)\n", p[0], p[1]);
        this->triangles(triIdx).print(true, true, false);
        assert(false);
        return false;
    }

    result = PlaneParam<ctype>::template linearInterpol<StaticVector<ctype,3> >(localCoords, iPos[tri[0]], iPos[tri[1]], iPos[tri[2]]);

    return true;
}



template <int dim, class ctype>
bool PSurface<dim,ctype>::directNormalMap(int triIdx, const StaticVector<ctype,2>& p, StaticVector<ctype,3>& result) const
{
    StaticVector<ctype,2> localCoords;
    std::tr1::array<int,3> tri;

    int status = map(triIdx, p, tri, localCoords);
    
    if (!status)
        return false;

    const StaticVector<ctype,3> a = iPos[tri[1]] - iPos[tri[0]];
    const StaticVector<ctype,3> b = iPos[tri[2]] - iPos[tri[0]];
    result = a.cross(b);
    result.normalize();

    assert(!std::isnan(result[0]) && !std::isnan(result[1]) && !std::isnan(result[2]));

    return true;
}

template <int dim, class ctype>
int PSurface<dim,ctype>::invertTriangles(int patch)
{
    
    int i;
    int count=0;
    
    for (i=int(0); i<int(this->getNumTriangles()); i++) 
        if (patch==-1 || this->triangles(i).patch==patch){
            
            this->triangles(i).flip();
            count++;
            
            if (hasUpToDatePointLocationStructure) {
                
                for (int j=0; j<this->triangles(i).nodes.size(); j++)
                    this->triangles(i).nodes[j].reverseNeighbors();
            }
            
        }
    
    return count;
}


template <int dim, class ctype>
NodeBundle PSurface<dim,ctype>::getNodeBundleAtVertex(int v) const
{
    NodeBundle result;
    std::vector<int> neighbors = this->getTrianglesPerVertex(v);

    result.resize(neighbors.size());

    for (int i=0; i<neighbors.size(); i++) {

        result[i].tri = neighbors[i];
        const DomainTriangle<ctype>& cT = this->triangles(neighbors[i]);
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


template <int dim, class ctype>
void PSurface<dim,ctype>::checkConsistency(const char* where) const
{
#ifndef NDEBUG
    int i, j;
    // first checks whether all triangles are consistent
    // sorts out invalid triangles
    std::vector<bool> isInvalid(this->triangleArray.size());
    std::fill(isInvalid.begin(), isInvalid.end(), false);

    for (i=0; i<this->freeTriangleStack.size(); i++)
        isInvalid[this->freeTriangleStack[i]] = true;

    for (i=0; i<this->getNumTriangles(); i++)
        if (!isInvalid[i]){
            //printf("i = %d\n", i);
            this->triangles(i).checkConsistency("where");

            for (j=0; j<3; j++) {
                const Edge& cE = this->edges(this->triangles(i).edges[j]);
                if (!(cE.from == this->triangles(i).vertices[j] && cE.to   == this->triangles(i).vertices[(j+1)%3]) &&
                    !(cE.to   == this->triangles(i).vertices[j] && cE.from == this->triangles(i).vertices[(j+1)%3])){
                    printf(where);
                    printf("inconsistent triangle edges\n");
                    assert(false);
                }
            }
        }

    // checks whether matching edgepoint arrays have the same size
    for (i=0; i<this->getNumEdges(); i++) {
        const Edge& cE = this->edges(i);

        if (cE.triangles.size()!=2)
            continue;

        const DomainTriangle<ctype>& tri1 = this->triangles(cE.triangles[0]);
        const DomainTriangle<ctype>& tri2 = this->triangles(cE.triangles[1]);

        if (tri1.edgePoints[tri1.getEdge(i)].size() != tri2.edgePoints[tri2.getEdge(i)].size()) {
            printf(where);
            printf("Nonmatching edgePoint arrays at edge %d  (%d vs. %d)!\n", i,
                   tri1.edgePoints[tri1.getEdge(i)].size(),
                   tri2.edgePoints[tri2.getEdge(i)].size());
            tri1.print(true);
            tri2.print(true);
            assert(false);
        }
    }


#endif
}



// ////////////////////////////////////////////////////////
//   Explicit template instantiations.
//   If you need more, you can add them here.
// ////////////////////////////////////////////////////////

template class PSURFACE_EXPORT PSurface<2,float>;
template class PSURFACE_EXPORT PSurface<2,double>;

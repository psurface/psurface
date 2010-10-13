#ifdef _MSC_VER
    // Required to make cmath define M_PI etc.
    #define _USE_MATH_DEFINES
#endif
#include <cmath>

#include <limits>

#include <psurface/PlaneParam.h>
#include <psurface/StaticMatrix.h>
#include <psurface/SparseMatrix.h>

#ifdef _WIN32
inline int random() {return rand();}
#endif

template <class ctype>
typename PlaneParam<ctype>::DirectedEdgeIterator PlaneParam<ctype>::BFLocate(const StaticVector<ctype,2> &p, int seed) const
{

    //printf("----- BFLocate -----\n");
    int abort=0;

    DirectedEdgeIterator cE;

    if (seed<0 || seed>nodes.size()-1)
        cE = firstDirectedEdge();
    else
        cE = firstDirectedEdge(seed);

    if (orientation(cE, p) == -1){
        cE.invert();
    }

    while (true){

        StaticVector<ctype,2> from = nodes[cE.from()].domainPos();
        StaticVector<ctype,2> to   = nodes[cE.to()].domainPos();

        //printf("cE:  %d (%f %f) --> %d (%f %f)\n", cE.from(), from.x, from.y, cE.to(), to.x, to.y);

        if (abort++ >20000){
            printf("loop found trying to map (%f %f)\n", p[0], p[1]);
            cE.fromNode = -1;
            return cE;

        } else{

            int whichop = 0;

            DirectedEdgeIterator Onext = cE.getONext();
            DirectedEdgeIterator Dprev = cE.getDPrev();

            // test for boundary edge
            if (Onext.to() != Dprev.from()){

                printf("cE:  %d --> %d\n", cE.from(), cE.to());
                printf("Onext (%d) != Dprev(%d)\n", Onext.to(), Dprev.from());
                printf("p = (%f %f)\n", p[0], p[1]);
                cE.fromNode = -1;
                return cE;

            }else{

//                 printf("Onext:  %d --> %d;    Dprev:  %d --> %d\n",
//                        Onext.from(), Onext.to(), Dprev.from(), Dprev.to());
//                  nodes[Onext.from()].print();
//                  nodes[Onext.to()].print();
//                  nodes[Dprev.from()].print();
//                  nodes[Dprev.to()].print();

                bool isEligibleOnext = !(nodes[Onext.from()].isOnEdge(0) && nodes[Onext.to()].isOnEdge(0)) &&
                    !(nodes[Onext.from()].isOnEdge(1) && nodes[Onext.to()].isOnEdge(1)) &&
                    !(nodes[Onext.from()].isOnEdge(2) && nodes[Onext.to()].isOnEdge(2));
                bool isEligibleDprev = !(nodes[Dprev.from()].isOnEdge(0) && nodes[Dprev.to()].isOnEdge(0)) &&
                    !(nodes[Dprev.from()].isOnEdge(1) && nodes[Dprev.to()].isOnEdge(1)) &&
                    !(nodes[Dprev.from()].isOnEdge(2) && nodes[Dprev.to()].isOnEdge(2));

//                 printf("Onext: %d  (%d %d)  (%d %d) (%d %d)\n", isEligibleOnext,
//                         nodes[Onext.from()].isOnEdge(0), nodes[Onext.to()].isOnEdge(0),
//                         nodes[Onext.from()].isOnEdge(1), nodes[Onext.to()].isOnEdge(1),
//                         nodes[Onext.from()].isOnEdge(2), nodes[Onext.to()].isOnEdge(2));
//                 printf("Dprev: %d  (%d %d)  (%d %d) (%d %d)\n", isEligibleDprev,
//                         nodes[Dprev.from()].isOnEdge(0), nodes[Dprev.to()].isOnEdge(0),
//                         nodes[Dprev.from()].isOnEdge(1), nodes[Dprev.to()].isOnEdge(1),
//                         nodes[Dprev.from()].isOnEdge(2), nodes[Dprev.to()].isOnEdge(2));


                if (orientation(Onext, p) != -1 &&
                    // handle degenerate triangles
                    isEligibleOnext)
                    whichop += 1;
                
                if (orientation(Dprev, p) != -1 &&
                    // handle degenerate triangles
                    isEligibleDprev)
                    whichop += 2;

                //printf("whichop: %d\n", whichop);
                switch (whichop){
                    
                case 0: 
                    return cE;
                    
                case 1:
                    cE = Onext;
                    break;
                    
                case 2:
                    cE = Dprev;
                    break;
                    
                case 3:
                    cE = (random() < RAND_MAX/2) ? Onext : Dprev;
                    break;
                }
            }
        }

    }

}


//////////////////////////////////////////////////////////////////
// This routine sorts the neighbors of a given vertex in a cyclic order.  It is not as
// robust as the topological algorithm, however, adjacent neighbors need not be connected
// by edges.  Therefore, the method can be called without calling insertExtraEdges first.
template <class ctype>
void PlaneParam<ctype>::makeCyclicGeometrically(Node<ctype>& center)
{
    if (center.degree()<=2)
        return;

    int i, j;
   
    std::vector<ctype> angles(center.degree());

    // compute angles
    StaticVector<ctype,2> edge0Vec = nodes[center.neighbors(0)].domainPos() - center.domainPos();
    StaticVector<ctype,2> normal   = StaticVector<ctype,2>(-edge0Vec[1], edge0Vec[0]);

    for (i=0; i<center.degree(); i++){

        StaticVector<ctype,2> cEVec  = nodes[center.neighbors(i)].domainPos() - center.domainPos();

        ctype x = cEVec.dot(edge0Vec);
        ctype y = cEVec.dot(normal);

        angles[i] = atan2(y, x);
        if (angles[i]<0)
            angles[i] += 2*M_PI;
    }

    // bubblesort

    for (i=center.degree(); i>1; i--){
        bool swapped = false;

        for (j=0; j<i-1; j++){

            if (angles[j] > angles[j+1]){
                swapped = true;
                std::swap(angles[j], angles[j+1]);
                center.swapNeighbors(j, j+1);
            }
        }

        if (!swapped)
            break;
    }

}

// makeCyclic and DFSVisit sort the star of a node, that is the list of all direct neighbors
// in a cyclic order.  This is done by performing a depth-first search on the graph of these
// neighbors and looking for a longest path.  See me for details
template <class ctype>
bool PlaneParam<ctype>::DFSVisit(const std::vector<typename Node<ctype>::NeighborReference> &star, 
                                 const typename Node<ctype>::NeighborReference& u, 
                                 std::vector<typename Node<ctype>::NeighborReference> &outStar)
{
    int i, j;
    
    for (i=0; i<star.size(); i++){
        if (!nodes[u].isConnectedTo(star[i])) continue;
        const typename Node<ctype>::NeighborReference& v = star[i];

        // a cycle?
        bool isNew = true;
        for (j=0; j<outStar.size(); j++)
            if (outStar[j]==v){
                isNew=false;
                break;
            }

        if (isNew){
            outStar.push_back(v);
            if (outStar.size()==star.size() && nodes[outStar.back()].isConnectedTo(outStar[0])) 
                return true;
            if (DFSVisit(star, v, outStar)) 
                return true;
            outStar.pop_back();
        }
    }

    return false;
}

// makeCyclic and DFSVisit sort the star of a neuron, that is the list of all direct neighbors
// in a cyclic order.  This is done by performing a depth-first search on the graph of these
// neighbors and looking for a longest path.  See me for details
// This is for the boundary case

// The parameter u needs to be handed over by value, because this method is used with
// u = outStar[0].  However the outStar.push_back in this method may lead to a relocation
// of the outStar content, and the reference to outStar[0] (in u) will dangle.
template <class ctype>
bool PlaneParam<ctype>::DFSBoundaryVisit(const std::vector<typename Node<ctype>::NeighborReference> &star, 
                                         typename Node<ctype>::NeighborReference u, int endNode,
                                         std::vector<typename Node<ctype>::NeighborReference> &outStar)
{
    for (int i=0; i<star.size(); i++){

        if (!nodes[u].isConnectedTo(star[i])) continue;
        const typename Node<ctype>::NeighborReference& v = star[i];

        // no cycle yet
        if (std::find(outStar.begin(), outStar.end(), v) == outStar.end()) {

            outStar.push_back(v);
            if (outStar.size()==star.size() && outStar.back()==endNode) 
                return true;
            if (DFSBoundaryVisit(star, v, endNode, outStar)) 
                return true;
            outStar.pop_back();

        }
    }

    return false;
}


template <class ctype>
void PlaneParam<ctype>::makeCyclicInteriorNode(Node<ctype> &center)
{
    std::vector<typename Node<ctype>::NeighborReference> outStar(1);
    outStar[0] = center.neighbors(0);

    if (!DFSVisit(center.nbs, center.neighbors(0), outStar)) {  // if not -> programming error
        printf("DFSVisit failed!\n");
        assert(false);
    }

    center.nbs = outStar;

    // the neighbors are now sorted in cyclic order.  But is the orientation correct?
    // The orientation needs to be consistent for all vertices in order to get
    // correctly oriented normals

    StaticVector<ctype,2> referenceVector = nodes[center.neighbors(0)].domainPos() - center.domainPos();
    StaticVector<ctype,2> normal          = StaticVector<ctype,2>(-referenceVector[1], referenceVector[0]);
    int i;
    int leastPosVector = -1;
    int mostPosVector  = -1;
    ctype maxDotProdukt = -std::numeric_limits<ctype>::max();
    ctype minDotProdukt = std::numeric_limits<ctype>::max();

    for (i=1; i<center.degree(); i++) {

        StaticVector<ctype,2> testVector = nodes[center.neighbors(i)].domainPos() - center.domainPos();
        if (testVector.dot(normal) > maxDotProdukt) {
            maxDotProdukt = testVector.dot(normal);
            mostPosVector = i;
        }

        if (testVector.dot(normal) < minDotProdukt) {
            minDotProdukt = testVector.dot(normal);
            leastPosVector = i;
        }

    }

    if (leastPosVector < mostPosVector)
        center.reverseNeighbors();

}

/////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////


template <class ctype>
StaticVector<ctype,2> PlaneParam<ctype>::computeBarycentricCoords(const StaticVector<ctype,2> &p, const StaticVector<ctype,2> &a, const StaticVector<ctype,2> &b, const StaticVector<ctype,2> &c)
{
    StaticVector<ctype,2> result;
    
    StaticMatrix<ctype,3> area0(p[0], b[0], c[0],  p[1], b[1], c[1], 1, 1, 1);
    StaticMatrix<ctype,3> area1(a[0], p[0], c[0],  a[1], p[1], c[1], 1, 1, 1);
    
    StaticMatrix<ctype,3> areaTotal(a[0], b[0], c[0], a[1], b[1], c[1], 1, 1, 1);
    ctype areaTotalDet = areaTotal.det();
    
    result[0] = area0.det()/areaTotalDet;
    result[1] = area1.det()/areaTotalDet;
    
    return result;
}

// This routine computes the barycentric coordinates of a point in space with respect to
// a triangle in space.  It tacitly assumes that the point is coplanar with the triangle.
template <class ctype>
StaticVector<ctype,2> PlaneParam<ctype>::computeBarycentricCoords(const StaticVector<ctype,3> &p, const StaticVector<ctype,3> &a, const StaticVector<ctype,3> &b, const StaticVector<ctype,3> &c)
{
    StaticVector<ctype,2> result;
    
    ctype area0 = (p-b).cross(p-c).length();
    ctype area1 = (p-a).cross(p-c).length();
    
    ctype areaTotal = (b-a).cross(c-a).length();
    
    result[0] = area0/areaTotal;
    result[1] = area1/areaTotal;
    
    if (std::isnan(result[1])) {
        printf("area0 %f   area1 %f    areaTotal %f   res  (%f %f)\n", area0, area1, areaTotal, 
               result[0], result[1]);
        assert(false);
    }

    return result;
}

template <class ctype>
int PlaneParam<ctype>::map(StaticVector<ctype,2> &domainCoord, std::tr1::array<NodeIdx, 3>& tri, StaticVector<ctype,2>& localBarycentricCoords,
                    int seed) const
{
    DirectedEdgeIterator e = BFLocate(domainCoord);

    if (!e.isValid()) {
        printf("[PlaneParam::map] An error occured when calling BFLocate\n");
        return false;
    }

    // test for boundary ParameterEdge
    DirectedEdgeIterator oNext = e.getONext();
    DirectedEdgeIterator dPrev = e.getDPrev();

    if (oNext.to() != dPrev.from()){
        e.invert();
        oNext = e.getONext();
    }

    tri[0] = e.from();
    tri[1] = e.to();
    tri[2] = oNext.to();
    

    localBarycentricCoords = computeBarycentricCoords(domainCoord, 
                                                      nodes[tri[0]].domainPos(),
                                                      nodes[tri[1]].domainPos(),
                                                      nodes[tri[2]].domainPos());
    
    if (localBarycentricCoords[0]<-0.05 || localBarycentricCoords[1] < -0.05 ||
        (localBarycentricCoords[0]+localBarycentricCoords[1] > 1.05)) {
        printf("There seems to be a self-intersection in your parametrization.\n");
        printf("You should try to smooth it and retry.\n");
        printf("localBarycentricCoords: (%f %f)\n", localBarycentricCoords[0], localBarycentricCoords[1]);
        return false;
    }

    return true;
}


template <class ctype>
void PlaneParam<ctype>::unflipTriangles(const std::vector<StaticVector<ctype,3> >& nodePositions)
{
    applyParametrization(nodePositions);
    return;

    // a plane triangulation contains flipped triangles if at least one
    // of its vertices is not a convex combination of its neighbors
//     Node* cN;
    
//     for (cN=nodes.first(); cN; cN=nodes.succ(cN))
//      if (cN->type==Node::INTERIOR_NODE && !cN->isConvexCombination()){
//          //printf("unflipping!\n");
//          applyParametrization(0);
//          return;
//      }
            
    //printf("NOT unflipping!\n");
}

////////////////////////////////////////////////////////////////
// this routine installs the shape-preserving parametrization 
// only INTERIOR_NODEs get moved
////////////////////////////////////////////////////////////////
template <class ctype>
void PlaneParam<ctype>::applyParametrization(const std::vector<StaticVector<ctype,3> >& nodePositions)
{
    int i;
    
    // compute lambdas
    SparseMatrix<ctype> lambda_ij(nodes.size());

    computeFloaterLambdas(lambda_ij, nodePositions);
    
    // build matrix 
    lambda_ij *= -1;

    for (i=0; i<lambda_ij.nRows(); i++)
        lambda_ij.setEntry(i, i, 1);
    
    // Compute the right side. We use complex numbers for solving the systems
    // for both x- and y-components in one pass.  This leads to a considerable
    // speedup.
    std::vector<std::complex<ctype> > b(nodes.size());
    
    std::fill(b.begin(), b.end(), 0);
    
    for (i=0; i<nodes.size(); i++) 
        if (!nodes[i].isINTERIOR_NODE()) {
            // not elegant
            b[i] = std::complex<ctype>(nodes[i].domainPos()[0], nodes[i].domainPos()[1]);
        }

    // solve the system
    int maxIter=3000;
    std::vector<std::complex<ctype> > residue;
    std::vector<std::complex<ctype> > result(nodes.size());
    
    for (i=0; i<nodes.size(); i++)
        result[i] = std::complex<ctype>(nodes[i].domainPos()[0], nodes[i].domainPos()[1]);
    
    lambda_ij.BiCGSTABC(b, result, residue, &maxIter, 1e-6);
    
    for (i=0; i<nodes.size(); i++)
        if (nodes[i].isINTERIOR_NODE())
            nodes[i].setDomainPos(StaticVector<ctype,2>(real(result[i]), imag(result[i])));

}


            
            
////////////////////////////////////////////////////////
// computes lambda_ij for the Floater-Parametrization
////////////////////////////////////////////////////////
template <class ctype>
void PlaneParam<ctype>::computeFloaterLambdas(SparseMatrix<ctype>& lambda_ij, 
                                       const std::vector<StaticVector<ctype,3> >& nodePositions)
{
    int i, k, l;
    int N = nodes.size();

    assert(lambda_ij.nRows()==N && lambda_ij.nCols()==N);

    // init lambda array

    // for all interiorPoints do
    for (i=0; i<nodes.size(); i++) {
        if (nodes[i].isINTERIOR_NODE()) {
            
            Node<ctype>& p = nodes[i];
            makeCyclicGeometrically(p);
            
            std::vector<int>    p_k(p.degree());
            std::vector<StaticVector<ctype,3> >  p_k_3DCoords(p.degree());
            std::vector<StaticVector<ctype,2> >  p_k_2DCoords(p.degree());
            std::vector<ctype>    angle;
            
            for (k=0; k<p.degree(); k++){
                p_k[k]          = (int)p.neighbors(k);
                //p_k_3DCoords[k] = nodes[p_k[k]].getImagePos(nodePositions);
                p_k_3DCoords[k] = nodePositions[nodes[p_k[k]].getNodeNumber()];
                if (std::isnan(p_k_3DCoords[k][0])) {
                    printf("iPos.size: %d,  nN: %d\n", nodePositions.size(), nodes[p_k[k]].getNodeNumber());
                    nodes[p_k[k]].print();
                }
                assert(!std::isnan(p_k_3DCoords[k][0]));
                assert(!std::isnan(p_k_3DCoords[k][1]));
                assert(!std::isnan(p_k_3DCoords[k][2]));
            }
            
            if (!polarMap(nodePositions[p.getNodeNumber()], p_k_3DCoords, p_k_2DCoords, angle )) {

                for (k=0; k<p.degree(); k++)
                    lambda_ij.setEntry(i, p.neighbors(k), 1/((ctype)p.degree()));

                continue;
            }
            
            if (p.degree()==3){
                
                StaticVector<ctype,2> lambdas = computeBarycentricCoords(StaticVector<ctype,2>(0,0), p_k_2DCoords[0], p_k_2DCoords[1], p_k_2DCoords[2]);
                
                StaticVector<ctype,3> l_ij;
                l_ij[0] = lambdas[0];
                l_ij[1] = lambdas[1];
                l_ij[2] = 1-lambdas[0]-lambdas[1];
                
                for (k=0; k<3; k++)
                    lambda_ij.setEntry(i, p_k[k], l_ij[k]);
                
            } else {
                
                std::vector<int> index(p.degree());
                for (l=0; l<p.degree(); l++) 
                    index[l] = p_k[l];
                
                for (l=0; l<p.degree(); l++) {
                    
                    int rlPlus1=0;
                    
                    ctype oppositeAngle = (angle[l]<M_PI) ? angle[l]+M_PI : angle[l]-M_PI;
                    
                    while (rlPlus1<angle.size() && angle[rlPlus1] < oppositeAngle) rlPlus1++;
                    
                    int rl = (rlPlus1 + p.degree()-1)%p.degree();
                    rlPlus1 = rlPlus1%p.degree();
                    
                    StaticVector<ctype,2> bCoords = computeBarycentricCoords(StaticVector<ctype,2>(0,0), p_k_2DCoords[l], p_k_2DCoords[rl], p_k_2DCoords[rlPlus1]);
                    
                    StaticVector<ctype,3> delta(bCoords[0], bCoords[1], 1-bCoords[0]-bCoords[1]);
                    
                    lambda_ij.addToEntry(i, index[l],       delta[0] / p.degree());
                    lambda_ij.addToEntry(i, index[rl],      delta[1] / p.degree());
                    lambda_ij.addToEntry(i, index[rlPlus1], delta[2] / p.degree());
                }
                
            }
            
        }
    }
}
            

template <class ctype>
bool PlaneParam<ctype>::polarMap(const StaticVector<ctype,3>& center, const std::vector<StaticVector<ctype,3> > &threeDStarVertices, 
                          std::vector<StaticVector<ctype,2> >& flattenedCoords, std::vector<ctype >& theta)
{
    /////////////////////////////////////
    // computes the flattened coordinates
    const int K = threeDStarVertices.size();

    flattenedCoords.resize(K);

    theta.resize(K+1);

    // compute the (accumulated) angles at the center point
    theta[0] = 0;

    int k;

    for (k=1; k<K+1; k++){
        StaticVector<ctype,3> pLeft  = threeDStarVertices[k-1];
        StaticVector<ctype,3> pRight = threeDStarVertices[k%K];

        if ( (pLeft-center).length()==0 || (pRight-center).length()==0){
            printf("vertex coincides with its neighbor, aborting polar map\n");
            return false;
        }

        theta[k] = theta[k-1] + (pLeft - center).angle(pRight - center);
        if (std::isnan(theta[k])){
            printf("center (%f %f %f)\n", center[0], center[1], center[2]);
            printf("pLeft - center (%f %f %f) pRight - center (%f %f %f)\n", 
                   pLeft[0] - center[0], pLeft[1] - center[1], pLeft[2] - center[2], 
                   pRight[0] - center[0], pRight[1] - center[1], pRight[2] - center[2]);
            printf("pLeft (%f %f %f)   pRight(%f %f %f)\n", pLeft[0], pLeft[1], pLeft[2], 
                   pRight[0], pRight[1], pRight[2]);

            printf("angle %f\n", (pLeft - center).angle(pRight - center));
            return false;
        }
        
    }

    ctype a = 2*M_PI/theta[K];

    // compute parameter domain coordinates
    for (k=0; k<K; k++){

        theta[k] *= a;
        ctype r = (threeDStarVertices[k] - center).length();
        ctype rPowA = powf(r, a);

        flattenedCoords[k] = StaticVector<ctype,2>(rPowA*cos(theta[k]), rPowA*sin(theta[k])); 
    }

    theta.pop_back();

    return true;
}

template <class ctype>
void PlaneParam<ctype>::removeExtraEdges()
{
    checkConsistency("before removing of extra edges");
    for (int i=0; i<nodes.size(); i++)
        for (int j=nodes[i].degree()-1; j>=0; j--)
            if (!nodes[i].neighbors(j).isRegular())
                nodes[i].removeNeighbor(j);

    checkConsistency("after removing of extra edges");

}

template <class ctype>
void PlaneParam<ctype>::makeCyclicBoundaryNode(Node<ctype>& center, int next, int previous)
{

//     printf("------------------------------------\n");
//     center.print();
//     printf("next %d    previous %d\n", next, previous);

    int i;
    std::vector<typename Node<ctype>::NeighborReference> outStar(1);

    // look for the correct NeighborReference pointing to #next#

    for (i=0; i<center.degree(); i++) {
        if (center.neighbors(i)==next) {
            outStar[0] = center.neighbors(i);  // #next# is an int, but center.neighbors are NeighborReferences!!
            break;
        }
    }

    assert(i<center.degree());

    if (!DFSBoundaryVisit(center.nbs, outStar[0], previous, outStar)) { // if not -> programming error
        printf("DFSBoundaryVisit failed!\n");
        center.print();printf("\n");
        for (i=0; i<center.degree(); i++){
            printf("### number %d\n", (int)center.neighbors(i));
            nodes[center.neighbors(i)].print();
        }
        
        //assert(false);
    }
    
    center.nbs = outStar;
}

template <class ctype>
void PlaneParam<ctype>::installWorldCoordinates(const StaticVector<ctype,2> &a, const StaticVector<ctype,2> &b, const StaticVector<ctype,2> &c)
{
    for (int i=0; i<nodes.size(); i++)
        nodes[i].setDomainPos(a*nodes[i].domainPos()[0] + b*nodes[i].domainPos()[1] + 
                              c*(1-nodes[i].domainPos()[0]-nodes[i].domainPos()[1]));
}


template <class ctype>
void PlaneParam<ctype>::installBarycentricCoordinates(const StaticVector<ctype,2> &a, const StaticVector<ctype,2> &b, const StaticVector<ctype,2> &c)
{
    for (int i=0; i<nodes.size(); i++) {
        //printf("node %d,  before (%f %f)  ", i, nodes[i].domainPos()[0], nodes[i].domainPos()[1]);
        nodes[i].setDomainPos(computeBarycentricCoords(nodes[i].domainPos(), a, b, c));
        //printf("after (%f %f) \n ", nodes[i].domainPos()[0], nodes[i].domainPos()[1]);
    }
}

template <class ctype>
void PlaneParam<ctype>::print(bool showNodes, bool showParamEdges, bool showExtraEdges) const 
{
    std::cout << "---------------------------------------------------------" << std::endl;
    std::cout << "parametrization contains " << nodes.size() << " nodes" << std::endl;
    
    if (showNodes){
        for (size_t i=0; i<nodes.size(); i++)
            nodes[i].print();
    }

    std::cout << "---------------------------------------------------------" << std::endl;
}   


template <class ctype>
void PlaneParam<ctype>::checkConsistency(const char* where) const
{
#ifndef NDEBUG
    for (size_t i=0; i<nodes.size(); i++) {

        const Node<ctype>& cN = nodes[i];
        if (cN.isInvalid())
            continue;

        if (std::isnan(cN.domainPos()[0]) || std::isnan(cN.domainPos()[1])) {
            printf(where);
            printf("\n****** A node mit NaN domainPos found! ******\n");
            cN.print();
            assert(false);
        }

        // make sure references are mutual
        for (size_t j=0; j<cN.degree(); j++)
            if (!nodes[cN.neighbors(j)].isConnectedTo(i)) {
                printf(where);
                printf("\n***** Neighbor relation is not mutual j=%d   k=%d *****\n", j, i);
                cN.print();
                nodes[cN.neighbors(j)].print();
                assert(false);
            }
        
        // make sure that no neighbor is invalid
        for (size_t j=0; j<cN.degree(); j++)
            if (nodes[cN.neighbors(j)].isInvalid()) {
                printf(where);
                printf("***** Node has an invalid neighbor *****\n");
                assert(false);
            }

        // check for double edges
        for (size_t l=0; l<cN.degree(); l++)
            for (size_t j=0; j<l; j++)
                if (cN.neighbors(l)==cN.neighbors(j)) {
                    printf(where);
                    printf("***** PlaneParam contains double edge! *****\n");
                    for (size_t k=0; k<cN.degree(); k++){
                        printf("   %d\n  ", (int)cN.neighbors(k));
                        nodes[cN.neighbors(k)].print();
                    }

                    cN.print();
                    assert(false);
                }

        if (!cN.degree() && !cN.isCORNER_NODE()){
            printf(where);
            cN.print();
            printf("NodeNumber = %d\n", i);
            printf("****** solitary Node found!\n");
            assert(false);
        }

        for (int i=0; i<cN.degree(); i++)
            assert(cN.neighbors(i)>=0 && cN.neighbors(i)<nodes.size());

    }

#endif
}


// ////////////////////////////////////////////////////////
//   Explicit template instantiations.
//   If you need more, you can add them here.
// ////////////////////////////////////////////////////////

template class PlaneParam<float>;
template class PlaneParam<double>;

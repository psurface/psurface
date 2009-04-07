#include <psurface/PlaneParam.h>
#include <psurface/StaticMatrix.h>
#include <mclib/McSparseMatrix.h>
#include <mclib/McDVector.h>

#include <limits>

#ifdef _WIN32
//#include <float.h>
inline int isnan(double x) {return _isnan(x);}
inline int random() {return rand();}
#endif

PlaneParam::~PlaneParam() {}

PlaneParam::DirectedEdgeIterator PlaneParam::BFLocate(const StaticVector<float,2> &p, int seed) const
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

        StaticVector<float,2> from = nodes[cE.from()].domainPos();
        StaticVector<float,2> to   = nodes[cE.to()].domainPos();

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
void PlaneParam::makeCyclicGeometrically(Node& center)
{
    if (center.degree()<=2)
        return;

    int i, j;
   
    McSmallArray<float, 12> angles(center.degree());

    // compute angles
    StaticVector<float,2> edge0Vec = nodes[center.neighbors(0)].domainPos() - center.domainPos();
    StaticVector<float,2> normal   = StaticVector<float,2>(-edge0Vec[1], edge0Vec[0]);

    for (i=0; i<center.degree(); i++){

        StaticVector<float,2> cEVec  = nodes[center.neighbors(i)].domainPos() - center.domainPos();

        float x = cEVec.dot(edge0Vec);
        float y = cEVec.dot(normal);

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
                angles.swap(j, j+1);
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
bool PlaneParam::DFSVisit(const McSmallArray<Node::NeighborReference, 6> &star, const Node::NeighborReference& u, 
                              McSmallArray<Node::NeighborReference, 6> &outStar)
{
    int i, j;
    
    for (i=0; i<star.size(); i++){
        if (!nodes[u].isConnectedTo(star[i])) continue;
        const Node::NeighborReference& v = star[i];

        // a cycle?
        bool isNew = true;
        for (j=0; j<outStar.size(); j++)
            if (outStar[j]==v){
                isNew=false;
                break;
            }

        if (isNew){
            outStar.append(v);
            if (outStar.size()==star.size() && nodes[outStar.last()].isConnectedTo(outStar[0])) return true;
            if (DFSVisit(star, v, outStar)) return true;
            outStar.removeLast();
        }
    }

    return false;
}

// makeCyclic and DFSVisit sort the star of a neuron, that is the list of all direct neighbors
// in a cyclic order.  This is done by performing a depth-first search on the graph of these
// neighbors and looking for a longest path.  See me for details
// This is for the boundary case
bool PlaneParam::DFSBoundaryVisit(const McSmallArray<Node::NeighborReference, 6> &star, 
                                  const Node::NeighborReference& u, int endNode,
                                  McSmallArray<Node::NeighborReference, 6> &outStar)
{
    int i, j;

//     for (i=0; i<outStar.size(); i++)
//      printf("  %d", (int)outStar[i]);

//     printf("\n");
    
    for (i=0; i<star.size(); i++){
        //printf("i = %d   star.size = %d\n", i, star.size());
        if (!nodes[u].isConnectedTo(star[i])) continue;
        const Node::NeighborReference& v = star[i];

        // a cycle?
        bool isNew = true;
        for (j=0; j<outStar.size(); j++)
            if (outStar[j]==v){
                isNew=false;
                break;
            }

        if (isNew){
            outStar.append(v);
            if (outStar.size()==star.size() && outStar.last()==endNode) return true;
            if (DFSBoundaryVisit(star, v, endNode, outStar)) return true;
            outStar.removeLast();
        }
    }

    return false;
}


void PlaneParam::makeCyclicInteriorNode(Node &center)
{
    McSmallArray<Node::NeighborReference, 6> outStar(1);
    outStar[0] = center.neighbors(0);

    if (!DFSVisit(center.nbs, center.neighbors(0), outStar)) {  // if not -> programming error
        printf("DFSVisit failed!\n");
        assert(false);
    }

    center.nbs = outStar;

    // the neighbors are now sorted in cyclic order.  But is the orientation correct?
    // The orientation needs to be consistent for all vertices in order to get
    // correctly oriented normals

    StaticVector<float,2> referenceVector = nodes[center.neighbors(0)].domainPos() - center.domainPos();
    StaticVector<float,2> normal          = StaticVector<float,2>(-referenceVector[1], referenceVector[0]);
    int i;
    int leastPosVector = -1;
    int mostPosVector  = -1;
    float maxDotProdukt = -std::numeric_limits<float>::max();
    float minDotProdukt = std::numeric_limits<float>::max();

    for (i=1; i<center.degree(); i++) {

        StaticVector<float,2> testVector = nodes[center.neighbors(i)].domainPos() - center.domainPos();
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


StaticVector<float,2> PlaneParam::computeBarycentricCoords(const StaticVector<float,2> &p, const StaticVector<float,2> &a, const StaticVector<float,2> &b, const StaticVector<float,2> &c)
{
    StaticVector<float,2> result;
    
    StaticMatrix<float,3> area0(p[0], b[0], c[0],  p[1], b[1], c[1], 1, 1, 1);
    StaticMatrix<float,3> area1(a[0], p[0], c[0],  a[1], p[1], c[1], 1, 1, 1);
    
    StaticMatrix<float,3> areaTotal(a[0], b[0], c[0], a[1], b[1], c[1], 1, 1, 1);
    float areaTotalDet = areaTotal.det();
    
    result[0] = area0.det()/areaTotalDet;
    result[1] = area1.det()/areaTotalDet;
    
    return result;
}

// This routine computes the barycentric coordinates of a point in space with respect to
// a triangle in space.  It tacitly assumes that the point is coplanar with the triangle.

StaticVector<float,2> PlaneParam::computeBarycentricCoords(const StaticVector<float,3> &p, const StaticVector<float,3> &a, const StaticVector<float,3> &b, const StaticVector<float,3> &c)
{
    StaticVector<float,2> result;
    
    float area0 = (p-b).cross(p-c).length();
    float area1 = (p-a).cross(p-c).length();
    
    float areaTotal = (b-a).cross(c-a).length();
    
    result[0] = area0/areaTotal;
    result[1] = area1/areaTotal;
    
    if (isnan(result[1])) {
        printf("area0 %f   area1 %f    areaTotal %f   res  (%f %f)\n", area0, area1, areaTotal, 
               result[0], result[1]);
        assert(false);
    }

    return result;
}

int PlaneParam::map(StaticVector<float,2> &domainCoord, McSArray<NodeIdx, 3>& tri, StaticVector<float,2>& localBarycentricCoords,
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


void PlaneParam::unflipTriangles(const std::vector<StaticVector<float,3> >& nodePositions)
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
void PlaneParam::applyParametrization(const std::vector<StaticVector<float,3> >& nodePositions)
{
    int i;
    
    // compute lambdas
    McSparseMatrix<float, false> lambda_ij(nodes.size());

    computeFloaterLambdas(lambda_ij, nodePositions);
    
    // build matrix 
    lambda_ij *= -1;

    for (i=0; i<lambda_ij.nRows(); i++)
        lambda_ij.setEntry(i, i, 1);
    
    // Compute the right side. We use complex numbers for solving the systems
    // for both x- and y-components in one pass.  This leads to a considerable
    // speedup.
    McDVector<MC_complex<float> > b(nodes.size());
    
    b.fill(0);
    
    for (i=0; i<nodes.size(); i++) 
        if (!nodes[i].isINTERIOR_NODE()) {
            // not elegant
            b[i] = MC_complex<float>(nodes[i].domainPos()[0], nodes[i].domainPos()[1]);
        }

    // solve the system
    int maxIter=3000;
    McDVector<MC_complex<float> > residue;
    McDVector<MC_complex<float> > result(nodes.size());
    
    for (i=0; i<nodes.size(); i++)
        result[i] = MC_complex<float>(nodes[i].domainPos()[0], nodes[i].domainPos()[1]);
    
    lambda_ij.BiCGSTABC(b, result, residue, &maxIter, 1e-6);
    
    for (i=0; i<nodes.size(); i++)
        if (nodes[i].isINTERIOR_NODE())
            nodes[i].setDomainPos(StaticVector<float,2>(real(result[i]), imag(result[i])));

}


            
            
////////////////////////////////////////////////////////
// computes lambda_ij for the Floater-Parametrization
void PlaneParam::computeFloaterLambdas(McSparseMatrix<float, false>& lambda_ij, 
                                       const std::vector<StaticVector<float,3> >& nodePositions)
{
    int i, k, l;
    int N = nodes.size();

    assert(lambda_ij.nRows()==N && lambda_ij.nCols()==N);

    // init lambda array

    // for all interiorPoints do
    for (i=0; i<nodes.size(); i++) {
        if (nodes[i].isINTERIOR_NODE()) {
            
            Node& p = nodes[i];
            makeCyclicGeometrically(p);
            
            McSmallArray<int, 15>    p_k(p.degree());
            McSmallArray<StaticVector<float,3>, 15>  p_k_3DCoords(p.degree());
            McSmallArray<StaticVector<float,2>, 15>  p_k_2DCoords(p.degree());
            McSmallArray<float, 15>    angle;
            
            for (k=0; k<p.degree(); k++){
                p_k[k]          = (int)p.neighbors(k);
                //p_k_3DCoords[k] = nodes[p_k[k]].getImagePos(nodePositions);
                p_k_3DCoords[k] = nodePositions[nodes[p_k[k]].getNodeNumber()];
                if (isnan(p_k_3DCoords[k][0])) {
                    printf("iPos.size: %d,  nN: %d\n", nodePositions.size(), nodes[p_k[k]].getNodeNumber());
                    nodes[p_k[k]].print();
                }
                assert(!isnan(p_k_3DCoords[k][0]));
                assert(!isnan(p_k_3DCoords[k][1]));
                assert(!isnan(p_k_3DCoords[k][2]));
            }
            
            if (!polarMap(nodePositions[p.getNodeNumber()], p_k_3DCoords, p_k_2DCoords, angle )) {

                for (k=0; k<p.degree(); k++)
                    lambda_ij.setEntry(i, p.neighbors(k), 1/((float)p.degree()));

                continue;
            }
            
            if (p.degree()==3){
                
                StaticVector<float,2> lambdas = computeBarycentricCoords(StaticVector<float,2>(0,0), p_k_2DCoords[0], p_k_2DCoords[1], p_k_2DCoords[2]);
                
                StaticVector<float,3> l_ij;
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
                    
                    float oppositeAngle = (angle[l]<M_PI) ? angle[l]+M_PI : angle[l]-M_PI;
                    
                    while (rlPlus1<angle.size() && angle[rlPlus1] < oppositeAngle) rlPlus1++;
                    
                    int rl = (rlPlus1 + p.degree()-1)%p.degree();
                    rlPlus1 = rlPlus1%p.degree();
                    
                    StaticVector<float,2> bCoords = computeBarycentricCoords(StaticVector<float,2>(0,0), p_k_2DCoords[l], p_k_2DCoords[rl], p_k_2DCoords[rlPlus1]);
                    
                    StaticVector<float,3> delta(bCoords[0], bCoords[1], 1-bCoords[0]-bCoords[1]);
                    
                    lambda_ij.addToEntry(i, index[l],       delta[0] / p.degree());
                    lambda_ij.addToEntry(i, index[rl],      delta[1] / p.degree());
                    lambda_ij.addToEntry(i, index[rlPlus1], delta[2] / p.degree());
                }
                
            }
            
        }
    }
}
            


bool PlaneParam::polarMap(const StaticVector<float,3>& center, const McSmallArray<StaticVector<float,3>, 15> &threeDStarVertices, 
                          McSmallArray<StaticVector<float,2>, 15>& flattenedCoords, McSmallArray<float, 15>& theta)
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
        StaticVector<float,3> pLeft  = threeDStarVertices[k-1];
        StaticVector<float,3> pRight = threeDStarVertices[k%K];

        if ( (pLeft-center).length()==0 || (pRight-center).length()==0){
            printf("vertex coincides with its neighbor, aborting polar map\n");
            return false;
        }

        theta[k] = theta[k-1] + (pLeft - center).angle(pRight - center);
        if (isnan(theta[k])){
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

    float a = 2*M_PI/theta[K];

    // compute parameter domain coordinates
    for (k=0; k<K; k++){

        theta[k] *= a;
        float r = (threeDStarVertices[k] - center).length();
        float rPowA = powf(r, a);

        flattenedCoords[k] = StaticVector<float,2>(rPowA*cos(theta[k]), rPowA*sin(theta[k])); 
    }

    theta.removeLast();

    return true;
}

void PlaneParam::removeExtraEdges()
{
    checkConsistency("before removing of extra edges");
    for (int i=0; i<nodes.size(); i++)
        for (int j=nodes[i].degree()-1; j>=0; j--)
            if (!nodes[i].neighbors(j).isRegular())
                nodes[i].removeNeighbor(j);

    checkConsistency("after removing of extra edges");

}

void PlaneParam::makeCyclicBoundaryNode(Node& center, int next, int previous)
{

//     printf("------------------------------------\n");
//     center.print();
//     printf("next %d    previous %d\n", next, previous);

    int i;
    McSmallArray<Node::NeighborReference, 6> outStar(1);

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

void PlaneParam::installWorldCoordinates(const StaticVector<float,2> &a, const StaticVector<float,2> &b, const StaticVector<float,2> &c)
{
    for (int i=0; i<nodes.size(); i++)
        nodes[i].setDomainPos(a*nodes[i].domainPos()[0] + b*nodes[i].domainPos()[1] + 
                              c*(1-nodes[i].domainPos()[0]-nodes[i].domainPos()[1]));
}


void PlaneParam::installBarycentricCoordinates(const StaticVector<float,2> &a, const StaticVector<float,2> &b, const StaticVector<float,2> &c)
{
    for (int i=0; i<nodes.size(); i++) {
        //printf("node %d,  before (%f %f)  ", i, nodes[i].domainPos()[0], nodes[i].domainPos()[1]);
        nodes[i].setDomainPos(computeBarycentricCoords(nodes[i].domainPos(), a, b, c));
        //printf("after (%f %f) \n ", nodes[i].domainPos()[0], nodes[i].domainPos()[1]);
    }
}

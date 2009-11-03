#include <psurface/PlaneParam.h>
#include <psurface/Domains.h>
#include <psurface/DomainPolygon.h>
#include <psurface/CircularPatch.h>

template <class ctype>
void Node<ctype>::print(bool showNeighbors) const
{
#ifndef NDEBUG
    printf("dom (%f %f) ", domainPos()[0], domainPos()[1]);

    switch(type){
    case INTERIOR_NODE:
        printf("INTERIOR_NODE");
        break;
    case TOUCHING_NODE:
        printf("TOUCHING_NODE");
        break;
    case INTERSECTION_NODE:
        printf("INTERSECTION_NODE");
        break;
    case CORNER_NODE:
        printf("CORNER_NODE");
        break;
    case GHOST_NODE:
        printf("GHOST_NODE");
        break;
    }

    printf(" number %d", nodeNumber);

    if (isOnEdge())
        std::cout << "  edge: " << getDomainEdge() << "  edgePos " << getDomainEdgePosition() << std::endl;
    else if (isOnCorner())
        printf("  corner: %d\n", getCorner());
    else
        printf("\n");

    if (showNeighbors)
        for (int i=0; i<degree(); i++)
            printf("   %d %s\n", (int)nbs[i], nbs[i].isRegular() ? " " : "c");

#endif
}
 
template <class ctype>
void PlaneParam<ctype>::print(bool showNodes, bool showParamEdges, bool showExtraEdges) const 
{
#ifndef NDEBUG
    std::cout << "---------------------------------------------------------" << std::endl;
    std::cout << "parametrization contains " << nodes.size() << " nodes" << std::endl;
    
    if (showNodes){
        for (size_t i=0; i<nodes.size(); i++)
            nodes[i].print();
    }

    std::cout << "---------------------------------------------------------" << std::endl;
#endif
}   

void DomainPolygon::print(bool showEdgePoints, bool showParamEdges, bool showNodes) const 
{
#ifndef NDEBUG
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
#endif
}

template <class ctype>
void DomainTriangle<ctype>::print(bool showEdgePoints, bool showParamEdges, bool showNodes) const
{
#ifndef NDEBUG
    int i, j;

    printf("--------------------------------------------------------\n");
    printf("--  Print Triangle  ------------------------------------\n");

    printf("vertices:  (%d %d %d)\n", vertices[0], vertices[1], vertices[2]);
    printf("edges:     (%d %d %d)\n", edges[0], edges[1], edges[2]);


    if (showEdgePoints){
        
        for (i=0; i<3; i++){
            printf("edgePoints %d:\n", i);
            for (j=0; j<edgePoints[i].size(); j++){
                printf("%d:   -- ", edgePoints[i][j]);
                this->nodes[edgePoints[i][j]].print();
            }
        }
        
        printf("\n");
    }

    if (showNodes){
        for (int cN=0; cN<this->nodes.size(); cN++){
            printf("%d  ", cN);
            this->nodes[cN].print(showParamEdges);
//             if (showParamEdges) 
//                 for (i=0; i<nodes[cN].degree(); i++)
//                     printf("    %d\n", (int)nodes[cN].neighbors(i));
        }
    }

    printf("---------------------------------------------------------\n\n");
#endif
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
                const McEdge& cE = this->edges(this->triangles(i).edges[j]);
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
        const McEdge& cE = this->edges(i);

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

template <class ctype>
void PlaneParam<ctype>::checkConsistency(const char* where) const
{
#ifndef NDEBUG
    for (size_t i=0; i<nodes.size(); i++) {

        const Node<ctype>& cN = nodes[i];
        if (cN.isInvalid())
            continue;

        if (isnan(cN.domainPos()[0]) || isnan(cN.domainPos()[1])) {
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
            for (size_t j=0; j<i; j++)
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

template <class ctype>
void DomainTriangle<ctype>::checkConsistency(const char* where) const
{
#ifndef NDEBUG
    if (this->nodes.size()<3){
        printf(where);
        //print(true, true, true);
        assert(this->nodes.size()>=3);
    }

    int i,j;

    // triangles should never contain obsolete nodes
    for (i=0; i<this->nodes.size(); i++)
        if (this->nodes[i].isInvalid()) {
            printf(where);
            printf("***** triangle contains invalid node *****\n");
            assert(false);
        }
        
    PlaneParam<ctype>::checkConsistency(where);

    // check whether all corner nodes are of type CORNER_NODE
    for (i=0; i<3; i++){

        assert(edgePoints[i].size()>=2);

        for (j=0; j<edgePoints[i].size(); j++)
            if (edgePoints[i][j]<0 || edgePoints[i][j]>=this->nodes.size()) {
                printf(where);
                printf("\n***** illegal node index %d in edgePoints array *****\n", edgePoints[i][j]);
                print();
                assert(false);
            }
        
        if (!this->nodes[edgePoints[i][0]].isCORNER_NODE()){
            printf(where);
            printf("***** corner node is not CORNER_NODE *****\n");
            assert(false);
        }

        for (j=1; j<edgePoints[i].size()-1; j++){
            if (this->nodes[edgePoints[i][j]].isCORNER_NODE()){
                printf(where);
                printf("******* corner node found in edgePoint array *** %d *****\n", j);
                assert(false);
            }
            if (this->nodes[edgePoints[i][j]].isINTERIOR_NODE()){
                printf(where);
                printf("******* interior node found in edgePoint array ********\n");
                printf("***     The node has the node number %ld      ***\n", this->nodes[edgePoints[i][j]].getNodeNumber());
                assert(false);
            }
        }
        // check if two subsequent TOUCHING_NODES are connected by an edge
        for (j=0; j<edgePoints[i].size()-1; j++){
            const Node<ctype>& nA = this->nodes[edgePoints[i][j]];
            const Node<ctype>& nB = this->nodes[edgePoints[i][j+1]];

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
            const Node<ctype>& nA = this->nodes[edgePoints[i][j]];
            const Node<ctype>& nB = this->nodes[edgePoints[i][j+2]];

            if (!nA.isInvalid() && !nB.isInvalid() &&
                (nA.isConnectedTo(edgePoints[i][j+2]) || nB.isConnectedTo(edgePoints[i][j]))) {

                printf(where);
                printf("Edge %d,  index %d\n", i, j);
                nA.print();
                nB.print();
                printf("****** two nonadjacent nodes are connected!! *******\n");
                assert(false);
            }
        }
    }

    // check whether all intersection nodes are pointed to from an edgePoint array    

    for (int k=0; k<this->nodes.size(); k++) {

        const Node<ctype>& cN = this->nodes[k];

        if (cN.isINTERSECTION_NODE()){

            bool isIn = false;
            
            for (i=0; i<3; i++)
                for (j=0; j<edgePoints[i].size(); j++)
                    if (edgePoints[i][j]==k)
                        isIn = true;
            
            if (!isIn){
                printf(where);
                printf("***** INTERSECTION NODE not in edgePoints array *****\n");
                cN.print();
                assert(false);
            }
        }
    }

#endif
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
//          display(counter++);
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


// ///////////////////////////////////////////////////////
//   explicit template instantiations
// ///////////////////////////////////////////////////////

template class DomainTriangle<float>;
template class DomainTriangle<double>;

template class PlaneParam<float>;
template class PlaneParam<double>;

template class Node<float>;
template class Node<double>;

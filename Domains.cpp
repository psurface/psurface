#include "config.h"

#include "Domains.h"

using namespace psurface;

template <class ctype>
void DomainTriangle<ctype>::insertExtraEdges()
{
    int i,j;

    // add the missing paramEdges on the base grid triangle edges
    for (i=0; i<3; i++){

        for (j=1; j<edgePoints[i].size(); j++){

            if ((this->nodes[edgePoints[i][j]].isINTERSECTION_NODE() ||
                 this->nodes[edgePoints[i][j]].isGHOST_NODE() ||
                 this->nodes[edgePoints[i][j-1]].isINTERSECTION_NODE() ||
                 this->nodes[edgePoints[i][j-1]].isGHOST_NODE()) &&
                !this->nodes[edgePoints[i][j]].isConnectedTo(edgePoints[i][j-1])){


                this->addEdge(edgePoints[i][j-1], edgePoints[i][j], true);
            }
        }

    }

    // turn quadrangles into two triangles
    for (i=0; i<3; i++){

        for (j=1; j<edgePoints[i].size(); j++)

            if (this->nodes[edgePoints[i][j]].isINTERSECTION_NODE()){

                //int interiorPoint = nodes[edgePoints[i][j]].neighbors(0);
                int interiorPoint = this->nodes[edgePoints[i][j]].theInteriorNode();


                if (!this->nodes[interiorPoint].isConnectedTo(edgePoints[i][j-1]))
                    this->addEdge(edgePoints[i][j-1], interiorPoint, true);

            }

    }

}



template <class ctype>
void DomainTriangle<ctype>::flip()
{
    // flip points
    std::swap(vertices[1], vertices[2]);

    // flip edges
    std::swap(edges[0], edges[2]);

    // flip edgePoints array
    std::swap(edgePoints[0], edgePoints[2]);

    for (int i=0; i<3; i++)
        std::reverse(edgePoints[i].begin(), edgePoints[i].end());

    // Change the pointers of intersection nodes to their respective positions
    // in the edgePoints arrays.  This is just in case that the pointLocationStructure
    // is intact.
    int i,j;
    for (i=0; i<3; i++)
        for (j=1; j<edgePoints[i].size()-1; j++){
            if (this->nodes[edgePoints[i][j]].isINTERSECTION_NODE()){
                this->nodes[edgePoints[i][j]].setDomainEdge(i);
                this->nodes[edgePoints[i][j]].setDomainEdgePosition(j);
            }
        }

    // turn the parametrization
    /** \todo This is slow and should be reprogrammed! */
    this->installWorldCoordinates(StaticVector<ctype,2>(0,0), StaticVector<ctype,2>(1,0), StaticVector<ctype,2>(0,1));
    PlaneParam<ctype>::installBarycentricCoordinates(StaticVector<ctype,2>(0,0), StaticVector<ctype,2>(0,1), StaticVector<ctype,2>(1,0));
}

template <class ctype>
void DomainTriangle<ctype>::rotate()
{
    std::rotate(vertices.begin(),   vertices.end()-1,   vertices.end());
    std::rotate(edges.begin(),      edges.end()-1,      edges.end());
    std::rotate(edgePoints.begin(), edgePoints.end()-1, edgePoints.end());

    // turn the parametrization
    /// \todo This is slow and should be replaced!
    this->installWorldCoordinates(StaticVector<ctype,2>(0,0), StaticVector<ctype,2>(1,0), StaticVector<ctype,2>(0,1));
    PlaneParam<ctype>::installBarycentricCoordinates(StaticVector<ctype,2>(0,1), StaticVector<ctype,2>(0,0), StaticVector<ctype,2>(1,0));
}


template <class ctype>
void DomainTriangle<ctype>::updateEdgePoints(int oldNode, int newNode)
{
    int i;
    for (i=0; i<3; i++){
        if (edgePoints[i][0]==oldNode)
            edgePoints[i][0] = newNode;
        if (edgePoints[i].back()==oldNode)
            edgePoints[i].back() = newNode;
    }
}


//////////////////////////////////////////////////////////////
// TOUCHING_NODES tend to 'drift away' from their original position
// due to repeated switching between barycentric and world coordinates
// we pull them back in place.  It'd be nice if there was a more elegant way
template <class ctype>
void DomainTriangle<ctype>::adjustTouchingNodes()
{
    int i;

    for (i=1; i<edgePoints[0].size()-1; i++)
        if (this->nodes[edgePoints[0][i]].isTOUCHING_NODE()
            || this->nodes[edgePoints[0][i]].isINTERSECTION_NODE()){
            StaticVector<ctype,2> tmp = this->nodes[edgePoints[0][i]].domainPos();
            ctype diff = (1.0f - tmp[0] - tmp[1]);
            tmp[0] += 0.5*diff;
            tmp[1] += 0.5*diff;
            this->nodes[edgePoints[0][i]].setDomainPos(tmp);
        }

    for (i=1; i<edgePoints[1].size()-1; i++)
        if (this->nodes[edgePoints[1][i]].isTOUCHING_NODE() || this->nodes[edgePoints[1][i]].isINTERSECTION_NODE())
            this->nodes[edgePoints[1][i]].setDomainPos(StaticVector<ctype,2>(0, this->nodes[edgePoints[1][i]].domainPos()[1]));

    for (i=1; i<edgePoints[2].size()-1; i++)
        if (this->nodes[edgePoints[2][i]].isTOUCHING_NODE() || this->nodes[edgePoints[2][i]].isINTERSECTION_NODE())
            this->nodes[edgePoints[2][i]].setDomainPos(StaticVector<ctype,2>(this->nodes[edgePoints[2][i]].domainPos()[0], 0));
}


template <class ctype>
void DomainTriangle<ctype>::createPointLocationStructure()
{
    //print(true, true, true);
    checkConsistency("BeforeCreate");

    for (int i=0; i<this->nodes.size(); i++)
        if (this->nodes[i].isINTERIOR_NODE())
            this->makeCyclicInteriorNode(this->nodes[i]);

    checkConsistency("AfterInterior");

    for (int i=0; i<3; i++) {

        this->makeCyclicBoundaryNode(this->nodes[edgePoints[i][0]],
                               edgePoints[i][1],
                               edgePoints[(i+2)%3][edgePoints[(i+2)%3].size()-2]);

        // should be setCorner()
        this->nodes[edgePoints[i][0]].setDomainEdge(i);

        for (int j=1; j<edgePoints[i].size()-1; j++){
            this->makeCyclicBoundaryNode(this->nodes[edgePoints[i][j]], edgePoints[i][j+1], edgePoints[i][j-1]);

            if (this->nodes[edgePoints[i][j]].isINTERSECTION_NODE() || this->nodes[edgePoints[i][j]].isTOUCHING_NODE()){
                this->nodes[edgePoints[i][j]].setDomainEdge(i);
                this->nodes[edgePoints[i][j]].setDomainEdgePosition(j);
            }
        }

        checkConsistency("AfterEdges");

    }

}

template <class ctype>
void DomainTriangle<ctype>::print(bool showEdgePoints, bool showParamEdges, bool showNodes) const
{
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
                std::cout << "******* interior node found in edgePoint array ********" << std::endl;
                std::cout << "***     The node has the node number "
                          << this->nodes[edgePoints[i][j]].getNodeNumber() << "      ***" << std::endl;
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

// ///////////////////////////////////////////////////////
//   explicit template instantiations
// ///////////////////////////////////////////////////////

template class PSURFACE_EXPORT DomainTriangle<float>;
template class PSURFACE_EXPORT DomainTriangle<double>;

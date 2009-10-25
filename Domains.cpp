#include <psurface/Domains.h>



void DomainTriangle::insertExtraEdges()
{
    int i,j;
    
    // add the missing paramEdges on the base grid triangle edges
    for (i=0; i<3; i++){
        
        for (j=1; j<edgePoints[i].size(); j++){
            
            if ((nodes[edgePoints[i][j]].isINTERSECTION_NODE() ||
                nodes[edgePoints[i][j]].isGHOST_NODE() ||
                nodes[edgePoints[i][j-1]].isINTERSECTION_NODE() ||
                nodes[edgePoints[i][j-1]].isGHOST_NODE()) &&
                !nodes[edgePoints[i][j]].isConnectedTo(edgePoints[i][j-1])){

                
                addEdge(edgePoints[i][j-1], edgePoints[i][j], true);
            }
        }

    }

    // turn quadrangles into two triangles
    for (i=0; i<3; i++){
    
        for (j=1; j<edgePoints[i].size(); j++)
            
            if (nodes[edgePoints[i][j]].isINTERSECTION_NODE()){

                //int interiorPoint = nodes[edgePoints[i][j]].neighbors(0);
                int interiorPoint = nodes[edgePoints[i][j]].theInteriorNode();


                if (!nodes[interiorPoint].isConnectedTo(edgePoints[i][j-1])){
                    //printf("i: %d   j: %d\n", i, j);
                    addEdge(edgePoints[i][j-1], interiorPoint, true);
                }
            }  

    }

}



void DomainTriangle::flip()
{
    // flip points
    //vertices.swap(1,2);
    int tmp = vertices[1];
    vertices[1] = vertices[2];
    vertices[2] = tmp;

    // flip edges
    //edges.swap(0,2);
    tmp = edges[0];
    edges[0] = edges[2];
    edges[2] = tmp;

    // flip edgePoints array
    //edgePoints.swap(0,2);
    std::vector<NodeIdx> tmp1 = edgePoints[0];
    edgePoints[0] = edgePoints[2];
    edgePoints[2] = tmp1;

    for (int i=0; i<3; i++)
        std::reverse(edgePoints[i].begin(), edgePoints[i].end());
    
    // Change the pointers of intersection nodes to their respective positions
    // in the edgePoints arrays.  This is just in case that the pointLocationStructure
    // is intact.
    int i,j;
    for (i=0; i<3; i++)
        for (j=1; j<edgePoints[i].size()-1; j++){
            if (nodes[edgePoints[i][j]].isINTERSECTION_NODE()){
                nodes[edgePoints[i][j]].setDomainEdge(i);
                nodes[edgePoints[i][j]].setDomainEdgePosition(j);
            }
        }

    // turn the parametrization
    /** \todo This is slow and should be reprogrammed! */
    installWorldCoordinates(StaticVector<float,2>(0,0), StaticVector<float,2>(1,0), StaticVector<float,2>(0,1));
    PlaneParam<float>::installBarycentricCoordinates(StaticVector<float,2>(0,0), StaticVector<float,2>(0,1), StaticVector<float,2>(1,0));
}

void DomainTriangle::rotate()
{
    rotate(vertices,1);
    rotate(edges,1);
    rotate(edgePoints,1);
    
    // turn the parametrization
    /// \todo This is slow and should be replaced!
    installWorldCoordinates(StaticVector<float,2>(0,0), StaticVector<float,2>(1,0), StaticVector<float,2>(0,1));
    PlaneParam<float>::installBarycentricCoordinates(StaticVector<float,2>(0,1), StaticVector<float,2>(0,0), StaticVector<float,2>(1,0));
}





void DomainTriangle::updateEdgePoints(int oldNode, int newNode)
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
void DomainTriangle::adjustTouchingNodes()
{
    int i;

    for (i=1; i<edgePoints[0].size()-1; i++)
        if (nodes[edgePoints[0][i]].isTOUCHING_NODE() || nodes[edgePoints[0][i]].isINTERSECTION_NODE()){
            StaticVector<float,2> tmp = nodes[edgePoints[0][i]].domainPos();
            float diff = (1.0f - tmp[0] - tmp[1]);
            tmp[0] += 0.5*diff;
            tmp[1] += 0.5*diff;
            nodes[edgePoints[0][i]].setDomainPos(tmp);
        }

    for (i=1; i<edgePoints[1].size()-1; i++)
        if (nodes[edgePoints[1][i]].isTOUCHING_NODE() || nodes[edgePoints[1][i]].isINTERSECTION_NODE())
            nodes[edgePoints[1][i]].setDomainPos(StaticVector<float,2>(0, nodes[edgePoints[1][i]].domainPos()[1]));
    
    for (i=1; i<edgePoints[2].size()-1; i++)
        if (nodes[edgePoints[2][i]].isTOUCHING_NODE() || nodes[edgePoints[2][i]].isINTERSECTION_NODE())
            nodes[edgePoints[2][i]].setDomainPos(StaticVector<float,2>(nodes[edgePoints[2][i]].domainPos()[0], 0));
}
        



void DomainTriangle::createPointLocationStructure() 
{
    int i, j;
    //print(true, true, true);
    checkConsistency("BeforeCreate");

    for (i=0; i<nodes.size(); i++){
        
        if (nodes[i].isINTERIOR_NODE()){
            makeCyclicInteriorNode(nodes[i]);
        }
    }

    checkConsistency("AfterInterior");

    for (i=0; i<3; i++) {
        
        makeCyclicBoundaryNode(nodes[edgePoints[i][0]], 
                               edgePoints[i][1], 
                               edgePoints[(i+2)%3][edgePoints[(i+2)%3].size()-2]);

        // should be setCorner()
        nodes[edgePoints[i][0]].setDomainEdge(i);

        for (j=1; j<edgePoints[i].size()-1; j++){
            makeCyclicBoundaryNode(nodes[edgePoints[i][j]], edgePoints[i][j+1], edgePoints[i][j-1]);

            if (nodes[edgePoints[i][j]].isINTERSECTION_NODE() || nodes[edgePoints[i][j]].isTOUCHING_NODE()){
                nodes[edgePoints[i][j]].setDomainEdge(i);
                nodes[edgePoints[i][j]].setDomainEdgePosition(j);
            }
        }

        checkConsistency("AfterEdges");

    }

}


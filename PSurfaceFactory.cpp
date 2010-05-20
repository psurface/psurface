#include <vector>

#include <psurface/PSurface.h>
#include <psurface/PSurfaceFactory.h>

template <int dim, class ctype>
void PSurfaceFactory<dim,ctype>::setTargetSurface(Surface* surface)
{
    psurface_->surface = surface;
}

template <int dim, class ctype>
void PSurfaceFactory<dim,ctype>::insertVertex(const StaticVector<ctype,dim+1>& position)
{
    psurface_->newVertex(position);
}

template <int dim, class ctype>
unsigned int  PSurfaceFactory<dim,ctype>::insertSimplex(const std::tr1::array<unsigned int, dim+1>& v)
{
    unsigned int idx = psurface_->createSpaceForTriangle(v[0], v[1], v[2]);
    psurface_->integrateTriangle(idx);
    return idx;
}

template <int dim, class ctype>
void PSurfaceFactory<dim,ctype>::insertTargetVertexMapping(unsigned int targetVertex, 
                                                           unsigned int domainTriangle, 
                                                           const StaticVector<ctype,dim>& domainLocalPosition,
                                                           NodeBundle& projectedTo,
                                                           int& domainVertex)
{
    const double eps = 0.0001;

    // determine which type of node to add
    typename Node<ctype>::NodeType newType = Node<ctype>::INTERIOR_NODE;
    int dir = -1;
    double mu;
    
    // if the normal projection hits a base grid vertex, this is the vertex
    domainVertex = -1;
    
    if (domainLocalPosition[0] < eps) {
        dir = 1;
        mu = 1-domainLocalPosition[1];
        if (domainLocalPosition[1] < eps) {
            newType = Node<ctype>::CORNER_NODE;
            domainVertex = psurface_->triangles(domainTriangle).vertices[2];
        } else if (domainLocalPosition[1] > 1-eps) {
            newType = Node<ctype>::CORNER_NODE;
            domainVertex = psurface_->triangles(domainTriangle).vertices[1];
        } else {
            newType = Node<ctype>::TOUCHING_NODE;
        }
    } else if (domainLocalPosition[1] < eps) {
        dir = 2;
        mu = domainLocalPosition[0];
        if (domainLocalPosition[0] < eps) {
            newType = Node<ctype>::CORNER_NODE;
            domainVertex = psurface_->triangles(domainTriangle).vertices[2];
        } else if (domainLocalPosition[0] > 1-eps) {
            newType = Node<ctype>::CORNER_NODE;
            domainVertex = psurface_->triangles(domainTriangle).vertices[0];
        } else {
            newType = Node<ctype>::TOUCHING_NODE;
        }
    } else if (1-domainLocalPosition[0]-domainLocalPosition[1] < eps) {
        dir = 0;
        mu = 1-domainLocalPosition[0];
        if (domainLocalPosition[1] < eps) {
            newType = Node<ctype>::CORNER_NODE;
            domainVertex = psurface_->triangles(domainTriangle).vertices[0];
        } else if (domainLocalPosition[1] > 1-eps) {
            newType = Node<ctype>::CORNER_NODE;
            domainVertex = psurface_->triangles(domainTriangle).vertices[1];
        } else {
            newType = Node<ctype>::TOUCHING_NODE;
        }
    }
    
    if (newType==Node<ctype>::TOUCHING_NODE) {
        
        // find the other triangle, if there is one
        int neighboringTri = psurface_->getNeighboringTriangle(domainTriangle, dir);
        
        if (neighboringTri == -1) {
            NodeIdx newNodeNumber = addTouchingNode(domainTriangle, domainLocalPosition, dir, targetVertex);
            projectedTo.resize(1);
            projectedTo[0].setValue(domainTriangle, newNodeNumber);
        } else {
            // find domain pos on other triangle
            int commonEdge = psurface_->triangles(domainTriangle).getCommonEdgeWith(psurface_->triangles(neighboringTri));
            int dir2 = psurface_->triangles(neighboringTri).getEdge(commonEdge);
            StaticVector<ctype,2> dP2((dir2==0)*(mu) + (dir2==2)*(1-mu), (dir2==0)*(1-mu) + (dir2==1)*(mu));
            
            // insert touching node pair
            NodeIdx newNodeNumber = addTouchingNodePair(domainTriangle, neighboringTri, domainLocalPosition, dP2, 
                                                                   dir, dir2, targetVertex);
            projectedTo.resize(2);
            projectedTo[0].setValue(domainTriangle, newNodeNumber);
            projectedTo[1].setValue(neighboringTri, psurface_->triangles(neighboringTri).nodes.size()-1);
        }
        
    } else if (newType==Node<ctype>::CORNER_NODE) {
        
        addCornerNodeBundle(domainVertex, targetVertex);
        std::vector<int> neighboringTris = psurface_->getTrianglesPerVertex(domainVertex);
        projectedTo.resize(neighboringTris.size());
        for (size_t j=0; j<neighboringTris.size(); j++) {
            projectedTo[j].setValue(neighboringTris[j],
                                                  psurface_->triangles(neighboringTris[j]).nodes.size()-1);
        }
        
    } else {
        
        NodeIdx newNodeNumber = addInteriorNode(domainTriangle, domainLocalPosition, targetVertex);
        
        projectedTo.resize(1);
        projectedTo[0].setValue(domainTriangle, newNodeNumber);
        
    }

}

template <int dim, class ctype>
void PSurfaceFactory<dim,ctype>::insertGhostNode(unsigned int domainVertex,
                                                 unsigned int targetTriangle,
                                                 const StaticVector<ctype,dim>& targetLocalPosition)
{
    std::vector<int> neighbors = psurface_->getTrianglesPerVertex(domainVertex);
    
    for (size_t i=0; i<neighbors.size(); i++) {
        
        int corner = psurface_->triangles(neighbors[i]).getCorner(domainVertex);
        addGhostNode(neighbors[i], corner, targetTriangle, targetLocalPosition);
        
    }
}

template <int dim, class ctype>
void PSurfaceFactory<dim,ctype>::insertEdge()
{
}

template <int dim, class ctype>
void PSurfaceFactory<dim,ctype>::addCornerNodeBundle(int domainVertex, int targetVertex)
{
    std::vector<int> neighbors = psurface_->getTrianglesPerVertex(domainVertex);

    for (size_t i=0; i<neighbors.size(); i++) {

        int corner = psurface_->triangles(neighbors[i]).getCorner(domainVertex);
        addCornerNode(neighbors[i], corner, targetVertex);

    }
        
}


template <int dim, class ctype>
NodeIdx PSurfaceFactory<dim,ctype>::addInteriorNode(int tri, const StaticVector<ctype,2>& dom, int nodeNumber)
{
    psurface_->triangles(tri).nodes.push_back(Node<ctype>(dom, nodeNumber, Node<ctype>::INTERIOR_NODE));
    return psurface_->triangles(tri).nodes.size()-1;
}

template <int dim, class ctype>
NodeIdx PSurfaceFactory<dim,ctype>::addGhostNode(int tri, int corner, int targetTri, const StaticVector<ctype,2>& localTargetCoords)
{
    psurface_->triangles(tri).nodes.push_back(Node<ctype>());
    psurface_->triangles(tri).nodes.back().makeGhostNode(corner, targetTri, localTargetCoords);
    return psurface_->triangles(tri).nodes.size()-1;
}

template <int dim, class ctype>
NodeIdx PSurfaceFactory<dim,ctype>::addCornerNode(int tri, int corner, int nodeNumber)
{
    DomainTriangle<ctype>& cT = psurface_->triangles(tri);

    cT.nodes.push_back(Node<ctype>());
    cT.nodes.back().makeCornerNode(corner, nodeNumber);
    return cT.nodes.size()-1;
}

// BUG: The node needs to be entered in the edgepoint arrays
template <int dim, class ctype>
NodeBundle PSurfaceFactory<dim,ctype>::addIntersectionNodePair(int tri1, int tri2,
                                                const StaticVector<ctype,2>& dP1, const StaticVector<ctype,2>& dP2, 
                                                int edge1, int edge2, const StaticVector<ctype,3>& range)
{
    // We return the pair of new nodes
    NodeBundle result(2);
    result[0].tri = tri1;
    result[1].tri = tri2;

    DomainTriangle<ctype>& cT1 = psurface_->triangles(tri1);
    DomainTriangle<ctype>& cT2 = psurface_->triangles(tri2);

    psurface_->iPos.push_back(range);
    int nodeNumber = psurface_->iPos.size()-1;

    cT1.nodes.push_back(Node<ctype>());
    cT2.nodes.push_back(Node<ctype>());

    result[0].idx = cT1.nodes.size()-1;
    result[1].idx = cT2.nodes.size()-1;
    
    cT1.nodes.back().setValue(dP1, nodeNumber, Node<ctype>::INTERSECTION_NODE);
    cT2.nodes.back().setValue(dP2, nodeNumber, Node<ctype>::INTERSECTION_NODE);

    cT1.nodes.back().setDomainEdge(edge1);
    cT2.nodes.back().setDomainEdge(edge2);

    return result;
}

// BUG: The node needs to be entered in the edgepoint arrays
template <int dim, class ctype>
NodeIdx PSurfaceFactory<dim,ctype>::addTouchingNode(int tri, const StaticVector<ctype,2>& dP, int edge, int nodeNumber)
{
    DomainTriangle<ctype>& cT = psurface_->triangles(tri);

    cT.nodes.push_back(Node<ctype>());
    
    cT.nodes.back().setValue(dP, nodeNumber, Node<ctype>::TOUCHING_NODE);
    cT.nodes.back().setDomainEdge(edge);
    return cT.nodes.size()-1;
}

// BUG: The node needs to be entered in the edgepoint arrays
template <int dim, class ctype>
NodeIdx PSurfaceFactory<dim,ctype>::addTouchingNodePair(int tri1, int tri2,
                                            const StaticVector<ctype,2>& dP1, const StaticVector<ctype,2>& dP2, 
                                            int edge1, int edge2, int nodeNumber)
{
    DomainTriangle<ctype>& cT1 = psurface_->triangles(tri1);
    DomainTriangle<ctype>& cT2 = psurface_->triangles(tri2);

    cT1.nodes.push_back(Node<ctype>());
    cT2.nodes.push_back(Node<ctype>());
    
    cT1.nodes.back().setValue(dP1, nodeNumber, Node<ctype>::TOUCHING_NODE);
    cT2.nodes.back().setValue(dP2, nodeNumber, Node<ctype>::TOUCHING_NODE);

    cT1.nodes.back().setDomainEdge(edge1);
    cT2.nodes.back().setDomainEdge(edge2);

    return cT1.nodes.size()-1;
}

template <int dim, class ctype>
void PSurfaceFactory<dim,ctype>::addParTriangle(int tri, const std::tr1::array<int,3>& p)
{
    DomainTriangle<ctype>& cT = psurface_->triangles(tri);

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


// ///////////////////////////////////////////////////////////////////////
//   Explicitly instantiate 'float' and 'double' versions of this code
// ///////////////////////////////////////////////////////////////////////

//template class PSurfaceFactory<1,float>;
//template class PSurfaceFactory<1,double>;

template class PSurfaceFactory<2,float>;
template class PSurfaceFactory<2,double>;

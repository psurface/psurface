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
            NodeIdx newNodeNumber = psurface_->addTouchingNode(domainTriangle, domainLocalPosition, dir, targetVertex);
            projectedTo.resize(1);
            projectedTo[0].setValue(domainTriangle, newNodeNumber);
        } else {
            // find domain pos on other triangle
            int commonEdge = psurface_->triangles(domainTriangle).getCommonEdgeWith(psurface_->triangles(neighboringTri));
            int dir2 = psurface_->triangles(neighboringTri).getEdge(commonEdge);
            StaticVector<ctype,2> dP2((dir2==0)*(mu) + (dir2==2)*(1-mu), (dir2==0)*(1-mu) + (dir2==1)*(mu));
            
            // insert touching node pair
            NodeIdx newNodeNumber = psurface_->addTouchingNodePair(domainTriangle, neighboringTri, domainLocalPosition, dP2, 
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
        
        NodeIdx newNodeNumber = psurface_->addInteriorNode(domainTriangle, domainLocalPosition, targetVertex);
        
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
    
    for (int i=0; i<neighbors.size(); i++) {
        
        int corner = psurface_->triangles(neighbors[i]).getCorner(domainVertex);
        psurface_->addGhostNode(neighbors[i], corner, targetTriangle, targetLocalPosition);
        
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

    for (int i=0; i<neighbors.size(); i++) {

        int corner = psurface_->triangles(neighbors[i]).getCorner(domainVertex);
        psurface_->addCornerNode(neighbors[i], corner, targetVertex);

    }
        
}

// ///////////////////////////////////////////////////////////////////////
//   Explicitly instantiate 'float' and 'double' versions of this code
// ///////////////////////////////////////////////////////////////////////

//template class PSurfaceFactory<1,float>;
//template class PSurfaceFactory<1,double>;

template class PSurfaceFactory<2,float>;
template class PSurfaceFactory<2,double>;

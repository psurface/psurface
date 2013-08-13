#include <cstddef>
#include "NormalProjector.h"
#include "PSurfaceFactory.h"
#include "DirectionFunction.h"

#include "StaticVector.h"
#include "StaticMatrix.h"

#include "NodeBundle.h"
#include "PathVertex.h"

#ifdef PSURFACE_STANDALONE
#include "TargetSurface.h"
#else
#include "hxsurface/Surface.h"
#endif

#include <stdexcept>
#include <exception>
#include <vector>
#include <set>

using namespace psurface;

/** \brief An Exception that is thrown whenever during the insertion of a target edge
 *         the inverse projection is not unique or does not exist at all.
 */
class WrongEdgeException: public std::exception
{
public:
    WrongEdgeException(std::string str) :
        str_(str)
    {}

    virtual const char* what() const throw()
    {
        return str_.c_str();
    }

    virtual ~WrongEdgeException() throw()
    {}

private:
    const std::string str_;
};

template <class ctype>
void NormalProjector<ctype>::project(const Surface* targetSurface,
                                     const DirectionFunction<3,ctype>* domainDirection,
                                     const DirectionFunction<3,ctype>* targetDirection)
{
    const double eps = 1e-4;

    PSurfaceFactory<2,ctype> factory(psurface_);

    const Surface* surf = psurface_->surface;

    // //////////////////////////
    // compute normals
    // //////////////////////////

    int nPoints = psurface_->getNumVertices();
        
    std::vector<StaticVector<ctype,3> > domainNormals(nPoints);
    computeDiscreteDomainDirections(domainDirection, domainNormals);

    // /////////////////////////////////////////////////////////////
    //   Compute the vertex normals of the target side
    // /////////////////////////////////////////////////////////////


    std::vector<StaticVector<ctype,3> > targetNormals(targetSurface->points.size());
    computeDiscreteTargetDirections(targetSurface, targetDirection, targetNormals);

    // /////////////////////////////////////////////////////////////////////////////////////
    // Insert the vertices of the contact boundary as nodes on the intermediate manifold
    // /////////////////////////////////////////////////////////////////////////////////////

    // TODO: Use an octree here.  The MultiDimOctree should get a method that returns
    // all verticess close to the cone spanned by the search directions at the domain
    // surface triangles.  Only those vertices should then be checked.
    // We need this quickly, because it is mentioned in the psurface papers, which
    // has already been submitted.

    // This array stores the preimages of each vertex in the target surface
    std::vector<NodeBundle> projectedTo(surf->points.size());

    // This bitfield marks whether base grid vertices already have a
    // corresponding image
    std::vector<bool> vertexHasBeenHandled(psurface_->getNumVertices(), false);

    // Loop over the vertices of the target surface
    for (size_t i=0; i<targetSurface->points.size(); i++) {

        StaticVector<ctype,2> bestDPos;
        int bestTri = -1;
        ctype bestDist = std::numeric_limits<ctype>::max();

        // magic to use a McVec3f as the argument
        StaticVector<ctype,3> targetVertex;
        for (int k=0; k<3; k++)
            targetVertex[k] = surf->points[i][k];

        for (size_t j=0; j<psurface_->getNumTriangles(); j++) {
            const StaticVector<ctype,3>& p0 = psurface_->vertices(psurface_->triangles(j).vertices[0]);
            const StaticVector<ctype,3>& p1 = psurface_->vertices(psurface_->triangles(j).vertices[1]);
            const StaticVector<ctype,3>& p2 = psurface_->vertices(psurface_->triangles(j).vertices[2]);

            const StaticVector<ctype,3>& n0 = domainNormals[psurface_->triangles(j).vertices[0]];
            const StaticVector<ctype,3>& n1 = domainNormals[psurface_->triangles(j).vertices[1]];
            const StaticVector<ctype,3>& n2 = domainNormals[psurface_->triangles(j).vertices[2]];

            StaticVector<ctype,3> x; // the unknown...
            if (computeInverseNormalProjection(p0, p1, p2, n0, n1, n2, targetVertex, x)) {

                // We want that the line from the domain surface to its projection
                // approaches the target surface from the front side, i.e., it should
                // not pass through the body represented by the target surface.
                // We do a simplified test by comparing the connecting segment
                // with the normal at the target surface and the normal at the
                // domain surface
                StaticVector<ctype,3> base       = p0*x[0] + p1*x[1] + (1-x[0]-x[1])*p2;
                StaticVector<ctype,3> baseNormal = n0*x[0] + n1*x[1] + (1-x[0]-x[1])*n2;
                StaticVector<ctype,3> segment(surf->points[i][0] - base[0],
                                surf->points[i][1] - base[1],
                                surf->points[i][2] - base[2]);
                
                ctype distance = segment.length() * segment.length();

                // if both conditions are not fulfilled we might want to allow some overlaps
                if(segment.dot(targetNormals[i]) > eps
                    && segment.dot(baseNormal) < -eps) {
                        if (distance > 0.01) // TODO this value should be set problem dependent
                            continue;
                } else if( segment.dot(targetNormals[i]) > eps
                    || segment.dot(baseNormal) < -eps)
                    continue;

                // There may be several inverse orthogonal projections.
                // We want the shortest one.

                if (distance < bestDist) {

                    bestDist = distance;
                    bestDPos[0] = x[0];
                    bestDPos[1] = x[1];
                    bestTri  = j;

                }
                
            }

        }
        
        if (bestTri != -1) {

            int domainVertex;
            factory.insertTargetVertexMapping(i, bestTri, bestDPos, projectedTo[i], domainVertex);
            if (domainVertex >= 0)
                vertexHasBeenHandled[domainVertex] = true;

        }
        
    }

    // ///////////////////////////////////////////////////////////////////
    //   Place ghost nodes at the vertices of the domain surface
    // ///////////////////////////////////////////////////////////////////
    
    // TODO: Use an octree here.  The MultiDimOctree should get a method that returns
    // all triangles close to a given ray.  Only those triangles should then be checked.
    // We need this quickly, because it is mentioned in the psurface papers, which
    // has already been submitted.

    for (int i=0; i<psurface_->getNumVertices(); i++) {

        // Has the vertex been hit by the projection of a target vertex already?
        if (vertexHasBeenHandled[i])
            continue;

        StaticVector<ctype,2> bestDPos;
        int bestTri = -1;
        ctype bestDist = std::numeric_limits<ctype>::max();

        const StaticVector<ctype,3>& basePoint = psurface_->vertices(i);
        StaticVector<ctype,3> normal;
        normal[0] = domainNormals[i][0];
        normal[1] = domainNormals[i][1];
        normal[2] = domainNormals[i][2];

        for (int j=0; j<targetSurface->triangles.size(); j++) {

            StaticVector<ctype,2> domainPos;
            ctype dist;

            // copy the coordinates, because they are stored in a McVec3f when compiled as part of Amira
            StaticVector<ctype,3> p0, p1, p2;
            for (int k=0; k<3; k++) {
                p0[k] = surf->points[targetSurface->triangles[j].points[0]][k];
                p1[k] = surf->points[targetSurface->triangles[j].points[1]][k];
                p2[k] = surf->points[targetSurface->triangles[j].points[2]][k];
            }

            if (rayIntersectsTriangle(basePoint, normal, p0, p1, p2, domainPos, dist, eps)) {

                if (dist<bestDist) {
                    bestTri = j;
                    bestDPos = domainPos;
                    bestDist = dist;
                }

            }
            
        }

        // Set ghost node mapping to the closest triangle intersected by the normal ray
        if (bestTri != -1)
            factory.insertGhostNode(i, bestTri, bestDPos);

    }

    // ////////////////////////////////////////////////////////////
    // Insert the edges.
    // We loop over all sides of all triangles.  For each side we
    // remember whether we have seen it before.
    // ////////////////////////////////////////////////////////////

    std::set<std::pair<int,int> > visitedEdges;

    for (int i=0; i<targetSurface->triangles.size(); i++) {

        for (int j=0; j<3; j++) {
            
            int from = targetSurface->triangles[i].points[j];
            int to   = targetSurface->triangles[i].points[(j+1)%3];
            
            // The two numbers max and min uniquely identify 'from' and 'to'.
            // However, we have always min<max, and hence we do not have to
            // worry about orientation anymore
            int max = std::max(from,to);
            int min = from + to - max;
            
            // Mark this edge as visited.  The return value is true, if this
            // edge has not been visited before.
            bool wasInserted = visitedEdges.insert(std::make_pair(min,max)).second;

            if (wasInserted) {

                // store the path so we don't have to compute it twice
                std::vector<PathVertex<ctype> > edgePath(1);
                if (edgeCanBeInserted(domainNormals, from, to, projectedTo, edgePath))
                    insertEdge(factory, from, to, edgePath);
                else {
                    //std::cout << "Skipping edge (" << from << ", " << to << ") ..." << std::endl;
                }
            }

        }

    }

    setupEdgePointArrays();
}


template <class ctype>
void NormalProjector<ctype>::computeDiscreteDomainDirections(const DirectionFunction<3,ctype>* direction,
                                                             std::vector<StaticVector<ctype,3> >& normals)
{
    int nPoints = psurface_->getNumVertices();
    int nTriangles = psurface_->getNumTriangles();
    
    normals.assign(nPoints, StaticVector<ctype,3>(0.0));
    
    if (direction) {
        
        for (int i=0; i<nPoints; i++) {
            
            if (dynamic_cast<const AnalyticDirectionFunction<3,ctype>*>(direction)) {
                normals[i] = (*dynamic_cast<const AnalyticDirectionFunction<3,ctype>*>(direction))(psurface_->vertices(i));
            } else if (dynamic_cast<const DiscreteDirectionFunction<3,ctype>*>(direction))
                normals[i] = (*dynamic_cast<const DiscreteDirectionFunction<3,ctype>*>(direction))(i);
            else {
                std::cerr << "Domain direction function not properly set!" << std::endl;
                abort();
            }
            
        }
        
    } else {
        
        for (int i=0; i<nTriangles; i++) {
            
            int p0 = psurface_->triangles(i).vertices[0];
            int p1 = psurface_->triangles(i).vertices[1];
            int p2 = psurface_->triangles(i).vertices[2];
            
            StaticVector<ctype,3> a = psurface_->vertices(p1) - psurface_->vertices(p0);
            StaticVector<ctype,3> b = psurface_->vertices(p2) - psurface_->vertices(p0);
            StaticVector<ctype,3> triNormal = a.cross(b);
            triNormal.normalize();
            
            normals[p0] += triNormal;
            normals[p1] += triNormal;
            normals[p2] += triNormal;
            
        }

        for (int i=0; i<nPoints; i++)
            normals[i].normalize();

    }

}

template <class ctype>
void NormalProjector<ctype>::computeDiscreteTargetDirections(const Surface* targetSurface,
                                                             const DirectionFunction<3,ctype>* direction,
                                                             std::vector<StaticVector<ctype,3> >& normals)
{
    int nPoints    = targetSurface->points.size();
    int nTriangles = targetSurface->triangles.size();

    normals.assign(nPoints, StaticVector<ctype,3>(0.0));

    if (direction) {
        
        for (int i=0; i<nPoints; i++) {
            
            if (dynamic_cast<const AnalyticDirectionFunction<3,ctype>*>(direction)) {
                StaticVector<ctype,3> p;
                for (int j=0; j<3; j++)
                    p[j] = targetSurface->points[i][j];
                normals[i] = (*dynamic_cast<const AnalyticDirectionFunction<3,ctype>*>(direction))(p);
            } else if (dynamic_cast<const DiscreteDirectionFunction<3,ctype>*>(direction))
                normals[i] = (*dynamic_cast<const DiscreteDirectionFunction<3,ctype>*>(direction))(i);
            else {
                std::cerr << "Target direction function not properly set!" << std::endl;
                abort();
            }
            
        }
        
    } else {
        
        for (int i=0; i<nTriangles; i++) {
            
            int p0 = targetSurface->triangles[i].points[0];
            int p1 = targetSurface->triangles[i].points[1];
            int p2 = targetSurface->triangles[i].points[2];
            
            StaticVector<ctype,3> a, b;
            
            for (int j=0; j<3; j++) {
                a[j] = targetSurface->points[p1][j] - targetSurface->points[p0][j];
                b[j] = targetSurface->points[p2][j] - targetSurface->points[p0][j];
            }        
            
            StaticVector<ctype,3> triNormal = a.cross(b);
            triNormal.normalize();
            
            normals[p0] += triNormal;
            normals[p1] += triNormal;
            normals[p2] += triNormal;
            
        }
        
        for (size_t i=0; i<targetSurface->points.size(); i++)
            normals[i].normalize();
	    
    }

}

template <class ctype>
void NormalProjector<ctype>::insertEdge(PSurfaceFactory<2,ctype>& factory,
                                        int from, int to, 
                                        std::vector<PathVertex<ctype> >& edgePath)
{
    
    // start inserting the edge path by starting at the "to" node
    while(edgePath.size() > 1) 
    {

        // Connect two nodes that are on the same triangle
        if (onSameTriangle(edgePath.back().bundle_, edgePath[0].bundle_)) {
            //assert(edgePath.size()==2);    

            // Get the common triangle
            std::vector<int> commonTris = getCommonTris(edgePath.back().bundle_, edgePath[0].bundle_);

            for (int i=0; i<commonTris.size(); i++) 
                psurface_->triangles(commonTris[i]).addEdge(edgePath.back().bundle_.triToIdx(commonTris[i]), 
                                                        edgePath[0].bundle_.triToIdx(commonTris[i]));
            break;
        }
        // insert next edge path segment and remove the last node
        insertEdgeSegment(factory, from, to, edgePath);
        edgePath.pop_back();

    }
}


template <class ctype>
void NormalProjector<ctype>::insertEdgeSegment(PSurfaceFactory<2,ctype>& factory, int from, int to,
                                               std::vector<PathVertex<ctype> >& edgePath)
{
        //assert(edgeNodes.size()>2);

        // the next node on that edge path
        PathVertex<ctype>& node = edgePath[edgePath.size()-2];
        PathVertex<ctype>& lastNode = edgePath.back();

        const Surface* surf = psurface_->surface;
        
        // the image point of that edge path point 
        StaticVector<ctype,3> image;

        // trick: within Amira this is assignment from a McVec3f
        for (int j=0; j<3; j++)
            image[j] = surf->points[from][j] + node.lambda_*(surf->points[to][j]-surf->points[from][j]);

        // get neighboring triangle
        int tri = node.tri_;

        // edge goes through an intersection node
        if (node.corner_==-1) {

            // edge the node lives on
            int edge = node.edge_;
            // neighboring triangle
            int neighboringTri = lastNode.tri_;

            // if no neighboring triangle --> error
            if (neighboringTri==-1) {
                std::cout << "[FromInsertEdgeSegment] Warning: Normal images leaves domain surface!" << std::endl;
                return;
            }

            // add intersection nodes on both sides

            // the neighboring edge is the entering edge of the next edge path node
            int e = lastNode.enteringEdge_;

            // the domain position of the new intersection node on the two triangles
            ctype mu = node.locEdge_;
            StaticVector<ctype,2> dom1((edge==0)*(1-mu) + (edge==2)*mu, (edge==0)*mu + (edge==1)*(1-mu));
            StaticVector<ctype,2> dom2((e==0)*mu + (e==2)*(1-mu), (e==0)*(1-mu) + (e==1)*mu);

            // now that we created the node, we can add the bundle to the edge path
            node.bundle_ = factory.addIntersectionNodePair(tri, neighboringTri,
                    dom1, dom2, edge, e, image);

            // insert new parameter edge
            psurface_->triangles(neighboringTri).addEdge(lastNode.bundle_.triToIdx(neighboringTri), node.bundle_[1].idx);

            // else edge goes through a ghost node
        } else {
            int corner = node.corner_;

            DomainTriangle<ctype>& cT = psurface_->triangles(tri);
            node.bundle_ = psurface_->getNodeBundleAtVertex(cT.vertices[corner]);

            // for touching nodes it may be the case that the ghost and touching node are on the same edge
            // case that the edge path is aligned with a domain edge
            if( ((lastNode.type_ == Node<ctype>::TOUCHING_NODE) 
                && ( node.bundle_.triToIdx(lastNode.tri_) == 
                        ((psurface_->triangles(lastNode.tri_).nodes[lastNode.bundle_.triToIdx(lastNode.tri_)].getDomainEdge()+2)%3))  ) 
                    || (lastNode.type_ == Node<ctype>::CORNER_NODE)
                        || (lastNode.type_ == Node<ctype>::GHOST_NODE) ) {

                std::vector<int> commonTris = getCommonTris(node.bundle_, lastNode.bundle_);
                for (size_t i=0; i<commonTris.size(); i++)
                    psurface_->triangles(commonTris[i]).addEdge(lastNode.bundle_.triToIdx(commonTris[i]), 
                                                          node.bundle_.triToIdx(commonTris[i]));
            } else
                psurface_->triangles(lastNode.tri_).addEdge(lastNode.bundle_.triToIdx(lastNode.tri_), 
                    node.bundle_.triToIdx(lastNode.tri_));
        }
}

template <class ctype>
bool NormalProjector<ctype>::edgeCanBeInserted(const std::vector<StaticVector<ctype,3> >& normals,
                                        int from, int to, 
                                        const std::vector<NodeBundle>& projectedTo, 
                                        std::vector<PathVertex<ctype> >& edgePath)
{
    if (projectedTo[from].size() == 0 || projectedTo[to].size() == 0)
        return false;

    int enteringEdge=-1;
    NodeBundle curr = projectedTo[from];
    
    // initialize first node on the edge path
    edgePath.resize(1);
    edgePath[0].bundle_ = curr;
    edgePath[0].type_ = psurface_->nodes(curr[0]).type;
    edgePath[0].tri_ = curr[0].tri;

    // If the two nodes are on the same triangle it is surely possible to enter the edge
    if (onSameTriangle(curr, projectedTo[to])) {

        edgePath.push_back(PathVertex<ctype>(projectedTo[to]));
        edgePath.back().type_ = psurface_->nodes(projectedTo[to][0]).type;

        // set triangle the both are contained in,
        // if there is more than one we'll handle it later
        for (int i=0; i<curr.size(); i++)
            for (int j=0; j<edgePath.back().bundle_.size(); j++)
                if (curr[i].tri==edgePath.back().bundle_[j].tri) {
                    edgePath.back().tri_ = curr[i].tri;
                    return true;
                } 
    }


    typename Node<ctype>::NodeType currType = psurface_->nodes(curr[0]).type;
    int currTri     = curr[0].tri;

    // parameter value for the edge to be inserted
    ctype lambda = 0;

    // sometimes more than one domain edge hits the target edge
    // in this case we might need to restart the search from an earlier 
    // intersection node and start looking at the other edge first
    int offset = 0;
    
    // avoid running into infinity-loops
    PathVertex<ctype> wrongEdgeNode;

    while (true) {

        // If the two nodes are on the same triangle it is surely possible to enter the edge
        if (onSameTriangle(currTri, projectedTo[to])) {
            edgePath.push_back(PathVertex<ctype>(projectedTo[to]));
            edgePath.back().tri_ = currTri;
            edgePath.back().type_ = psurface_->nodes(projectedTo[to][0]).type;
            edgePath.back().enteringEdge_ = enteringEdge;
            return true;
        }

        if ((currType==Node<ctype>::GHOST_NODE || currType==Node<ctype>::CORNER_NODE)
                && onSameTriangle(curr, projectedTo[to])) {
                
            edgePath.push_back(PathVertex<ctype>(projectedTo[to]));
            edgePath.back().type_ = psurface_->nodes(projectedTo[to][0]).type;
            edgePath.back().enteringEdge_ = enteringEdge;
            // set triangle the both are contained in,
            // if there is more than one we'll handle it later
            for (int i=0; i<curr.size(); i++)
                for (int j=0; j<edgePath.back().bundle_.size(); j++)
                    if (curr[i].tri==edgePath.back().bundle_[j].tri) {
                        edgePath.back().tri_ = curr[i].tri;
                        return true;
                } 
        }


        switch (currType) {
        case Node<ctype>::TOUCHING_NODE:

            if (!testInsertEdgeFromTouchingNode(normals, from, to, lambda, curr, currType, 
                        currTri, enteringEdge, edgePath, offset))
                return false;
            break;

        case Node<ctype>::GHOST_NODE:
        case Node<ctype>::CORNER_NODE:

            try {
                if (!testInsertEdgeFromCornerNode(normals, from, to, lambda,
                                              curr, currType, currTri, enteringEdge,edgePath, offset))
                return false;
            } catch (WrongEdgeException e) {

                std::cout<<"Exception caught! "<<e.what()<<std::endl;
                
                // if we were already ran into that node, then the projection of the path does not exist
                if (edgePath.back()==wrongEdgeNode)
                    return false;

                // this edge led to a wrong triangle so check if another edge is also feasible
                offset = (edgePath.back().edge_+1)%3;
                currTri = edgePath.back().tri_;
                currType = edgePath.back().type_;
                lambda = edgePath[edgePath.size()-2].lambda_;
                enteringEdge = edgePath.back().enteringEdge_;
                wrongEdgeNode = edgePath.back();
                edgePath.pop_back();
            }
            break;

        case Node<ctype>::INTERSECTION_NODE:

            try {
                if (!testInsertEdgeFromIntersectionNode(normals, from, to, lambda, 
                        curr, currType, currTri, enteringEdge, edgePath, offset))
                return false;
            } catch (WrongEdgeException e) {

                std::cout<<"Exception caught! "<<e.what()<<std::endl;

                // if we were already ran into that node, then the projection of the path does not exist
                if (edgePath.back()==wrongEdgeNode)
                    return false;

                // this edge led to a wrong triangle so check if another edge is also feasible
                offset = (edgePath.back().edge_+1)%3;
                currTri = edgePath.back().tri_;
                currType = edgePath.back().type_;
                lambda = edgePath[edgePath.size()-2].lambda_;
                enteringEdge = edgePath.back().enteringEdge_;
                wrongEdgeNode = edgePath.back();
                edgePath.pop_back();
            }
            break;
            
        case Node<ctype>::INTERIOR_NODE:

            if (!testInsertEdgeFromInteriorNode(normals, from, to, lambda, 
                        curr, currType, currTri, enteringEdge, edgePath, offset))
                return false;

            break;

        default:
            std::cout << "ERROR: unknown node type found!" << std::endl;
            abort();
        }
    }
    
    //std::cout << "should not occur" << std::endl;
    return true;
}


template <class ctype>
bool NormalProjector<ctype>::testInsertEdgeFromInteriorNode(const std::vector<StaticVector<ctype,3> >& normals,
                                                     int from, int to, ctype &lambda,
                                                     NodeBundle& curr,
                                                     typename Node<ctype>::NodeType& currType, int& currTri,
                                                     int& enteringEdge, std::vector<PathVertex<ctype> >& edgePath, 
                                                     int offset)
{
    // loop over the three edges of the current triangle (except for the entering edge) and
    // check whether the paramPolyEdge leaves the triangle via this edge
    ctype eps = 1e-5;
    int i=offset;
    for (int j=0; j<3; j++,i=(i+1)%3) {
            
        if (i==enteringEdge)
            continue;
            
        StaticVector<ctype,3> x;
        int p = psurface_->triangles(currTri).vertices[i];
        int q = psurface_->triangles(currTri).vertices[(i+1)%3];

        const Surface* surf = psurface_->surface;

        StaticVector<ctype,3> targetFrom, targetTo;
        for (int j=0; j<3; j++) {
            targetFrom[j] = surf->points[from][j];
            targetTo[j]   = surf->points[to][j];
        }

        if (edgeIntersectsNormalFan(targetFrom, targetTo,
                                    psurface_->vertices(p), psurface_->vertices(q),
                                    normals[p], normals[q], x)) {

            const ctype& newLambda = x[1];
            const ctype& mu        = x[0];
                
            if (newLambda < lambda) {
                // Error: the normal projection is not continuous!
                return false;
            }

            int corner = -1;
            if (mu<eps) 
                corner = i;
            else if (mu>1-eps)
                corner = (i+1)%3;

            if (corner==-1) {
                // get neighboring triangle
                int neighboringTri = psurface_->getNeighboringTriangle(currTri, i);
                    
                // if no neighboring triangle --> error
                if (neighboringTri==-1) {
                    // Error: Normal images leaves domain surface!
                    return false;
                }
                    
                // add intersection nodes on both sides
                    
                // the domain position of the new intersection node on the next triangle
                // better: using getEdge()
                int e = psurface_->triangles(neighboringTri).getCorner(q);
                    
                currType = Node<ctype>::INTERSECTION_NODE;

                // add node on the path to the array
                edgePath.push_back(PathVertex<ctype>(currTri,i,mu,currType,NodeBundle(),newLambda,enteringEdge));

                currTri  = neighboringTri;
                lambda   = newLambda;
                enteringEdge = e;

    
                return true;
                
            } else {
             
                // parameter polyedge is leaving base grid triangle through a ghost node
                
                // get all ghost nodes for the base grid vertex
                int vertex = psurface_->triangles(currTri).vertices[corner];
                std::vector<int> neighbors = psurface_->getTrianglesPerVertex(vertex);
                
                curr.resize(0);
                for (int k=0; k<neighbors.size(); k++) {
                    
                    int cornerOnNeighbor = psurface_->triangles(neighbors[k]).getCorner(vertex);
                    
                    /** \todo Linear search: pretty slow */
                    for (int l=0; l<psurface_->triangles(neighbors[k]).nodes.size(); l++) {
                        
                        if (psurface_->triangles(neighbors[k]).nodes[l].isGHOST_NODE()
                            && psurface_->triangles(neighbors[k]).nodes[l].getCorner() == cornerOnNeighbor){
                            
                            curr.push_back(GlobalNodeIdx(neighbors[k], l));
                            break;
                            
                        }
                        
                    }
                    
                }
                
                currType = Node<ctype>::GHOST_NODE;
                assert(currType==type(curr));
                lambda = newLambda;
                edgePath.push_back(PathVertex<ctype>(currTri,i,mu,currType,curr, newLambda, enteringEdge, corner));
                return true;
   
            }

        }
            
    }

    return false;
}



template <class ctype>
bool NormalProjector<ctype>::testInsertEdgeFromIntersectionNode(const std::vector<StaticVector<ctype,3> >& normals,
                                                         int from, int to, ctype &lambda,
                                                         NodeBundle& curr,
                                                         typename Node<ctype>::NodeType& currType, int& currTri,
                                                         int& enteringEdge,  
                                                         std::vector<PathVertex<ctype> >& edgePath, 
                                                         int offset)
{
    // loop over the three edges of the current triangle (except for the entering edge) and
    // check whether the paramPolyEdge leaves the triangle via this edge
    ctype eps = 1e-5;
    
    int i=offset;
    for (int j=0; j<3; j++,i=(i+1)%3) {
            
        if (i==enteringEdge)
            continue;
            
        StaticVector<ctype,3> x;
        int p = psurface_->triangles(currTri).vertices[i];
        int q = psurface_->triangles(currTri).vertices[(i+1)%3];
        
        const Surface* surf = psurface_->surface;

        StaticVector<ctype,3> targetFrom, targetTo;
        for (int j=0; j<3; j++) {
            targetFrom[j] = surf->points[from][j];
            targetTo[j]   = surf->points[to][j];
        }

        if (edgeIntersectsNormalFan(targetFrom, targetTo,
                                    psurface_->vertices(p), psurface_->vertices(q),
                                    normals[p], normals[q], x)) {

            const ctype& newLambda = x[1];
            const ctype& mu        = x[0];
                
            if (newLambda < lambda) {
                // Error: the normal projection is not continuous!
                return false;                    
            }

            int corner = -1;
            if (mu<eps) 
                corner = i;
            else if (mu>1-eps)
                corner = (i+1)%3;

            if (corner==-1) {
                // get neighboring triangle
                int neighboringTri = psurface_->getNeighboringTriangle(currTri, i);
                    
                // if no neighboring triangle --> error
                if (neighboringTri==-1) {
                    // Error: Normal images leaves domain surface!
                    return false;
                }
                    
                // add intersection nodes on both sides
                    
                // the domain position of the new intersection node on the next triangle
                // better: using getEdge()
                int e = psurface_->triangles(neighboringTri).getCorner(q);

                currType = Node<ctype>::INTERSECTION_NODE;
                edgePath.push_back(PathVertex<ctype>(currTri, i, mu, currType,
                                NodeBundle(), newLambda, enteringEdge));
                currTri  = neighboringTri;
                lambda   = newLambda;
                enteringEdge = e;
                

                return true;

            } else {

                // parameter polyedge is leaving base grid triangle through a ghost node
                
                // get all ghost nodes for the base grid vertex
                int vertex = psurface_->triangles(currTri).vertices[corner];
                curr = psurface_->getNodeBundleAtVertex(vertex);
                currType = Node<ctype>::GHOST_NODE;
                assert(currType==type(curr));
                lambda = newLambda;
                
                edgePath.push_back(PathVertex<ctype>(currTri,i,mu,currType,curr, 
                                newLambda, enteringEdge, corner));
                return true;

            }

        }
    }

    throw(WrongEdgeException("No intersection found (in testInsertEdgeFromIntersectionNode)!\n"));
    return false;
}



template <class ctype>
bool NormalProjector<ctype>::testInsertEdgeFromTouchingNode(const std::vector<StaticVector<ctype,3> >& normals,
                                                     int from, int to, ctype &lambda,
                                                     NodeBundle& curr,
                                                     typename Node<ctype>::NodeType& currType, int& currTri,
                                                     int& enteringEdge, std::vector<PathVertex<ctype> >& edgePath,
                                                     int offset)
{
    const Surface* surf = psurface_->surface;
    ctype eps = 1e-5;
    
    // The other end of the edge is *not* on this triangle
    for (int i=0; i<curr.size(); i++) {

        const DomainTriangle<ctype>& cT = psurface_->triangles(curr[i].tri);
        int currentEdge = cT.nodes[curr[i].idx].getDomainEdge();

        int j=offset;
        for (int k=0; k<3; k++,j=(j+1)%3) {
        //for (int j=0; j<3; j++) {
            if (j==currentEdge)
                continue;
        
            StaticVector<ctype,3> x;
            int p = cT.vertices[j];
            int q = cT.vertices[(j+1)%3];

            StaticVector<ctype,3> targetFrom, targetTo;
            for (int l=0; l<3; l++) {
                targetFrom[l] = surf->points[from][l];
                targetTo[l]   = surf->points[to][l];
            }

            if (edgeIntersectsNormalFan(targetFrom, targetTo,
                                        psurface_->vertices(p), psurface_->vertices(q),
                                        normals[p], normals[q], x)) {
                
                const ctype& newLambda = x[1];
                
                if (newLambda < lambda) {
                    // Edge insertion not possible: the normal projection is not continuous!
                    return false;                    
                }
                
                int corner = -1;
                if (x[0]<eps) 
                    corner = j;
                else if (x[0]>1-eps)
                    corner = (j+1)%3;
                
                if (corner==-1) {
                    // parameter polyedge is leaving basegrid triangle through an edge
                    
                    // get neighboring triangle
                    int neighboringTri = psurface_->getNeighboringTriangle(curr[i].tri, j);

                    // if no neighboring triangle --> error
                    if (neighboringTri==-1)
                        return false;
                    
                    // add intersection nodes on both sides
                
                    // better: using getEdge()
                    int e = psurface_->triangles(neighboringTri).getCorner(q);
                     
                    currType = Node<ctype>::INTERSECTION_NODE;
                    edgePath.push_back(PathVertex<ctype>(curr[i].tri,j,x[0],currType,NodeBundle(), 
                                            newLambda, enteringEdge));
                    currTri  = neighboringTri;
                    lambda   = newLambda;
                    enteringEdge = e;
                    

                    return true;
                
                } else {
                    // parameter polyedge is leaving base grid triangle through a ghost node

                    // get all ghost nodes for the base grid vertex
                    int vertex = psurface_->triangles(curr[i].tri).vertices[corner];
                    int copyTri = curr[i].tri;
                    curr = psurface_->getNodeBundleAtVertex(vertex); 
                    currType = Node<ctype>::GHOST_NODE;
                    assert(currType==type(curr));
                    lambda = newLambda;

                    edgePath.push_back(PathVertex<ctype>(copyTri,j,x[0],currType,curr,
                                    newLambda, enteringEdge,corner));
                    return true;

                }
                
            }
            
        }

    }

    return false;

}



template <class ctype>
bool NormalProjector<ctype>::testInsertEdgeFromCornerNode(const std::vector<StaticVector<ctype,3> >& normals, 
                                                   int from, int to, ctype &lambda,
                                                   NodeBundle& curr, 
                                                   typename Node<ctype>::NodeType& currType, int& currTri,
                                                   int& leavingEdge, std::vector<PathVertex<ctype> >& edgePath,
                                                   int offset)
{
    const Surface* surf = psurface_->surface;
    ctype eps = 1e-5;

    // The other end of the edge is *not* on this triangle
    for (int i=0; i<curr.size(); i++) {

        int cT = curr[i].tri;

        int thisCorner = psurface_->triangles(cT).nodes[curr[i].idx].getCorner();
        int oppEdge = (thisCorner+1)%3;

        StaticVector<ctype,3> x;
        int p = psurface_->triangles(cT).vertices[(thisCorner+1)%3];
        int q = psurface_->triangles(cT).vertices[(thisCorner+2)%3];

        StaticVector<ctype,3> targetFrom, targetTo;
        for (int j=0; j<3; j++) {
            targetFrom[j] = surf->points[from][j];
            targetTo[j]   = surf->points[to][j];
        }

        if (edgeIntersectsNormalFan(targetFrom, targetTo,
                                    psurface_->vertices(p), psurface_->vertices(q),
                                    normals[p], normals[q], x)) {
            
            const ctype& newLambda = x[1];
            
            if (newLambda < lambda) {
                // Shouldn't this rather be a 'return false' here?
                continue;
            }
            
            int corner = -1;
            if (x[0]<eps) 
                corner = (thisCorner+1)%3;
            else if (x[0]>1-eps)
                corner = (thisCorner+2)%3;
            
            if (corner==-1) {
                // parameter polyedge is leaving basegrid triangle 
                // through the opposite edge
                
                // get neighboring triangle
                int neighboringTri = psurface_->getNeighboringTriangle(cT, oppEdge);

                // if no neighboring triangle --> error
                if (neighboringTri==-1)
                    return false;

                // add intersection nodes on both sides
                
                // better: using getEdge()
                int e = psurface_->triangles(neighboringTri).getCorner(q);

                currType = Node<ctype>::INTERSECTION_NODE;
                edgePath.push_back(PathVertex<ctype>(cT,oppEdge, x[0],currType,NodeBundle(),
                                    leavingEdge,newLambda));
                currTri  = neighboringTri;
                lambda   = newLambda;
                leavingEdge = e;


                return true;

            } else {

                // parameter polyedge is leaving base grid triangle through a ghost node

                // get all ghost nodes for the base grid vertex
                int vertex = psurface_->triangles(curr[i].tri).vertices[corner];
                curr = psurface_->getNodeBundleAtVertex(vertex);
                currType = Node<ctype>::GHOST_NODE;
                assert(currType == type(curr));
                lambda = newLambda;
                
                edgePath.push_back(PathVertex<ctype>(cT,oppEdge,x[0],currType,curr, leavingEdge, newLambda, corner));
                return true;
                
            }
            
        }
        
    }

    // TODO think about the case below
    throw(WrongEdgeException("No intersection found (in testInsertEdgeFromCornerNode)!\n"));

    // If we get to here this means that the intersection of the edge from 'from' to 'to' with
    // the star around the vertex corresponding to the corner node we are testing consists
    // only of a single point.  This can happen if the star is not a full circle.  In this
    // case the edge cannot be inserted and consequently we return 'false'.
    return false;
}



template <class ctype>
bool NormalProjector<ctype>::onSameTriangle(const NodeBundle& a, const NodeBundle& b) const
{
    for (int i=0; i<a.size(); i++)
        for (int j=0; j<b.size(); j++)
            if (a[i].tri==b[j].tri)
                return true;

    return false;
}


template <class ctype>
bool NormalProjector<ctype>::onSameTriangle(const int& tri, const NodeBundle& b) const
{
    for (size_t j=0; j<b.size(); j++)
        if (tri==b[j].tri)
            return true;

    return false;
}


template <class ctype>
bool NormalProjector<ctype>::computeInverseNormalProjection(const StaticVector<ctype,3>& p0, const StaticVector<ctype,3>& p1, const StaticVector<ctype,3>& p2,
                                                     const StaticVector<ctype,3>& n0, const StaticVector<ctype,3>& n1, const StaticVector<ctype,3>& n2,
                                                     const StaticVector<ctype,3>& target, StaticVector<ctype,3>& x)
{
    const ctype eps = 1e-6;

    // try to solve a cubic equation for the distance parameter nu, then compute the barycentric coordinates from it

    // cubic coefficient
    StaticVector<ctype,3> n02 = n0 - n2;
    StaticVector<ctype,3> n12 = n1 - n2;
    StaticVector<ctype,3> n02n12 = n02.cross(n12);

    double cubic = n2.dot(n02n12);

    // quadratic coefficient

    StaticVector<ctype,3> p02 = p0 - p2;
    StaticVector<ctype,3> p12 = p1 - p2;
    StaticVector<ctype,3> p2q = p2 -target;
    StaticVector<ctype,3> p02n12 = p02.cross(n12);
    StaticVector<ctype,3> n02p12 = n02.cross(p12);

    ctype quadratic = n2.dot(p02n12)+n2.dot(n02p12)+p2q.dot(n02n12);

    // constant coefficient
    StaticVector<ctype,3> p02p12 = p02.cross(p12);
    ctype constant = p2q.dot(p02p12);

    // linear coefficient
    ctype linear = n2.dot(p02p12) + p2q.dot(n02p12) + p2q.dot(p02n12);

    // save all zeros we find
    std::vector<ctype> zeros;

    if (std::fabs(cubic) <1e-10 && std::fabs(quadratic)<1e-10 && std::fabs(linear)<1e-10) {
        return false;
    } else if (std::fabs(cubic) <1e-10 && std::fabs(quadratic)<1e-10) {

        // problem is linear
        zeros.push_back(-constant/linear);

    } else if(std::fabs(cubic)<1e-10) {

        // problem is quadratic
        ctype p = linear/quadratic;
        ctype q = constant/quadratic;

        ctype sqt = 0.25*p*p -q;

        // no real solution
        if (sqt<-1e-10)
            return false;

        zeros.push_back(-0.5*p + std::sqrt(sqt));
        zeros.push_back(-0.5*p -std::sqrt(sqt));

    } else {

        // problem is cubic
        quadratic /= cubic;
        linear /= cubic;
        constant /= cubic;

        // Transform to reduced form z^3 + p*z + q = 0 where x = z-quadratic/3
        ctype p= linear - quadratic*quadratic/3;
        ctype q=quadratic*(2*quadratic*quadratic/27 - linear/3) + constant;

        // use Cardano's method to solve the problem
        ctype D = 0.25*q*q + std::pow(p,3)/27;

        if (D>1e-10) {
            // one real zero

            // be careful when computing the cubic roots
            ctype nu = -q/2+std::sqrt(D);
            ctype zer = std::pow(std::fabs(nu),1.0/3.0) * ((nu<-1e-10) ? -1 : 1);

            nu = -q/2-std::sqrt(D);
            zer += std::pow(std::fabs(nu),1.0/3.0) * ((nu<-1e-10) ? -1 : 1);

            zeros.push_back(zer-quadratic/3);

        } else if (D<-1e-10) {

            // three real zeros, using trigonometric functions to compute them
            ctype a = std::sqrt(-4*p/3);
            ctype b = std::acos(-0.5*q*std::sqrt(-27/(std::pow(p,3))));

            for (int i=0;i<3; i++)
                zeros.push_back(std::pow(-1,i+1)*a*std::cos((b+(1-i)*M_PI)/3) -quadratic/3);


        } else {
            // one single and one double zero

            if (std::fabs(q)<1e-10) {
                zeros.push_back(-quadratic/3);

                if (p<-1e-10)
                    zeros.push_back(std::sqrt(-p)-quadratic/3);

            } else if (std::fabs(p)<1e-10) { // is this case correct?

                double nu = std::pow(std::fabs(q),1.0/3.0) * ((q<-eps) ? -1 : 1);
                zeros.push_back(nu-quadratic/3);

            } else {
                zeros.push_back(3*q/p - quadratic/3);
                zeros.push_back(-1.5*q/p - quadratic/3);
            }
        }
    }

    int index = -1;
    StaticVector<ctype,3> r;
    std::vector<StaticVector<ctype,3> > lamb(zeros.size());
    for (int i=0;i<zeros.size();i++) {

        ctype nu=zeros[i];
        // only look in the direction of the outer normals
        if (nu<-1e-1) // allowed overlaps
            continue;

        if (index != -1)
            if (nu > zeros[index]) // is this one really closer ?
                continue;

        r[2] = nu;

        // the computation of the other components might lead to nan or inf
        // if this happens use a different equation to compute them
        StaticVector<ctype,3> c = (p2q+nu*n2).cross(p02+nu*(n02));
        StaticVector<ctype,3> d = (p02 +nu*n02).cross(p12+nu*n12);
        StaticVector<ctype,3> e = p2q+nu*n2;
        StaticVector<ctype,3> f = p12+nu*n12;
        StaticVector<ctype,3> g = p02+nu*n02;

        // computation of the other components is unstable
        for (int j=0;j<3; j++) {

            r[1] = c[j]/d[j];

            if (isnan(r[1]) || isinf(r[1]))
                continue;

            r[0] = -(e[(j+1)%3]+r[1]*f[(j+1)%3])/g[(j+1)%3];

            if (!(isnan(r[0]) || isinf(r[0])) && (p2q +r[0]*p02 + r[1]*p12 + r[2]*r[0]*n02+r[2]*r[1]*n12+r[2]*n2).length()<1e-3)
                break;

            r[0] = -(e[(j+2)%3]+r[1]*f[(j+2)%3])/g[(j+2)%3];


            if (!(isnan(r[0]) || isinf(r[0])) && (p2q + r[0]*p02 + r[1]*p12 + r[2]*r[0]*n02+r[2]*r[1]*n12+r[2]*n2).length()<1e-3)
                break;

        }
        lamb[i] = r;
        if (r[0] > -eps && r[1]> -eps && (r[0]+r[1] < 1+eps)) {
            index = i;
            x = r;
        }
    }

    StaticVector<ctype,3> res = p2q + x[0]*p02 + x[1]*p12 + x[2]*x[0]*n02+x[2]*x[1]*n12+x[2]*n2;

    if (res.length()<1e-6) {
        if (index >= 0)
            return true;

        return false;
    }

    //std::cout<<"Direct solution failed, use Newton method\n";

    StaticVector<ctype,3> oldX = x;

    // Fix some initial value
    // Some problems have two solutions and the Newton converges to the wrong one
    x.assign(0.5);

    for (int i=0; i<30; i++) {

        // compute Newton correction
        StaticVector<ctype,3> Fxk = x[0]*(p0-p2) + x[1]*(p1-p2) + x[2]*x[0]*(n0-n2) + x[2]*x[1]*(n1-n2) + x[2]*n2 + p2 - target;

        StaticMatrix<ctype,3> FPrimexk(p0 - p2 + x[2]*(n0-n2),
                         p1 - p2 + x[2]*(n1-n2),
                         x[0]*(n0-n2) + x[1]*(n1-n2) + n2);

        StaticMatrix<ctype,3> FPrimexkInv = FPrimexk.inverse();

        StaticVector<ctype,3> newtonCorrection; // = (-1) * FPrimexk.inverse() * Fxk;

        FPrimexkInv.multMatrixVec(-Fxk, newtonCorrection);

        x += newtonCorrection;

    }

    StaticVector<ctype,3> res2 = p2q + x[0]*p02 + x[1]*p12 + x[2]*x[0]*n02+x[2]*x[1]*n12+x[2]*n2;

    if (x[0]>-eps && x[1]>-eps && (x[0]+x[1] <1+eps)) {

        // Newton did not converge either
        if (res2.length()>1e-6)
            return false;

        return true;
    } 

    return false;
}


template <class ctype>
bool NormalProjector<ctype>::edgeIntersectsNormalFan(const StaticVector<ctype,3>& q0, const StaticVector<ctype,3>& q1,
                                              const StaticVector<ctype,3>& p0, const StaticVector<ctype,3>& p1,
                                              const StaticVector<ctype,3>& n0, const StaticVector<ctype,3>& n1,
                                              StaticVector<ctype,3>& x)
{
    int i;

    ctype eps = 1e-6;
    // solve a quadratic scalar equation for the distance parameter eta, then compute the barycentric coordinates from it

    StaticVector<ctype,3> n10 = n1 - n0;
    StaticVector<ctype,3> p10 = p1 - p0;
    StaticVector<ctype,3> q10 = q1 - q0;
    StaticVector<ctype,3> q10n10 = q10.cross(n10);
    StaticVector<ctype,3> q10p10 = q10.cross(p10);
    StaticVector<ctype,3> p0q0 = p0 - q0;

    // quadratic coefficient
    ctype quadratic = n0.dot(q10n10);

    // linear coefficient
    ctype linear = n0.dot(q10p10) + p0q0.dot(q10n10);

    // constant coefficient
    ctype constant = p0q0.dot(q10p10);

    // save all zeros we find
    std::vector<ctype> zeros;

    if (std::fabs(quadratic)<1e-10 && std::fabs(linear)<1e-10) {
        return false;
    } else if (std::fabs(quadratic)<1e-10) {

        // problem is linear
        zeros.push_back(-constant/linear);

    } else {

        // problem is quadratic
        ctype p = linear/quadratic;
        ctype q = constant/quadratic;

        ctype sqt = 0.25*p*p -q;

        // no real solution
        if (sqt<-1e-10)
            return false;

        zeros.push_back(-0.5*p + std::sqrt(sqt));
        zeros.push_back(-0.5*p -std::sqrt(sqt));

    }

    int index = -1;
    StaticVector<ctype,3> r;
    std::vector<StaticVector<ctype,3> > lamb(zeros.size());
    for (int i=0;i<zeros.size();i++) {

        ctype eta=zeros[i];

        // only look in the direction of the outer normals
        if (eta<-1e-1)
            continue;

        r[2] = eta;

        // the computation of the other components might lead to nan or inf
        // if this happens use a different equation to compute them
        StaticVector<ctype,3> c =(p0q0+eta*n0).cross(q10);
        StaticVector<ctype,3> d =q10.cross(p10+eta*n10);
        for (int j=0;j<3; j++) {

            r[0] = c[j]/d[j];
            if (isnan(r[0]) || isinf(r[0]))
                continue;

            r[1] = (p0q0[(j+1)%3]+eta*n0[(j+1)%3] + r[0]*(p10[(j+1)%3]+eta*n10[(j+1)%3]))/q10[(j+1)%3];

            // computation of the other components can be instable
            if (!(isnan(r[1]) || isinf(r[1])) && (p0q0 + r[0]*p10 + r[2]*n0 +r[2]*r[0]*n10 -r[1]*q10).length()<1e-3  )
                break;

            r[1] = (p0q0[(j+2)%3]+eta*n0[(j+2)%3] + r[0]*(p10[(j+2)%3]+eta*n10[(j+2)%3]))/q10[(j+2)%3];

            // computation of the other components can be instable
            if (!(isnan(r[1]) || isinf(r[1])) && (p0q0 + r[0]*p10 + r[2]*n0 + r[2]*r[0]*n10 -r[1]*q10).length()<1e-3)
                break;

        }
        lamb[i] = r;
        if (r[0] >= -eps && r[1]>= -eps && (r[0]<=1+eps)  && (r[1] <= 1+eps)) {
            index = i;
            x = r;
        }

    }

    StaticVector<ctype,3> res = p0q0 + x[0]*p10 + x[2]*n0 + x[2]*x[0]*n10 -x[1]*q10;
    if (res.length()<eps)
    {
        if (index >= 0)
            return true;

        return false;
    }

    // if the direct compuation failed, use a Newton method to compute at least one zero
    StaticVector<ctype,3> oldX = x;

    // Fix some initial value
    // sometimes it only works when the initial value is an intersection...
    x[0] = x[1] = 0.5;
    x[2] = 0.5;
    StaticVector<ctype,3> newtonCorrection;

    for (i=0; i<30; i++) {

        // compute Newton correction

        StaticVector<ctype,3> Fxk = p0-q0 + x[0]*(p1-p0) + x[2]*n0 + x[2]*x[0]*(n1-n0) - x[1]*(q1-q0);

        StaticMatrix<ctype,3> FPrimexk(p1-p0 + x[2]*(n1-n0),
                         q0-q1,
                         n0 + x[0]*(n1-n0));

        StaticMatrix<ctype,3> FPrimexkInv = FPrimexk.inverse();

        FPrimexkInv.multMatrixVec(-Fxk, newtonCorrection);

        x += newtonCorrection;

    }

    StaticVector<ctype,3> res2 = p0q0 + x[0]*p10 + x[2]*n0 + x[2]*x[0]*n10 -x[1]*q10;
    if (res2.length()<=eps) {

        if (x[0]>=-eps && x[0]<=(1+eps) && x[1]>=-eps && x[1]<=(1+eps))
            return true;

        return false;

    }

    std::cout<<"Newton did not converge either!\n";

    return false;
}


template <class ctype>
bool NormalProjector<ctype>::rayIntersectsTriangle(const StaticVector<ctype,3>& basePoint, const StaticVector<ctype,3>& direction,
                                            const StaticVector<ctype,3>& a, const StaticVector<ctype,3>& b, const StaticVector<ctype,3>& c,
                                            StaticVector<ctype,2>& localCoords, ctype& normalDist, ctype eps)
{
    const StaticVector<ctype,3> &p = basePoint;

    StaticVector<ctype,3> e1 = b-a;
    StaticVector<ctype,3> e2 = c-a;
    e1.normalize();
    e2.normalize();
    bool parallel = fabs(StaticMatrix<ctype,3>(e1, e2, direction).det()) <eps;

    // Cramer's rule

    if (!parallel){

        ctype det = StaticMatrix<ctype,3>(b-a, c-a, direction).det();

        // triangle and edge are not parallel
        ctype nu = StaticMatrix<ctype,3>(b-a, c-a, p-a).det() / det;
       
        // only allow a certain overlaps 
        if (nu>1e-2)
            return false;

        ctype lambda = StaticMatrix<ctype,3>(p-a, c-a, direction).det() / det;
        if (lambda<-eps) return false;

        ctype mu = StaticMatrix<ctype,3>(b-a, p-a, direction).det() / det;
        if (mu<-eps) return false;

        if (lambda + mu > 1+eps)
            return false;
        else {
            localCoords[0] = 1-lambda-mu;
            localCoords[1] = lambda;
            normalDist     = -nu;

            return true;
        }

    } else {

        // triangle and edge are parallel
        ctype alpha = StaticMatrix<ctype,3>(b-a, c-a, p-a).det();
        if (alpha<-eps || alpha>eps)
            return false;
        else {
            printf("ray and triangle are parallel!\n");
            return false;

        }

    }


}


template <class ctype>
NodeIdx NormalProjector<ctype>::getCornerNode(const DomainTriangle<ctype>& cT, int corner)
{
    assert(corner>=0 && corner<3);

    for (size_t i=0; i<cT.nodes.size(); i++)
        if ((cT.nodes[i].isCORNER_NODE() || cT.nodes[i].isGHOST_NODE()) &&
                cT.nodes[i].getCorner()==corner)
            return i;

    return -1;
}

template <class ctype>
typename Node<ctype>::NodeType NormalProjector<ctype>::type(const NodeBundle& b) const
{
    assert(b.size()>0);
    for (size_t i=0; i<b.size(); i++)
        assert(psurface_->nodes(b[i]).type == psurface_->nodes(b[0]).type);
    return psurface_->nodes(b[0]).type;
}


template <class ctype>
int NormalProjector<ctype>::getCommonTri(const NodeBundle& a, const NodeBundle& b)
{
    for (size_t i=0; i<a.size(); i++)
        for (size_t j=0; j<b.size(); j++)
            if (a[i].tri==b[j].tri)
                return a[i].tri;

    return -1;
}


template <class ctype>
std::vector<int> NormalProjector<ctype>::getCommonTris(const NodeBundle& a, const NodeBundle& b)
{
    std::vector<int> result;

    for (size_t i=0; i<a.size(); i++)
        for (size_t j=0; j<b.size(); j++)
            if (a[i].tri==b[j].tri)
                result.push_back(a[i].tri);

    return result;
}



template <class ctype>
void NormalProjector<ctype>::setupEdgePointArrays()
{
    int i, j;

    for (i=0; i<psurface_->getNumTriangles(); i++) {

        DomainTriangle<ctype>& cT = psurface_->triangles(i);

        cT.edgePoints[0].clear();
        cT.edgePoints[1].clear();
        cT.edgePoints[2].clear();

        for (j=0; j<cT.nodes.size(); j++) {

            Node<ctype>& cN = cT.nodes[j];

            if (cN.isINTERIOR_NODE())
                continue;

            if (cN.isCORNER_NODE() || cN.isGHOST_NODE()) {
                int corner = cN.getCorner();
                cT.edgePoints[corner].insert(cT.edgePoints[corner].begin(), j);
                cT.edgePoints[(corner+2)%3].push_back(j);
                continue;
            } 

            ctype lambda = cN.getDomainEdgeCoord();
            int domainEdge = cN.getDomainEdge();
            std::vector<int>& cEP = cT.edgePoints[domainEdge];

            int idx = 0;
            while (idx<cEP.size() && cT.nodes[cEP[idx]].getDomainEdgeCoord(domainEdge)<lambda) {
                idx++;
            }                

            cEP.insert(cEP.begin()+idx, j);

        }

    }   

}

// ///////////////////////////////////////////////////////////////
//   A few static methods for the 1d-in-2d case.
// ///////////////////////////////////////////////////////////////

template <class ctype>
bool NormalProjector<ctype>::computeInverseNormalProjection(const StaticVector<ctype,2>& p0,
                                                            const StaticVector<ctype,2>& p1,
                                                            const StaticVector<ctype,2>& n0,
                                                            const StaticVector<ctype,2>& n1,
                                                            const StaticVector<ctype,2>& q,
                                                            ctype& local)
{
    ctype a = (p1[1]-p0[1])*(n1[0]-n0[0]) - (p1[0]-p0[0])*(n1[1]-n0[1]);
    ctype b = -(q[1]-p0[1])*(n1[0]-n0[0]) + (p1[1]-p0[1])*n0[0] + (q[0]-p0[0])*(n1[1]-n0[1]) - (p1[0]-p0[0])*n0[1];
    ctype c = -(q[1]-p0[1])*n0[0] + (q[0]-p0[0])*n0[1];

    // Is the quadratic formula degenerated to a linear one?
    if (std::abs(a) < 1e-10) {
        local = -c/b;
        //printf("mu:  %g,  old local %g\n", mu, ((q[0]-p0[0]) / (p1[0]-p0[0])));

        return local >= 0 && local <= 1;
    }

    // The abc-formula
    ctype mu_0 = (-b + std::sqrt(b*b - 4*a*c))/(2*a);
    ctype mu_1 = (-b - std::sqrt(b*b - 4*a*c))/(2*a);

    if (mu_0 >= 0 && mu_0 <= 1) {
        local = mu_0;
        return true;
    } else if (mu_1 >= 0 && mu_1 <= 1) {
        local = mu_1;
        return true;
    }
    return false;
}

template <class ctype>
bool NormalProjector<ctype>::normalProjection(const StaticVector<ctype,2>& base,
                                              const StaticVector<ctype,2>& direction,
                                              int& bestSegment,
                                              ctype& rangeLocalPosition,
                                              const std::vector<std::tr1::array<int,2> >& targetSegments,
                                              const std::vector<std::tr1::array<ctype, 2> >& coords)
{
    bestSegment = -1;
    int nTargetSegments = targetSegments.size();
    ctype bestDistance = std::numeric_limits<ctype>::max();

    for (int i=0; i<nTargetSegments; i++) {

        StaticVector<ctype,2> p0, p1;
        p0[0] = coords[targetSegments[i][0]][0];
        p0[1] = coords[targetSegments[i][0]][1];

        p1[0] = coords[targetSegments[i][1]][0];
        p1[1] = coords[targetSegments[i][1]][1];

        ctype distance, targetLocal;
        if (!rayIntersectsLine(base, direction, p0, p1, distance, targetLocal))
            continue;

        if (distance < bestDistance) {
            bestDistance = distance;
            bestSegment  = i;
            rangeLocalPosition = targetLocal;
        }

    }

    return bestSegment != -1;
}

template <class ctype>
bool NormalProjector<ctype>::
rayIntersectsLine(const StaticVector<ctype, 2>& basePoint, 
                  const StaticVector<ctype, 2>& direction,
                  const StaticVector<ctype, 2>& a, 
                  const StaticVector<ctype, 2>& b, 
                  ctype& distance, ctype& targetLocal)
{
    // we solve the equation basePoint + x_0 * normal = a + x_1 * (b-a)

    StaticMatrix<ctype,2> mat;
    mat[0][0] = direction[0];
    mat[1][0] = direction[1];
    mat[0][1] = a[0]-b[0];
    mat[1][1] = a[1]-b[1];

    /** \todo Easier with expression templates */
    StaticVector<ctype,2> rhs;
    rhs[0] = a[0]-basePoint[0];
    rhs[1] = a[1]-basePoint[1];

    StaticVector<ctype,2> x;

    // Solve the system.  If it is singular the normal and the segment
    // are parallel and there is no intersection

    ctype detinv = mat[0][0]*mat[1][1]-mat[0][1]*mat[1][0];
    if (std::abs(detinv)<1e-80)
        return false;
    detinv = 1/detinv;

    x[0] = detinv*(mat[1][1]*rhs[0]-mat[0][1]*rhs[1]);
    x[1] = detinv*(mat[0][0]*rhs[1]-mat[1][0]*rhs[0]);

    // x[0] is the distance, x[1] is the intersection point 
    // in local coordinates on the segment
    if (x[1]<0 || x[1] > 1)
        return false;

    distance    = x[0];
    targetLocal = x[1];

    return true;

}


// ////////////////////////////////////////////////////////
//   Explicit template instantiations.
//   If you need more, you can add them here.
// ////////////////////////////////////////////////////////

template class PSURFACE_EXPORT NormalProjector<float>;
template class PSURFACE_EXPORT NormalProjector<double>;

#include <psurface/NormalProjector.h>
#include <psurface/PSurfaceFactory.h>
#include <psurface/DirectionFunction.h>

#include <psurface/StaticVector.h>
#include <psurface/StaticMatrix.h>

#include <psurface/NodeBundle.h>

#ifdef PSURFACE_STANDALONE
#include "TargetSurface.h"
#else
#include "hxsurface/Surface.h"
#endif

#include <stdexcept>
#include <vector>


template <class ctype>
void NormalProjector<ctype>::project(const Surface* targetSurface,
                                     const DirectionFunction<3,ctype>* domainDirection,
                                     const DirectionFunction<3,ctype>* targetDirection)
{
    const double eps = 0.0001;

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
        
    // Insert the vertices of the contact boundary as nodes on the intermediate manifold
    // /////////////////////////////////////////////////////////////////////////////////////

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

        for (size_t j=0; j<psurface_->getNumTriangles(); j++) {

            const StaticVector<ctype,3>& p0 = psurface_->vertices(psurface_->triangles(j).vertices[0]);
            const StaticVector<ctype,3>& p1 = psurface_->vertices(psurface_->triangles(j).vertices[1]);
            const StaticVector<ctype,3>& p2 = psurface_->vertices(psurface_->triangles(j).vertices[2]);

            const StaticVector<ctype,3>& n0 = domainNormals[psurface_->triangles(j).vertices[0]];
            const StaticVector<ctype,3>& n1 = domainNormals[psurface_->triangles(j).vertices[1]];
            const StaticVector<ctype,3>& n2 = domainNormals[psurface_->triangles(j).vertices[2]];

            StaticVector<ctype,3> x; // the unknown...

            // magic to use a McVec3f as the argument
            StaticVector<ctype,3> targetVertex;
            for (int k=0; k<3; k++)
                targetVertex[k] = surf->points[i][k];

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

                if (segment.dot(targetNormals[i]) > -eps
                    && segment.dot(baseNormal) > 0
                    && distance > 1e-10) {
                    continue;
                }

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
    // Insert the edges
    // ////////////////////////////////////////////////////////////

    for (int i=0; i<targetSurface->triangles.size(); i++) {

        for (int j=0; j<3; j++) {
            
            int from = targetSurface->triangles[i].points[j];
            int to   = targetSurface->triangles[i].points[(j+1)%3];

            // If (from, to) is an inner edge we pass it twice, but want to
            // test it only once.  That's before the ||.  The second clause
            // tests for boundary edges
            if (from < to || containsEdge(targetSurface, from, to)==1) {
                
                if (edgeCanBeInserted(domainNormals, from, to, projectedTo))
                    insertEdge(factory, domainNormals, from, to, projectedTo);
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
int NormalProjector<ctype>::containsEdge(const Surface* surface, int from, int to) const
{
    int counter=0;
    for (size_t i=0; i<surface->triangles.size(); i++)
        for (int j=0; j<3; j++)
	    if ((surface->triangles[i].points[j]==from && surface->triangles[i].points[(j+1)%3]==to) ||
                (surface->triangles[i].points[j]==to   && surface->triangles[i].points[(j+1)%3]==from))
                counter++;

    return counter;
}

template <class ctype>
void NormalProjector<ctype>::insertEdge(PSurfaceFactory<2,ctype>& factory,
                                        const std::vector<StaticVector<ctype,3> >& normals,
                                 int from, int to, 
                                 const std::vector<NodeBundle>& projectedTo)
{
    int enteringEdge=-1;
    NodeBundle curr = projectedTo[from];

    // parameter value for the edge to be inserted
    ctype lambda = 0;

    while (curr!=projectedTo[to]) {

        // Connect two nodes that are on the same triangle
        if (onSameTriangle(curr, projectedTo[to])) {
            
            // Get the common triangle
            std::vector<int> commonTris = getCommonTris(curr, projectedTo[to]);
            for (int i=0; i<commonTris.size(); i++) {
                psurface_->triangles(commonTris[i]).addEdge(curr.triToIdx(commonTris[i]), 
                                                      projectedTo[to].triToIdx(commonTris[i]));
            }
            
            break;
        }

        try {
            
            typename Node<ctype>::NodeType currType = psurface_->nodes(curr[0]).type;
            switch (currType) {
            case Node<ctype>::TOUCHING_NODE:
                
                insertEdgeFromTouchingNode(factory, normals, from, to, lambda, projectedTo, curr, enteringEdge);
                break;
                
            case Node<ctype>::GHOST_NODE:
            case Node<ctype>::CORNER_NODE:
                insertEdgeFromCornerNode(factory, normals, from, to, lambda, projectedTo, curr, enteringEdge);
                break;
                
            case Node<ctype>::INTERSECTION_NODE:
                insertEdgeFromIntersectionNode(factory, normals, from, to, lambda, projectedTo, curr, enteringEdge);
                break;
                
            case Node<ctype>::INTERIOR_NODE:
                insertEdgeFromInteriorNode(factory, normals, from, to, lambda, projectedTo, curr, enteringEdge);
                break;
                
            }
            
        } catch (std::runtime_error e) {
            std::cout << "Exception caught!" << std::endl;
            return;
        }
    }
    
}


template <class ctype>
void NormalProjector<ctype>::insertEdgeFromInteriorNode(PSurfaceFactory<2,ctype>& factory,
                                                        const std::vector<StaticVector<ctype,3> >& normals,
                                                 int from, int to, ctype &lambda,
                                                 const std::vector<NodeBundle>& projectedTo,
                                                 NodeBundle& curr, int& enteringEdge)
{
    int i;
    ctype eps = 1e-5;
    int cT = curr[0].tri;

    // loop over the three edges of the current triangle (except for the entering edge) and
    // check whether the paramPolyEdge leaves the triangle via this edge
    for (i=0; i<3; i++) {
            
        if (i==enteringEdge)
            continue;
            
        StaticVector<ctype,3> x;
        int p = psurface_->triangles(curr[0].tri).vertices[i];
        int q = psurface_->triangles(curr[0].tri).vertices[(i+1)%3];

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
                throw(std::runtime_error("[FromInteriorNode] Error: the normal projection is not continuous!"));
            }

            int corner = -1;
            if (mu<eps) 
                corner = i;
            else if (mu>1-eps)
                corner = (i+1)%3;

            if (corner==-1) {

                // get neighboring triangle
                int neighboringTri = psurface_->getNeighboringTriangle(curr[0].tri, i);
                    
                // if no neighboring triangle --> error
                if (neighboringTri==-1) {
                        
                    printf("[FromInteriorNode] Warning: Normal images leaves domain surface!\n");
                    curr = projectedTo[to];
                    return;
                        
                }
                    
                // add intersection nodes on both sides
                    
                // the domain position of the new intersection node on the next triangle
                // better: using getEdge()
                int e = psurface_->triangles(neighboringTri).getCorner(q);
                    
                // the domain position of the new intersection node on this triangle
                StaticVector<ctype,2> dom1((i==0)*(1-mu) + (i==2)*mu, (i==0)*mu + (i==1)*(1-mu));
                StaticVector<ctype,2> dom2((e==0)*mu + (e==2)*(1-mu), (e==0)*(1-mu) + (e==1)*mu);
                    
                StaticVector<ctype,3> image;

                // hack: within Amira this is assignment from a McVec3f
                for (int j=0; j<3; j++)
                    image[j] = surf->points[from][j] + newLambda*(surf->points[to][j]-surf->points[from][j]);

                StaticVector<ctype,2> dom1Ctype(dom1[0], dom1[1]);
                StaticVector<ctype,2> dom2Ctype(dom2[0], dom2[1]);
                NodeBundle newNodePair  = factory.addIntersectionNodePair(curr[0].tri, neighboringTri,
                                                                          dom1Ctype, dom2Ctype, i, e, image);
                    
                // insert new parameter edge
                psurface_->triangles(curr[0].tri).addEdge(curr[0].idx, newNodePair[0].idx);

                //
                curr.resize(1);
                curr[0] = newNodePair[1];
                lambda  = newLambda;
                enteringEdge = e;
                
            } else {
                enteringEdge = curr[0].tri;
                psurface_->triangles(curr[0].tri).addEdge(curr[0].idx, getCornerNode(psurface_->triangles(cT), corner));
                curr = psurface_->getNodeBundleAtVertex(psurface_->triangles(cT).vertices[corner]);
                lambda  = newLambda;
            }
            break;
        }
            
    }
    if (i==3) {
        printf("No intersection found (in insertEdgeFromInteriorNode)!\n");
        curr = projectedTo[to];
        assert(false);
        return;
    }
        
}


template <class ctype>
void NormalProjector<ctype>::insertEdgeFromIntersectionNode(PSurfaceFactory<2,ctype>& factory,
                                                            const std::vector<StaticVector<ctype,3> >& normals,
                                                 int from, int to, ctype &lambda,
                                                 const std::vector<NodeBundle>& projectedTo,
                                                 NodeBundle& curr, int& enteringEdge)
{
    int i;
    ctype eps = 1e-5;
    int cTIdx = curr[0].tri;

    // loop over the three edges of the current triangle (except for the entering edge) and
    // check whether the paramPolyEdge leaves the triangle via this edge
    //int currentEdge = cT.nodes[curr[0].idx].getDomainEdge();
    //printf("currentEdge: %d\n", currentEdge);

    for (i=0; i<3; i++) {
            
        if (i==enteringEdge)
            continue;
            
        StaticVector<ctype,3> x;
        int p = psurface_->triangles(curr[0].tri).vertices[i];
        int q = psurface_->triangles(curr[0].tri).vertices[(i+1)%3];

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
                
            if (newLambda < lambda)
                throw(std::runtime_error("[FromIntersectionNode] Error: the normal projection is not continuous!"));

            int corner = -1;
            if (mu<eps) 
                corner = i;
            else if (mu>1-eps)
                corner = (i+1)%3;

            if (corner==-1) {
                // get neighboring triangle
                int neighboringTri = psurface_->getNeighboringTriangle(curr[0].tri, i);
                    
                // if no neighboring triangle --> error
                if (neighboringTri==-1) {
                        
                    std::cout << "[FromIntersectionNode] Warning: Normal images leaves domain surface!" << std::endl;
                    curr = projectedTo[to];
                    //assert(false);
                    return;
                        
                }
                    
                // add intersection nodes on both sides
                    
                // the domain position of the new intersection node on the next triangle
                // better: using getEdge()
                int e = psurface_->triangles(neighboringTri).getCorner(q);
                    
                // the domain position of the new intersection node on this triangle
                StaticVector<ctype,2> dom1((i==0)*(1-mu) + (i==2)*mu, (i==0)*mu + (i==1)*(1-mu));
                StaticVector<ctype,2> dom2((e==0)*mu + (e==2)*(1-mu), (e==0)*(1-mu) + (e==1)*mu);
                    
                StaticVector<ctype,3> image;

                // trick: within Amira this is assignment from a McVec3f
                for (int j=0; j<3; j++)
                    image[j] = surf->points[from][j] + newLambda*(surf->points[to][j]-surf->points[from][j]);
                    
                NodeBundle newNodePair = factory.addIntersectionNodePair(curr[0].tri, neighboringTri,
                                                                         dom1, dom2, i, e, image);
                    
                // insert new parameter edge
                psurface_->triangles(curr[0].tri).addEdge(curr[0].idx, newNodePair[0].idx);

                //
                curr.resize(1);
                curr[0] = newNodePair[1];
                lambda  = newLambda;
                enteringEdge = e;
                
            } else {
                enteringEdge = curr[0].tri;
                psurface_->triangles(curr[0].tri).addEdge(curr[0].idx, getCornerNode(psurface_->triangles(cTIdx), corner));
                curr = psurface_->getNodeBundleAtVertex(psurface_->triangles(cTIdx).vertices[corner]);
                lambda = newLambda;
            }
            break;
        }
            
    }
    if (i==3) {
        printf("No intersection found (in insertEdgeFromIntersectionNode)!\n");
        psurface_->triangles(curr[0].tri).nodes.erase(psurface_->triangles(curr[0].tri).nodes.begin() + curr[0].idx);
        curr = projectedTo[to];
        assert(false);
        return;
    }
        
}


template <class ctype>
void NormalProjector<ctype>::insertEdgeFromTouchingNode(PSurfaceFactory<2,ctype>& factory,
                                                        const std::vector<StaticVector<ctype,3> >& normals,
                                                 int from, int to, ctype &lambda,
                                                 const std::vector<NodeBundle>& projectedTo,
                                                 NodeBundle& curr, int& enteringTri)
{
    const Surface* surf = psurface_->surface;

    // The other end of the edge is *not* on this triangle
    for (int i=0; i<curr.size(); i++) {

        // I don't think the following if-clause is ever true
        assert(enteringTri==-1);
        if (curr[i].tri == enteringTri)
            continue;

        DomainTriangle<ctype>& cT = psurface_->triangles(curr[i].tri);
        int currentEdge = cT.nodes[curr[i].idx].getDomainEdge();
        
        for (int j=0; j<3; j++) {
            
            if (j==currentEdge)
                continue;
            
            StaticVector<ctype,3> x;
            int p = cT.vertices[j];
            int q = cT.vertices[(j+1)%3];

            StaticVector<ctype,3> targetFrom, targetTo;
            for (int k=0; k<3; k++) {
                targetFrom[k] = surf->points[from][k];
                targetTo[k]   = surf->points[to][k];
            }

            if (edgeIntersectsNormalFan(targetFrom, targetTo,
                                        psurface_->vertices(p), psurface_->vertices(q),
                                        normals[p], normals[q], x)) {
                
                if (x[1] < lambda) {
                    throw(std::runtime_error("[FromTouchingNode] Error: the normal projection is not continuous!"));
                }
                
                lambda = x[1];
                const ctype& mu = x[0];

                int corner = -1;
                if (x[0]<0.00001) 
                    corner = j;
                else if (x[0]>0.99999)
                    corner = (j+1)%3;

                if (corner==-1) {
                    // parameter polyedge is leaving basegrid triangle through an edge
                    
                    // get neighboring triangle
                    int neighboringTri = psurface_->getNeighboringTriangle(curr[i].tri, j);

                    // if no neighboring triangle --> error
                    if (neighboringTri==-1) {
                        
                        printf("[FromTouchingNode] Error: Normal images leaves intermediate surface!\n");
                        abort();
                        
                    }
                    // add intersection nodes on both sides
                    
                    // better: using getEdge()
                    int e = psurface_->triangles(neighboringTri).getCorner(q);

                    // the domain position of the new intersection node on this triangle
                    StaticVector<ctype,2> dom1((j==0)*(1-mu) + (j==2)*mu, (j==0)*mu + (j==1)*(1-mu));
                    StaticVector<ctype,2> dom2((e==0)*mu + (e==2)*(1-mu), (e==0)*(1-mu) + (e==1)*mu);
                    
                    StaticVector<ctype,3> image;

                    for (int k=0; k<3; k++)
                        image[k] = surf->points[from][k] + lambda*(surf->points[to][k]-surf->points[from][k]);
                    
                    NodeBundle newNodePair = factory.addIntersectionNodePair(curr[i].tri, neighboringTri,
                                                                             dom1, dom2, j, e, image);
                    
                    // insert new parameter edge
                    psurface_->triangles(curr[i].tri).addEdge(curr[i].idx, newNodePair[0].idx);

                    curr.resize(1);
                    curr[0] = newNodePair[1];
                    enteringTri = e;
                    
                } else if (corner== ((currentEdge+2)%3)) {

                    // parameter polyedge leaves BG triangle through the opposite vertex
                    psurface_->triangles(curr[i].tri).addEdge(curr[i].idx, getCornerNode(cT, corner));
                    enteringTri = curr[i].tri;
                    curr = psurface_->getNodeBundleAtVertex(cT.vertices[corner]);

                } else {
                    // parameter polyedge leaves BG triangle through an adjacent vertex
                    NodeBundle target = psurface_->getNodeBundleAtVertex(cT.vertices[corner]);
                    for (int k=0; k<curr.size(); k++)
                        psurface_->triangles(curr[k].tri).addEdge(curr[k].idx, target.triToIdx(curr[k].tri));

                    curr = target;
                    
                }
                
                return;
                
            }
            
        }
        
    }

    printf("No Intersection found (in insertEdgeFromTouchingNode)!\n");
    curr = projectedTo[to];
    assert(false);
    return;
}


template <class ctype>
void NormalProjector<ctype>::insertEdgeFromCornerNode(PSurfaceFactory<2,ctype>& factory,
                                                      const std::vector<StaticVector<ctype,3> >& normals, int from, int to, ctype &lambda,
                                                 const std::vector<NodeBundle>& projectedTo,
                                                 NodeBundle& curr, int& enteringEdge)
{
    const Surface* surf = psurface_->surface;

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
            
            if (x[1] < lambda) {
                continue;
            }
            
            lambda = x[1];
            const ctype& mu = x[0];
            int corner = -1;
            if (x[0]<0.00001) 
                corner = (thisCorner+1)%3;
            else if (x[0]>0.99999)
                corner = (thisCorner+2)%3;
            
            if (corner==-1) {
                // parameter polyedge is leaving basegrid triangle 
                // through the opposite edge
                
                // get neighboring triangle
                int neighboringTri = psurface_->getNeighboringTriangle(cT, oppEdge);
                // if no neighboring triangle --> error
                if (neighboringTri==-1) {
                    
                    throw(std::runtime_error("[FromCornerNode] Error: Normal images leaves intermediate surface!"));
                    
                }
                // add intersection nodes on both sides
                
                // better: using getEdge()
                int e = psurface_->triangles(neighboringTri).getCorner(q);

                // the domain position of the new intersection node on this triangle
                StaticVector<ctype,2> dom1((oppEdge==0)*(1-mu) + (oppEdge==2)*mu, (oppEdge==0)*mu + (oppEdge==1)*(1-mu));
                StaticVector<ctype,2> dom2((e==0)*mu + (e==2)*(1-mu), (e==0)*(1-mu) + (e==1)*mu);
                
                StaticVector<ctype,3> image;

                for (int j=0; j<3; j++)
                    image[j] = surf->points[from][j] + lambda*(surf->points[to][j]-surf->points[from][j]);
                
                NodeBundle newNodePair = factory.addIntersectionNodePair(cT, neighboringTri,
                                                                         dom1, dom2, oppEdge, e, image);
                
                // insert new parameter edge
                psurface_->triangles(cT).addEdge(curr[i].idx, newNodePair[0].idx);

                curr.resize(1);
                curr[0] = newNodePair[1];
                enteringEdge = e;
                
            } else {

                // edge leaves through a corner
                NodeBundle target = psurface_->getNodeBundleAtVertex(psurface_->triangles(cT).vertices[corner]);
                std::vector<int> commonTris = getCommonTris(curr, target);
                for (i=0; i<commonTris.size(); i++) {
                    psurface_->triangles(commonTris[i]).addEdge(curr.triToIdx(commonTris[i]), 
                                                          target.triToIdx(commonTris[i]));
                }
                
                curr = target;
                
            }
            
            return;
            
        }
        
    }

    printf("no intersection found (in insertEdgeFromCornerNode)!\n");
    curr = projectedTo[to];
    assert(false);
}



template <class ctype>
bool NormalProjector<ctype>::edgeCanBeInserted(const std::vector<StaticVector<ctype,3> >& normals,
                                        int from, int to, 
                                        const std::vector<NodeBundle>& projectedTo)
{
    if (projectedTo[from].size() == 0 || projectedTo[to].size() == 0)
        return false;

    int enteringEdge=-1;
    NodeBundle curr = projectedTo[from];

    // If the two nodes are on the same triangle it is surely possible to enter the edge
    if (onSameTriangle(curr, projectedTo[to]))
        return true;

    typename Node<ctype>::NodeType currType = psurface_->nodes(curr[0]).type;
    int currTri     = curr[0].tri;

    // parameter value for the edge to be inserted
    ctype lambda = 0;

    while (true) {

        // If the two nodes are on the same triangle it is surely possible to enter the edge
        if (onSameTriangle(currTri, projectedTo[to])
            || ((currType==Node<ctype>::GHOST_NODE || currType==Node<ctype>::CORNER_NODE)
                && onSameTriangle(curr, projectedTo[to])))
            return true;

        switch (currType) {
        case Node<ctype>::TOUCHING_NODE:

            if (!testInsertEdgeFromTouchingNode(normals, 
                                                from, to, 
                                                lambda,
                                                curr, currType, currTri, enteringEdge))
                return false;

            break;
            
        case Node<ctype>::GHOST_NODE:
        case Node<ctype>::CORNER_NODE:

            if (!testInsertEdgeFromCornerNode(normals, 
                                              from, to, 
                                              lambda,
                                              curr, currType, currTri, enteringEdge))
                return false;
            break;

        case Node<ctype>::INTERSECTION_NODE:
            if (!testInsertEdgeFromIntersectionNode(normals, from, to, lambda, curr, currType, currTri, enteringEdge))
                return false;
            break;
            
        case Node<ctype>::INTERIOR_NODE:
            if (!testInsertEdgeFromInteriorNode(normals, from, to, lambda, curr, currType, currTri, enteringEdge))
                return false;
            break;
            
        default:
            std::cout << "ERROR: unknown node type found!" << std::endl;
            abort();
        }
            
    }
    
    std::cout << "should not occur" << std::endl;
    return true;
}



template <class ctype>
bool NormalProjector<ctype>::testInsertEdgeFromInteriorNode(const std::vector<StaticVector<ctype,3> >& normals,
                                                     int from, int to, ctype &lambda,
                                                     NodeBundle& curr,
                                                            typename Node<ctype>::NodeType& currType, int& currTri,
                                                     int& enteringEdge)
{
    // loop over the three edges of the current triangle (except for the entering edge) and
    // check whether the paramPolyEdge leaves the triangle via this edge
    ctype eps = 1e-5;

    for (int i=0; i<3; i++) {
            
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

                return true;
   
            }

        }
            
    }

    printf("No intersection found (in testInsertEdgeFromInteriorNode)!\n");
    return false;
}



template <class ctype>
bool NormalProjector<ctype>::testInsertEdgeFromIntersectionNode(const std::vector<StaticVector<ctype,3> >& normals,
                                                         int from, int to, ctype &lambda,
                                                                NodeBundle& curr,
                                                         typename Node<ctype>::NodeType& currType, int& currTri,
                                                         int& enteringEdge)
{
    // loop over the three edges of the current triangle (except for the entering edge) and
    // check whether the paramPolyEdge leaves the triangle via this edge
    ctype eps = 1e-5;

    for (int i=0; i<3; i++) {
            
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
                
                return true;

            }

        }
            
    }

    printf("No intersection found (in testInsertEdgeFromIntersectionNode)!\n");
    return false;
}



template <class ctype>
bool NormalProjector<ctype>::testInsertEdgeFromTouchingNode(const std::vector<StaticVector<ctype,3> >& normals,
                                                     int from, int to, ctype &lambda,
                                                     NodeBundle& curr,
                                                     typename Node<ctype>::NodeType& currType, int& currTri,
                                                     int& enteringEdge)
{
    const Surface* surf = psurface_->surface;
    ctype eps = 1e-5;

    // The other end of the edge is *not* on this triangle
    for (int i=0; i<curr.size(); i++) {

        const DomainTriangle<ctype>& cT = psurface_->triangles(curr[i].tri);
        int currentEdge = cT.nodes[curr[i].idx].getDomainEdge();
        
        for (int j=0; j<3; j++) {
            
            if (j==currentEdge)
                continue;
        
            StaticVector<ctype,3> x;
            int p = cT.vertices[j];
            int q = cT.vertices[(j+1)%3];

            StaticVector<ctype,3> targetFrom, targetTo;
            for (int k=0; k<3; k++) {
                targetFrom[k] = surf->points[from][k];
                targetTo[k]   = surf->points[to][k];
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
                    currTri  = neighboringTri;
                    lambda   = newLambda;
                    enteringEdge = e;
                    
                    return true;
                
                } else {
                    // parameter polyedge is leaving base grid triangle through a ghost node

                    // get all ghost nodes for the base grid vertex
                    int vertex = psurface_->triangles(curr[i].tri).vertices[corner];
                    curr = psurface_->getNodeBundleAtVertex(vertex); 
                    currType = Node<ctype>::GHOST_NODE;
                    assert(currType==type(curr));
                    lambda = newLambda;

                    return true;

                }
                
            }
            
        }

    }
    
    printf("No intersection found (in testInsertEdgeFromTouchingNode)!\n");
    assert(false);
    return false;

}



template <class ctype>
bool NormalProjector<ctype>::testInsertEdgeFromCornerNode(const std::vector<StaticVector<ctype,3> >& normals, 
                                                   int from, int to, ctype &lambda,
                                                   NodeBundle& curr, 
                                                   typename Node<ctype>::NodeType& currType, int& currTri,
                                                   int& leavingEdge)
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
                
                return true;
                
            }
            
        }
        
    }

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
    // Fix some initial value
    x.assign(1.0);

    for (int i=0; i<10; i++) {

        // compute Newton correction
        StaticVector<ctype,3> Fxk = x[0]*(p0-p2) + x[1]*(p1-p2) + x[2]*x[0]*(n0-n2) + x[2]*x[1]*(n1-n2) + x[2]*n2 + p2;// - target;
        Fxk[0] -= target[0];
        Fxk[1] -= target[1];
        Fxk[2] -= target[2];

        StaticMatrix<ctype,3> FPrimexk(p0 - p2 + x[2]*(n0-n2),
                         p1 - p2 + x[2]*(n1-n2),
                         x[0]*(n0-n2) + x[1]*(n1-n2) + n2);

        StaticMatrix<ctype,3> FPrimexkInv = FPrimexk.inverse();

        StaticVector<ctype,3> newtonCorrection; // = (-1) * FPrimexk.inverse() * Fxk;
        
        FPrimexkInv.multMatrixVec(-Fxk, newtonCorrection);

        x += newtonCorrection;

    }

    if (x[0]>=-eps && x[1]>=-eps && (x[0]+x[1] <=1+eps)){
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

    // Fix some initial value
    // sometimes it only works when the initial value is an intersection...
    x[0] = x[1] = 0.5;
    x[2] = 1;
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

    if (x[0]>=0 && x[0]<=1 && x[1]>=0 && x[1]<=1 && newtonCorrection.length()<1e-4){
        return true;
    } 
           
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

template class NormalProjector<float>;
template class NormalProjector<double>;

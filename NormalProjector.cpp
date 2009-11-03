#include <psurface/ContactBoundary.h>
#include <psurface/NormalProjector.h>

#include <psurface/StaticVector.h>
#include <psurface/StaticMatrix.h>

#include <psurface/NodeBundle.h>

#include <stdexcept>
#include <vector>


template <class ctype>
void NormalProjector<ctype>::handleSide(PSurface<2,ctype>* par, const ContactBoundary& contactPatch,
                                 void (*obsDirections)(const double* pos, double* dir))
{
    int i, j;
    const double eps = 0.0001;
    const Surface* surf = par->surface;

    // //////////////////////////
    // compute normals
    // //////////////////////////
    int nPoints = par->getNumVertices();
    int nTriangles = par->getNumTriangles();
        
    std::vector<StaticVector<double,3> > normals(nPoints);
    normals.assign(nPoints, StaticVector<double,3>(0.0));
    
    std::vector<unsigned char> nTriPerVertex(nPoints);

    nTriPerVertex.assign(nTriPerVertex.size(), 0);

    if (obsDirections) {

        for (i=0; i<nPoints; i++) {
            double pos[3];
            double dir[3];

            for (int j=0; j<3; j++)
                pos[j] = par->vertices(i)[j];
            (*obsDirections)(pos, dir);
            for (int j=0; j<3; j++)
                normals[i][j] = dir[j];

        }

    } else {

        for (i=0; i<nTriangles; i++) {
            
            int p0 = par->triangles(i).vertices[0];
            int p1 = par->triangles(i).vertices[1];
            int p2 = par->triangles(i).vertices[2];
            
            StaticVector<ctype,3> a_ = par->vertices(p1) - par->vertices(p0);
            StaticVector<ctype,3> b_ = par->vertices(p2) - par->vertices(p0);
            
            StaticVector<double,3> a(a_[0], a_[1], a_[2]);
            StaticVector<double,3> b(b_[0], b_[1], b_[2]);
            StaticVector<double,3> triNormal = a.cross(b);
            triNormal.normalize();
            
            normals[p0] += triNormal;
            normals[p1] += triNormal;
            normals[p2] += triNormal;
            
            nTriPerVertex[p0]++;
            nTriPerVertex[p1]++;
            nTriPerVertex[p2]++;
        }

        for (i=0; i<nPoints; i++)
            if (nTriPerVertex[i])
                normals[i].normalize();

    }
    
    // /////////////////////////////////////////////////////////////
    //   Compute the vertex normals of the target side
    // /////////////////////////////////////////////////////////////
    int nTargetPoints = contactPatch.surf->points.size();
    int nTargetTriangles = contactPatch.triIdx.size();
        
    targetNormals.assign(nTargetPoints, StaticVector<double,3>(0.0,0.0,0.0));
    std::vector<bool> hasTargetNormal;
    hasTargetNormal.assign(nTargetPoints, false);

    for (i=0; i<nTargetTriangles; i++) {
        
        int p0 = contactPatch.triangles(i).points[0];
        int p1 = contactPatch.triangles(i).points[1];
        int p2 = contactPatch.triangles(i).points[2];
        
        StaticVector<ctype,3> a_, b_;

        for (int j=0; j<3; j++) {
            a_[j] = contactPatch.surf->points[p1][j] - contactPatch.surf->points[p0][j];
            b_[j] = contactPatch.surf->points[p2][j] - contactPatch.surf->points[p0][j];
        }        

        StaticVector<double,3> a(a_[0], a_[1], a_[2]);
        StaticVector<double,3> b(b_[0], b_[1], b_[2]);
        StaticVector<double,3> triNormal = a.cross(b);
        triNormal.normalize();
        
        targetNormals[p0] += triNormal;
        targetNormals[p1] += triNormal;
        targetNormals[p2] += triNormal;
             
        hasTargetNormal[p0] = true;
        hasTargetNormal[p1] = true;
        hasTargetNormal[p2] = true;

    }
    
    for (i=0; i<contactPatch.vertices.size(); i++)
        if (hasTargetNormal[contactPatch.vertices[i]])
            targetNormals[contactPatch.vertices[i]].normalize();

    // /////////////////////////////////////////////////////////////////////////////////////
    // Insert the vertices of the contact boundary as nodes on the intermediate manifold
    // /////////////////////////////////////////////////////////////////////////////////////

    // This array stores the preimages of each vertex in the target surface
    std::vector<NodeBundle> projectedTo(surf->points.size());

    // This bitfield marks whether base grid vertices already have a
    // corresponding image
    std::vector<bool> vertexHasBeenHandled(par->getNumVertices(), false);
    
    // Loop over the vertices of the target surface
    for (i=0; i<contactPatch.vertices.size(); i++) {

        StaticVector<double,2> bestDPos;
        int bestTri = -1;
        double bestDist = std::numeric_limits<double>::max();

        for (j=0; j<par->getNumTriangles(); j++) {

            const StaticVector<ctype,3>& p0 = par->vertices(par->triangles(j).vertices[0]);
            const StaticVector<ctype,3>& p1 = par->vertices(par->triangles(j).vertices[1]);
            const StaticVector<ctype,3>& p2 = par->vertices(par->triangles(j).vertices[2]);

            const StaticVector<double,3>& n0 = normals[par->triangles(j).vertices[0]];
            const StaticVector<double,3>& n1 = normals[par->triangles(j).vertices[1]];
            const StaticVector<double,3>& n2 = normals[par->triangles(j).vertices[2]];

            StaticVector<double,3> x; // the unknown...

            if (computeInverseNormalProjection(p0, p1, p2, n0, n1, n2, 
                                               // magic to use a McVec3f as the argument
                                               *(StaticVector<ctype,3>*)&surf->points[contactPatch.vertices[i]][0], 
                                               x)) {

                // We want that the line from the domain surface to its projection
                // approaches the target surface from the front side, i.e., it should
                // not pass through the body represented by the target surface.
                // We do a simplified test by comparing the connecting segment
                // with the normal at the target surface and the normal at the
                // domain surface
                StaticVector<ctype,3> base       = p0*x[0] + p1*x[1] + (1-x[0]-x[1])*p2;
                StaticVector<double,3> baseNormal = n0*x[0] + n1*x[1] + (1-x[0]-x[1])*n2;
                StaticVector<double,3> segment(surf->points[contactPatch.vertices[i]][0] - base[0],
                                surf->points[contactPatch.vertices[i]][1] - base[1],
                                surf->points[contactPatch.vertices[i]][2] - base[2]);
                
                double distance = segment.length() * segment.length();

                if (segment.dot(targetNormals[contactPatch.vertices[i]]) > -eps
                    && segment.dot(baseNormal) > 0
                    && distance > 1e-10) {
                    continue;
                }

                // There may be several inverse orthogonal projections.
                // We want the shortest one.

                if (distance < bestDist) {

                    bestDist = distance;
                    bestDPos = StaticVector<double,2>(x[0], x[1]);
                    bestTri  = j;

                }
                
            }

        }
        
        if (bestTri != -1) {

            // determine which type of node to add
            typename Node<ctype>::NodeType newType = Node<ctype>::INTERIOR_NODE;
            int dir = -1;
            double mu;

            // if the normal projection hits a base grid vertex, this is the vertex
            int v = -1;

            if (bestDPos[0] < eps) {
                dir = 1;
                mu = 1-bestDPos[1];
                if (bestDPos[1] < eps) {
                    newType = Node<ctype>::CORNER_NODE;
                    v       = par->triangles(bestTri).vertices[2];
                } else if (bestDPos[1] > 1-eps) {
                    newType = Node<ctype>::CORNER_NODE;
                    v       = par->triangles(bestTri).vertices[1];
                } else {
                    newType = Node<ctype>::TOUCHING_NODE;
                }
            } else if (bestDPos[1] < eps) {
                dir = 2;
                mu = bestDPos[0];
                if (bestDPos[0] < eps) {
                    newType = Node<ctype>::CORNER_NODE;
                    v       = par->triangles(bestTri).vertices[2];
                } else if (bestDPos[0] > 1-eps) {
                    newType = Node<ctype>::CORNER_NODE;
                    v       = par->triangles(bestTri).vertices[0];
                } else {
                    newType = Node<ctype>::TOUCHING_NODE;
                }
            } else if (1-bestDPos[0]-bestDPos[1] < eps) {
                dir = 0;
                mu = 1-bestDPos[0];
                if (bestDPos[1] < eps) {
                    newType = Node<ctype>::CORNER_NODE;
                    v       = par->triangles(bestTri).vertices[0];
                } else if (bestDPos[1] > 1-eps) {
                    newType = Node<ctype>::CORNER_NODE;
                    v       = par->triangles(bestTri).vertices[1];
                } else {
                    newType = Node<ctype>::TOUCHING_NODE;
                }
            }
                    
            StaticVector<ctype,2> bestDPosCtype(bestDPos[0], bestDPos[1]);

            if (newType==Node<ctype>::TOUCHING_NODE) {

                // find the other triangle, if there is one
                int neighboringTri = par->getNeighboringTriangle(bestTri, dir);

                if (neighboringTri == -1) {
                    NodeIdx newNodeNumber = par->addTouchingNode(bestTri, bestDPosCtype, dir, contactPatch.vertices[i]);
                    projectedTo[contactPatch.vertices[i]].resize(1);
                    projectedTo[contactPatch.vertices[i]][0].setValue(bestTri, newNodeNumber);
                } else {
                    // find domain pos on other triangle
                    int commonEdge = par->triangles(bestTri).getCommonEdgeWith(par->triangles(neighboringTri));
                    int dir2 = par->triangles(neighboringTri).getEdge(commonEdge);
                    StaticVector<ctype,2> dP2((dir2==0)*(mu) + (dir2==2)*(1-mu), (dir2==0)*(1-mu) + (dir2==1)*(mu));

                    // insert touching node pair
                    NodeIdx newNodeNumber = par->addTouchingNodePair(bestTri, neighboringTri, bestDPosCtype, dP2, 
                                                                     dir, dir2, contactPatch.vertices[i]);
                    projectedTo[contactPatch.vertices[i]].resize(2);
                    projectedTo[contactPatch.vertices[i]][0].setValue(bestTri, newNodeNumber);
                    projectedTo[contactPatch.vertices[i]][1].setValue(neighboringTri, par->triangles(neighboringTri).nodes.size()-1);
                }
                
            } else if (newType==Node<ctype>::CORNER_NODE) {

                addCornerNodeBundle(par, v, contactPatch.vertices[i]);
                vertexHasBeenHandled[v] = true;
                std::vector<int> neighboringTris = par->getTrianglesPerVertex(v);
                projectedTo[contactPatch.vertices[i]].resize(neighboringTris.size());
                for (j=0; j<neighboringTris.size(); j++) {
                    projectedTo[contactPatch.vertices[i]][j].setValue(neighboringTris[j],
                                                                      par->triangles(neighboringTris[j]).nodes.size()-1);
                }

            } else {

                NodeIdx newNodeNumber = par->addInteriorNode(bestTri, bestDPosCtype, contactPatch.vertices[i]);
            
                projectedTo[contactPatch.vertices[i]].resize(1);
                projectedTo[contactPatch.vertices[i]][0].setValue(bestTri, newNodeNumber);
            
            }
            
        }
        
    }

#if 0
    // ////////////////////////////////////////////////////////////
    // Compute the reduced contact boundary, i.e., the image of the
    // intermediate surface under the normal projection

    ContactBoundary reducedContactPatch(surf);
    reducedContactPatch.vertices = contactPatch.vertices;

    for (i=0; i<contactPatch.triIdx.size(); i++) {

        if (projectedTo[contactPatch.triangles(i).points[0]].size() &&
            projectedTo[contactPatch.triangles(i).points[1]].size() &&
            projectedTo[contactPatch.triangles(i).points[2]].size()) {

            reducedContactPatch.triIdx.append(contactPatch.triIdx[i]);

        }

    }
#endif
    // ///////////////////////////////////////////////////////////////////
    //   Place ghost nodes at the vertices of the intermediate manifold
    // ///////////////////////////////////////////////////////////////////
    for (i=0; i<par->getNumVertices(); i++) {

        // Has the vertex been hit by the projection of a target vertex already?
        if (vertexHasBeenHandled[i])
            continue;

        StaticVector<double,2> bestDPos;
        int bestTri = -1;
        double bestDist = std::numeric_limits<double>::max();

        const StaticVector<ctype,3>& basePointCtype = par->vertices(i);
        const StaticVector<double,3> basePoint(basePointCtype[0], basePointCtype[1], basePointCtype[2]);
        const StaticVector<double,3>& normal    = normals[i];

        //for (j=0; j<reducedContactPatch.triIdx.size(); j++) {
        for (j=0; j<contactPatch.triIdx.size(); j++) {

            StaticVector<double,2> domainPos;
            double dist;

            // magic to be able to take a reference of a McVec3f when compiled within Amira
            const StaticVector<ctype,3>& p0 = *(StaticVector<ctype,3>*)&surf->points[contactPatch.triangles(j).points[0]][0];
            const StaticVector<ctype,3>& p1 = *(StaticVector<ctype,3>*)&surf->points[contactPatch.triangles(j).points[1]][0];
            const StaticVector<ctype,3>& p2 = *(StaticVector<ctype,3>*)&surf->points[contactPatch.triangles(j).points[2]][0];

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
            insertGhostNodeAtVertex(par, i, contactPatch.triIdx[bestTri], bestDPos);

    }

    // ////////////////////////////////////////////////////////////
    // Insert the edges
    // ////////////////////////////////////////////////////////////

    for (i=0; i<contactPatch.triIdx.size(); i++) {

        for (j=0; j<3; j++) {
            
            int from = contactPatch.triangles(i).points[j];
            int to   = contactPatch.triangles(i).points[(j+1)%3];

            // If (from, to) is an inner edge we pass it twice, but want to
            // test it only once.  That's before the ||.  The second clause
            // tests for boundary edges
            if (from < to || contactPatch.containsEdge(from, to)==1) {
                
                if (edgeCanBeInserted(par, normals, from, to, projectedTo))
                    insertEdge(par, normals, from, to, projectedTo);
                else {
                    //std::cout << "Skipping edge (" << from << ", " << to << ") ..." << std::endl;
                }
            } 
                
        }
        
    }

    setupEdgePointArrays(par);
}


template <class ctype>
void NormalProjector<ctype>::insertEdge(PSurface<2,ctype>* par, 
                                 const std::vector<StaticVector<double,3> >& normals,
                                 int from, int to, 
                                 const std::vector<NodeBundle>& projectedTo)
{
    int enteringEdge=-1;
    NodeBundle curr = projectedTo[from];

    // parameter value for the edge to be inserted
    double lambda = 0;

    while (curr!=projectedTo[to]) {

        // Connect two nodes that are on the same triangle
        if (onSameTriangle(curr, projectedTo[to])) {
            
            // Get the common triangle
            std::vector<int> commonTris = getCommonTris(curr, projectedTo[to]);
            for (int i=0; i<commonTris.size(); i++) {
                par->triangles(commonTris[i]).addEdge(curr.triToIdx(commonTris[i]), 
                                                      projectedTo[to].triToIdx(commonTris[i]));
            }
            
            break;
        }

        try {
            
            typename Node<ctype>::NodeType currType = par->nodes(curr[0]).type;
            switch (currType) {
            case Node<ctype>::TOUCHING_NODE:
                
                insertEdgeFromTouchingNode(par, normals, from, to, lambda, projectedTo, curr, enteringEdge);
                break;
                
            case Node<ctype>::GHOST_NODE:
            case Node<ctype>::CORNER_NODE:
                insertEdgeFromCornerNode(par, normals, from, to, lambda, projectedTo, curr, enteringEdge);
                break;
                
            case Node<ctype>::INTERSECTION_NODE:
                insertEdgeFromIntersectionNode(par, normals, from, to, lambda, projectedTo, curr, enteringEdge);
                break;
                
            case Node<ctype>::INTERIOR_NODE:
                insertEdgeFromInteriorNode(par, normals, from, to, lambda, projectedTo, curr, enteringEdge);
                break;
                
            }
            
        } catch (std::runtime_error e) {
            std::cout << "Exception caught!" << std::endl;
            return;
        }
    }
    
}


template <class ctype>
void NormalProjector<ctype>::insertEdgeFromInteriorNode(PSurface<2,ctype>* par, 
                                                 const std::vector<StaticVector<double,3> >& normals,
                                                 int from, int to, double &lambda,
                                                 const std::vector<NodeBundle>& projectedTo,
                                                 NodeBundle& curr, int& enteringEdge)
{
    int i;
    int cT = curr[0].tri;

    // loop over the three edges of the current triangle (except for the entering edge) and
    // check whether the paramPolyEdge leaves the triangle via this edge
    for (i=0; i<3; i++) {
            
        if (i==enteringEdge)
            continue;
            
        StaticVector<double,3> x;
        int p = par->triangles(curr[0].tri).vertices[i];
        int q = par->triangles(curr[0].tri).vertices[(i+1)%3];

        const Surface* surf = par->surface;


        if (edgeIntersectsNormalFan(*(StaticVector<ctype,3>*)&surf->points[from][0], 
                                    *(StaticVector<ctype,3>*)&surf->points[to][0],
                                    par->vertices(p), par->vertices(q),
                                    normals[p], normals[q], x)) {

            const double& newLambda = x[1];
            const double& mu        = x[0];
                
            if (newLambda < lambda) {
                throw(std::runtime_error("[FromInteriorNode] Error: the normal projection is not continuous!"));
            }

            int corner = -1;
            if (mu<0.00001) 
                corner = i;
            else if (mu>9.9999)
                corner = (i+1)%3;

            if (corner==-1) {

                // get neighboring triangle
                int neighboringTri = par->getNeighboringTriangle(curr[0].tri, i);
                    
                // if no neighboring triangle --> error
                if (neighboringTri==-1) {
                        
                    printf("[FromInteriorNode] Warning: Normal images leaves domain surface!\n");
                    curr = projectedTo[to];
                    return;
                        
                }
                    
                // add intersection nodes on both sides
                    
                // the domain position of the new intersection node on the next triangle
                // better: using getEdge()
                int e = par->triangles(neighboringTri).getCorner(q);
                    
                // the domain position of the new intersection node on this triangle
                StaticVector<double,2> dom1((i==0)*(1-mu) + (i==2)*mu, (i==0)*mu + (i==1)*(1-mu));
                StaticVector<double,2> dom2((e==0)*mu + (e==2)*(1-mu), (e==0)*(1-mu) + (e==1)*mu);
                    
                StaticVector<ctype,3> image;

                // hack: within Amira this is assignment from a McVec3f
                for (int j=0; j<3; j++)
                    image[j] = surf->points[from][j] + newLambda*(surf->points[to][j]-surf->points[from][j]);

                StaticVector<ctype,2> dom1Ctype(dom1[0], dom1[1]);
                StaticVector<ctype,2> dom2Ctype(dom2[0], dom2[1]);
                NodeIdx newNodeIn  = par->addIntersectionNodePair(curr[0].tri, neighboringTri,
                                                                  dom1Ctype, dom2Ctype, i, e, image);
                NodeIdx newNodeOut = par->triangles(neighboringTri).nodes.size()-1;
                    
                // insert new parameter edge
                par->triangles(curr[0].tri).addEdge(curr[0].idx, newNodeIn);

                //
                curr.resize(1);
                curr[0].setValue(neighboringTri, newNodeOut);
                lambda   = newLambda;
                enteringEdge = e;
                
            } else {
                enteringEdge = curr[0].tri;
                par->triangles(curr[0].tri).addEdge(curr[0].idx, getCornerNode(par->triangles(cT), corner));
                curr = par->getNodeBundleAtVertex(par->triangles(cT).vertices[corner]);
            }
            break;
        }
            
    }
    if (i==3) {
#ifndef NDEBUG
        printf("No intersection found!\n");
#endif
        curr = projectedTo[to];
        return;
        assert(false);
    }
        
}


template <class ctype>
void NormalProjector<ctype>::insertEdgeFromIntersectionNode(PSurface<2,ctype>* par, 
                                                 const std::vector<StaticVector<double,3> >& normals,
                                                 int from, int to, double &lambda,
                                                 const std::vector<NodeBundle>& projectedTo,
                                                 NodeBundle& curr, int& enteringEdge)
{
    int i;
    int cTIdx = curr[0].tri;

    // loop over the three edges of the current triangle (except for the entering edge) and
    // check whether the paramPolyEdge leaves the triangle via this edge
    //int currentEdge = cT.nodes[curr[0].idx].getDomainEdge();
    //printf("currentEdge: %d\n", currentEdge);

    for (i=0; i<3; i++) {
            
        if (i==enteringEdge)
            continue;
            
        StaticVector<double,3> x;
        int p = par->triangles(curr[0].tri).vertices[i];
        int q = par->triangles(curr[0].tri).vertices[(i+1)%3];

        const Surface* surf = par->surface;
        if (edgeIntersectsNormalFan(*(StaticVector<ctype,3>*)&surf->points[from][0], 
                                    *(StaticVector<ctype,3>*)&surf->points[to][0],
                                    par->vertices(p), par->vertices(q),
                                    normals[p], normals[q], x)) {

            const double& newLambda = x[1];
            const double& mu        = x[0];
                
            if (newLambda < lambda)
                throw(std::runtime_error("[FromIntersectionNode] Error: the normal projection is not continuous!"));

            int corner = -1;
            if (mu<0.00001) 
                corner = i;
            else if (mu>9.9999)
                corner = (i+1)%3;

            if (corner==-1) {
                // get neighboring triangle
                int neighboringTri = par->getNeighboringTriangle(curr[0].tri, i);
                    
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
                int e = par->triangles(neighboringTri).getCorner(q);
                    
                // the domain position of the new intersection node on this triangle
                StaticVector<ctype,2> dom1((i==0)*(1-mu) + (i==2)*mu, (i==0)*mu + (i==1)*(1-mu));
                StaticVector<ctype,2> dom2((e==0)*mu + (e==2)*(1-mu), (e==0)*(1-mu) + (e==1)*mu);
                    
                StaticVector<ctype,3> image;

                // trick: within Amira this is assignment from a McVec3f
                for (int j=0; j<3; j++)
                    image[j] = surf->points[from][j] + newLambda*(surf->points[to][j]-surf->points[from][j]);
                    
                NodeIdx newNodeIn  = par->addIntersectionNodePair(curr[0].tri, neighboringTri,
                                                                  dom1, dom2, i, e, image);
                NodeIdx newNodeOut = par->triangles(neighboringTri).nodes.size()-1;
                    
                // insert new parameter edge
                par->triangles(curr[0].tri).addEdge(curr[0].idx, newNodeIn);

                //
                curr.resize(1);
                curr[0].setValue(neighboringTri, newNodeOut);
                lambda   = newLambda;
                enteringEdge = e;
                
            } else {
                enteringEdge = curr[0].tri;
                par->triangles(curr[0].tri).addEdge(curr[0].idx, getCornerNode(par->triangles(cTIdx), corner));
                curr = par->getNodeBundleAtVertex(par->triangles(cTIdx).vertices[corner]);
            }
            break;
        }
            
    }
    if (i==3) {
#ifndef NDEBUG
        printf("No intersection found!\n");
#endif
        par->triangles(curr[0].tri).nodes.erase(par->triangles(curr[0].tri).nodes.begin() + curr[0].idx);
        curr = projectedTo[to];
        return;
        assert(false);
    }
        
}


template <class ctype>
void NormalProjector<ctype>::insertEdgeFromTouchingNode(PSurface<2,ctype>* par,
                                                 const std::vector<StaticVector<double,3> >& normals,
                                                 int from, int to, double &lambda,
                                                 const std::vector<NodeBundle>& projectedTo,
                                                 NodeBundle& curr, int& enteringTri)
{
    const Surface* surf = par->surface;

    // The other end of the edge is *not* on this triangle
    for (int i=0; i<curr.size(); i++) {

        // I don't think the following if-clause is ever true
        assert(enteringTri==-1);
        if (curr[i].tri == enteringTri)
            continue;

        DomainTriangle<ctype>& cT = par->triangles(curr[i].tri);
        int currentEdge = cT.nodes[curr[i].idx].getDomainEdge();
        
        for (int j=0; j<3; j++) {
            
            if (j==currentEdge)
                continue;
            
            StaticVector<double,3> x;
            int p = cT.vertices[j];
            int q = cT.vertices[(j+1)%3];

            if (edgeIntersectsNormalFan(*(StaticVector<ctype,3>*)&surf->points[from][0],
                                        *(StaticVector<ctype,3>*)&surf->points[to][0],
                                        par->vertices(p), par->vertices(q),
                                        normals[p], normals[q], x)) {
                
                if (x[1] < lambda) {
                    throw(std::runtime_error("[FromTouchingNode] Error: the normal projection is not continuous!"));
                }
                
                lambda = x[1];
                const double& mu = x[0];

                int corner = -1;
                if (x[0]<0.00001) 
                    corner = j;
                else if (x[0]>0.99999)
                    corner = (j+1)%3;

                if (corner==-1) {
                    // parameter polyedge is leaving basegrid triangle through an edge
                    
                    // get neighboring triangle
                    int neighboringTri = par->getNeighboringTriangle(curr[i].tri, j);

                    // if no neighboring triangle --> error
                    if (neighboringTri==-1) {
                        
                        printf("[FromTouchingNode] Error: Normal images leaves intermediate surface!\n");
                        abort();
                        
                    }
                    // add intersection nodes on both sides
                    
                    // better: using getEdge()
                    int e = par->triangles(neighboringTri).getCorner(q);

                    // the domain position of the new intersection node on this triangle
                    StaticVector<ctype,2> dom1((j==0)*(1-mu) + (j==2)*mu, (j==0)*mu + (j==1)*(1-mu));
                    StaticVector<ctype,2> dom2((e==0)*mu + (e==2)*(1-mu), (e==0)*(1-mu) + (e==1)*mu);
                    
                    StaticVector<ctype,3> image;

                    for (int k=0; k<3; k++)
                        image[k] = surf->points[from][k] + lambda*(surf->points[to][k]-surf->points[from][k]);
                    
                    NodeIdx newNodeIn  = par->addIntersectionNodePair(curr[i].tri, neighboringTri,
                                                                        dom1, dom2, j, e, image);
                    NodeIdx newNodeOut = par->triangles(neighboringTri).nodes.size()-1;
                    
                    // insert new parameter edge
                    par->triangles(curr[i].tri).addEdge(curr[i].idx, newNodeIn);

                    curr.resize(1);
                    curr[0].setValue(neighboringTri, newNodeOut);
                    enteringTri = e;
                    
                } else if (corner== ((currentEdge+2)%3)) {

                    // parameter polyedge leaves BG triangle through the opposite vertex
                    par->triangles(curr[i].tri).addEdge(curr[i].idx, getCornerNode(cT, corner));
                    enteringTri = curr[i].tri;
                    curr = par->getNodeBundleAtVertex(cT.vertices[corner]);

                } else {
                    // parameter polyedge leaves BG triangle through an adjacent vertex
                    NodeBundle target = par->getNodeBundleAtVertex(cT.vertices[corner]);
                    for (int k=0; k<curr.size(); k++)
                        par->triangles(curr[k].tri).addEdge(curr[k].idx, target.triToIdx(curr[k].tri));

                    curr = target;
                    
                }
                
                return;
                
            }
            
        }
        
    }

#ifndef NDEBUG
    printf("No Intersection found!\n");
#endif
    curr = projectedTo[to];
}


template <class ctype>
void NormalProjector<ctype>::insertEdgeFromCornerNode(PSurface<2,ctype>* par,
                                               const std::vector<StaticVector<double,3> >& normals, int from, int to, double &lambda,
                                                 const std::vector<NodeBundle>& projectedTo,
                                                 NodeBundle& curr, int& enteringEdge)
{
    const Surface* surf = par->surface;

    // The other end of the edge is *not* on this triangle
    for (int i=0; i<curr.size(); i++) {
        
        int cT = curr[i].tri;

        // it is called enteringEdge, but the incoming value is the entering*Tri*!
        if (cT == enteringEdge)
            continue;

        int thisCorner = par->triangles(cT).nodes[curr[i].idx].getCorner();
        int oppEdge = (thisCorner+1)%3;

        StaticVector<double,3> x;
        int p = par->triangles(cT).vertices[(thisCorner+1)%3];
        int q = par->triangles(cT).vertices[(thisCorner+2)%3];

        if (edgeIntersectsNormalFan(*(StaticVector<ctype,3>*)&surf->points[from][0],
                                    *(StaticVector<ctype,3>*)&surf->points[to][0],
                                    par->vertices(p), par->vertices(q),
                                    normals[p], normals[q], x)) {
            
            if (x[1] < lambda) {
                continue;
            }
            
            lambda = x[1];
            const double& mu = x[0];
            int corner = -1;
            if (x[0]<0.00001) 
                corner = (thisCorner+1)%3;
            else if (x[0]>0.99999)
                corner = (thisCorner+2)%3;
            
            if (corner==-1) {
                // parameter polyedge is leaving basegrid triangle 
                // through the opposite edge
                
                // get neighboring triangle
                int neighboringTri = par->getNeighboringTriangle(cT, oppEdge);
                // if no neighboring triangle --> error
                if (neighboringTri==-1) {
                    
                    throw(std::runtime_error("[FromCornerNode] Error: Normal images leaves intermediate surface!"));
                    
                }
                // add intersection nodes on both sides
                
                // better: using getEdge()
                int e = par->triangles(neighboringTri).getCorner(q);

                // the domain position of the new intersection node on this triangle
                StaticVector<ctype,2> dom1((oppEdge==0)*(1-mu) + (oppEdge==2)*mu, (oppEdge==0)*mu + (oppEdge==1)*(1-mu));
                StaticVector<ctype,2> dom2((e==0)*mu + (e==2)*(1-mu), (e==0)*(1-mu) + (e==1)*mu);
                
                StaticVector<ctype,3> image;

                for (int j=0; j<3; j++)
                    image[j] = surf->points[from][j] + lambda*(surf->points[to][j]-surf->points[from][j]);
                
                NodeIdx newNodeIn  = par->addIntersectionNodePair(cT, neighboringTri,
                                                                  dom1, dom2, oppEdge, e, image);
                NodeIdx newNodeOut = par->triangles(neighboringTri).nodes.size()-1;
                
                // insert new parameter edge
                par->triangles(cT).addEdge(curr[i].idx, newNodeIn);

                curr.resize(1);
                curr[0].setValue(neighboringTri, newNodeOut);
                enteringEdge = e;
                
            } else {

                // edge leaves through a corner
                NodeBundle target = par->getNodeBundleAtVertex(par->triangles(cT).vertices[corner]);
                std::vector<int> commonTris = getCommonTris(curr, target);
                for (i=0; i<commonTris.size(); i++) {
                    par->triangles(commonTris[i]).addEdge(curr.triToIdx(commonTris[i]), 
                                                          target.triToIdx(commonTris[i]));
                }
                
                curr = target;
                
            }
            
            return;
            
        }
        
    }

#ifndef NDEBUG
    printf("no intersection found!\n");
#endif
    curr = projectedTo[to];
}



template <class ctype>
bool NormalProjector<ctype>::edgeCanBeInserted(const PSurface<2,ctype>* par, 
                                        const std::vector<StaticVector<double,3> >& normals,
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

    typename Node<ctype>::NodeType currType = par->nodes(curr[0]).type;
    int currTri     = curr[0].tri;

    // parameter value for the edge to be inserted
    double lambda = 0;

    while (true) {

        // If the two nodes are on the same triangle it is surely possible to enter the edge
        if (onSameTriangle(currTri, projectedTo[to]))
            return true;

        switch (currType) {
        case Node<ctype>::TOUCHING_NODE:

            if (!testInsertEdgeFromTouchingNode(par, normals, 
                                                from, to, 
                                                lambda, projectedTo, 
                                                curr, currType, currTri, enteringEdge))
                return false;

            break;
            
        case Node<ctype>::GHOST_NODE:
        case Node<ctype>::CORNER_NODE:

            if (!testInsertEdgeFromCornerNode(par, normals, 
                                              from, to, 
                                              lambda, projectedTo, 
                                              curr, currType, currTri, enteringEdge))
                return false;
            break;

        case Node<ctype>::INTERSECTION_NODE:
            if (!testInsertEdgeFromIntersectionNode(par, normals, from, to, lambda, projectedTo, currType, currTri, enteringEdge))
                return false;
            break;
            
        case Node<ctype>::INTERIOR_NODE:
            if (!testInsertEdgeFromInteriorNode(par, normals, from, to, lambda, projectedTo, currType, currTri, enteringEdge))
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
bool NormalProjector<ctype>::testInsertEdgeFromInteriorNode(const PSurface<2,ctype>* par, 
                                                     const std::vector<StaticVector<double,3> >& normals,
                                                     int from, int to, double &lambda,
                                                     const std::vector<NodeBundle>& projectedTo,
                                                            typename Node<ctype>::NodeType& currType, int& currTri,
                                                     int& enteringEdge)
{
    // loop over the three edges of the current triangle (except for the entering edge) and
    // check whether the paramPolyEdge leaves the triangle via this edge
    int i;
    for (i=0; i<3; i++) {
            
        if (i==enteringEdge)
            continue;
            
        StaticVector<double,3> x;
        int p = par->triangles(currTri).vertices[i];
        int q = par->triangles(currTri).vertices[(i+1)%3];

        const Surface* surf = par->surface;


        if (edgeIntersectsNormalFan(*(StaticVector<ctype,3>*)&surf->points[from][0],
                                    *(StaticVector<ctype,3>*)&surf->points[to][0],
                                    par->vertices(p), par->vertices(q),
                                    normals[p], normals[q], x)) {

            const double& newLambda = x[1];
            const double& mu        = x[0];
                
            if (newLambda < lambda) {
                // Error: the normal projection is not continuous!
                return false;
            }

            int corner = -1;
            if (mu<0.00001) 
                corner = i;
            else if (mu>9.9999)
                corner = (i+1)%3;

            if (corner==-1) {
                // get neighboring triangle
                int neighboringTri = par->getNeighboringTriangle(currTri, i);
                    
                // if no neighboring triangle --> error
                if (neighboringTri==-1) {
                    // Error: Normal images leaves domain surface!
                    return false;
                }
                    
                // add intersection nodes on both sides
                    
                // the domain position of the new intersection node on the next triangle
                // better: using getEdge()
                int e = par->triangles(neighboringTri).getCorner(q);
                    
                currType = Node<ctype>::INTERSECTION_NODE;
                currTri  = neighboringTri;
                lambda   = newLambda;
                enteringEdge = e;
                
            } else {
                assert(false);
            }
            break;
        }
            
    }
    if (i==3) {
        printf("No intersection found!\n");
        return false;
    }
        
    return true;
}



template <class ctype>
bool NormalProjector<ctype>::testInsertEdgeFromIntersectionNode(const PSurface<2,ctype>* par, 
                                                         const std::vector<StaticVector<double,3> >& normals,
                                                         int from, int to, double &lambda,
                                                         const std::vector<NodeBundle>& projectedTo,
                                                         typename Node<ctype>::NodeType& currType, int& currTri,
                                                         int& enteringEdge)
{
    // loop over the three edges of the current triangle (except for the entering edge) and
    // check whether the paramPolyEdge leaves the triangle via this edge
    int i;
    for (i=0; i<3; i++) {
            
        if (i==enteringEdge)
            continue;
            
        StaticVector<double,3> x;
        int p = par->triangles(currTri).vertices[i];
        int q = par->triangles(currTri).vertices[(i+1)%3];
        
        const Surface* surf = par->surface;
        if (edgeIntersectsNormalFan(*(StaticVector<ctype,3>*)&surf->points[from][0],
                                    *(StaticVector<ctype,3>*)&surf->points[to][0],
                                    par->vertices(p), par->vertices(q),
                                    normals[p], normals[q], x)) {

            const double& newLambda = x[1];
            const double& mu        = x[0];
                
            if (newLambda < lambda) {
                // Error: the normal projection is not continuous!
                return false;                    
            }

            int corner = -1;
            if (mu<0.00001) 
                corner = i;
            else if (mu>9.9999)
                corner = (i+1)%3;

            if (corner==-1) {
                // get neighboring triangle
                int neighboringTri = par->getNeighboringTriangle(currTri, i);
                    
                // if no neighboring triangle --> error
                if (neighboringTri==-1) {
                    // Error: Normal images leaves domain surface!
                    return false;
                }
                    
                // add intersection nodes on both sides
                    
                // the domain position of the new intersection node on the next triangle
                // better: using getEdge()
                int e = par->triangles(neighboringTri).getCorner(q);

                currType = Node<ctype>::INTERSECTION_NODE;
                currTri  = neighboringTri;
                lambda   = newLambda;
                enteringEdge = e;
                
            } else {
                assert(false);
            }
            break;
        }
            
    }
    if (i==3) {
        printf("No intersection found!\n");
        return false;
    }
    return true;
}



template <class ctype>
bool NormalProjector<ctype>::testInsertEdgeFromTouchingNode(const PSurface<2,ctype>* par,
                                                     const std::vector<StaticVector<double,3> >& normals,
                                                     int from, int to, double &lambda,
                                                     const std::vector<NodeBundle>& projectedTo,
                                                     NodeBundle& curr,
                                                     typename Node<ctype>::NodeType& currType, int& currTri,
                                                     int& enteringEdge)
{
    const Surface* surf = par->surface;

    // The other end of the edge is *not* on this triangle
    for (int i=0; i<curr.size(); i++) {

        const DomainTriangle<ctype>& cT = par->triangles(curr[i].tri);
        int currentEdge = cT.nodes[curr[i].idx].getDomainEdge();
        
        for (int j=0; j<3; j++) {
            
            if (j==currentEdge)
                continue;
        
            StaticVector<double,3> x;
            int p = cT.vertices[j];
            int q = cT.vertices[(j+1)%3];

            if (edgeIntersectsNormalFan(*(StaticVector<ctype,3>*)&surf->points[from][0],
                                        *(StaticVector<ctype,3>*)&surf->points[to][0],
                                        par->vertices(p), par->vertices(q),
                                        normals[p], normals[q], x)) {
                
                const double& newLambda = x[1];
                
                if (newLambda < lambda) {
                    // Edge insertion not possible: the normal projection is not continuous!
                    return false;                    
                }
                
                int corner = -1;
                if (x[0]<0.00001) 
                    corner = j;
                else if (x[0]>0.99999)
                    corner = (j+1)%3;
                
                if (corner==-1) {
                    // parameter polyedge is leaving basegrid triangle through an edge
                    
                    // get neighboring triangle
                    int neighboringTri = par->getNeighboringTriangle(curr[i].tri, j);

                    // if no neighboring triangle --> error
                    if (neighboringTri==-1)
                        return false;
                    
                    // add intersection nodes on both sides
                
                    // better: using getEdge()
                    int e = par->triangles(neighboringTri).getCorner(q);
                
                    currType = Node<ctype>::INTERSECTION_NODE;
                    currTri  = neighboringTri;
                    lambda   = newLambda;
                    enteringEdge = e;
                    
                    return true;
                
                } else {
                    // parameter polyedge is leaving base grid triangle through a ghost node

                    // get all ghost nodes for the base grid vertex
                    int vertex = par->triangles(currTri).vertices[corner];
                    std::vector<int> neighbors = par->getTrianglesPerVertex(vertex);

                    curr.resize(0);
                    for (int k=0; k<neighbors.size(); k++) {
                        
                        int cornerOnNeighbor = par->triangles(neighbors[k]).getCorner(vertex);

                        /** \todo Linear search: pretty slow */
                        for (int l=0; l<par->triangles(neighbors[k]).nodes.size(); l++) {

                            if (par->triangles(neighbors[k]).nodes[l].isGHOST_NODE()
                                && par->triangles(neighbors[k]).nodes[l].getCorner() == cornerOnNeighbor){

                                curr.push_back(GlobalNodeIdx(neighbors[k], l));
                                break;

                            }
                        
                        }
                        
                    }

                    currType = Node<ctype>::GHOST_NODE;
                    //currTri = ???;
                    //enteringEdge = -1;

                    lambda = newLambda;

                    return true;

                }
                
            }
            
        }

    }
    
    printf("No intersection found!\n");
    return false;

}



template <class ctype>
bool NormalProjector<ctype>::testInsertEdgeFromCornerNode(const PSurface<2,ctype>* par,
                                                   const std::vector<StaticVector<double,3> >& normals, 
                                                   int from, int to, double &lambda,
                                                   const std::vector<NodeBundle>& projectedTo,
                                                   NodeBundle& curr, 
                                                   typename Node<ctype>::NodeType& currType, int& currTri,
                                                   int& leavingEdge)
{
    const Surface* surf = par->surface;

    // The other end of the edge is *not* on this triangle
    for (int i=0; i<curr.size(); i++) {
        
        int cT = curr[i].tri;

        int thisCorner = par->triangles(cT).nodes[curr[i].idx].getCorner();
        int oppEdge = (thisCorner+1)%3;

        StaticVector<double,3> x;
        int p = par->triangles(cT).vertices[(thisCorner+1)%3];
        int q = par->triangles(cT).vertices[(thisCorner+2)%3];

        if (edgeIntersectsNormalFan(*(StaticVector<ctype,3>*)&surf->points[from][0],
                                    *(StaticVector<ctype,3>*)&surf->points[to][0],
                                    par->vertices(p), par->vertices(q),
                                    normals[p], normals[q], x)) {
            
            const double& newLambda = x[1];
            
            if (newLambda < lambda) {
                // Shouldn't this rather be a 'return false' here?
                continue;
            }
            
            int corner = -1;
            if (x[0]<0.00001) 
                corner = (thisCorner+1)%3;
            else if (x[0]>0.99999)
                corner = (thisCorner+2)%3;
            
            if (corner==-1) {
                // parameter polyedge is leaving basegrid triangle 
                // through the opposite edge
                
                // get neighboring triangle
                int neighboringTri = par->getNeighboringTriangle(cT, oppEdge);

                // if no neighboring triangle --> error
                if (neighboringTri==-1)
                    return false;

                // add intersection nodes on both sides
                
                // better: using getEdge()
                int e = par->triangles(neighboringTri).getCorner(q);

                currType = Node<ctype>::INTERSECTION_NODE;
                currTri  = neighboringTri;
                lambda   = newLambda;
                leavingEdge = e;

                return true;

            } else {

                // parameter polyedge is leaving base grid triangle through a ghost node

                // get all ghost nodes for the base grid vertex
                int vertex = par->triangles(currTri).vertices[corner];
                std::vector<int> neighbors = par->getTrianglesPerVertex(vertex);
                
                curr.resize(0);
                for (int k=0; k<neighbors.size(); k++) {
                    
                    int cornerOnNeighbor = par->triangles(neighbors[k]).getCorner(vertex);
                    
                    /** \todo Linear search: pretty slow */
                    for (int l=0; l<par->triangles(neighbors[k]).nodes.size(); l++) {
                        
                        if (par->triangles(neighbors[k]).nodes[l].isGHOST_NODE()
                            && par->triangles(neighbors[k]).nodes[l].getCorner() == cornerOnNeighbor){
                            
                            curr.push_back(GlobalNodeIdx(neighbors[k], l));
                            break;
                            
                        }
                        
                    }
                    
                }
                
                currType = Node<ctype>::GHOST_NODE;
                //currTri = ???;
                //enteringEdge = -1;
                
                lambda = newLambda;
                
                return true;
                
            }
            
        }
        
    }

    printf("no intersection found!\n");
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
    for (int j=0; j<b.size(); j++)
        if (tri==b[j].tri)
            return true;

    return false;
}


template <class ctype>
void NormalProjector<ctype>::insertGhostNodeAtVertex(PSurface<2,ctype>* par, int v, 
                                              int targetTri, const StaticVector<double,2>& localTargetCoords)
{
    std::vector<int> neighbors = par->getTrianglesPerVertex(v);

    for (int i=0; i<neighbors.size(); i++) {

        int corner = par->triangles(neighbors[i]).getCorner(v);
        StaticVector<ctype,2> lTC_Ctype(localTargetCoords[0], localTargetCoords[1]);
        par->addGhostNode(neighbors[i], corner, targetTri, lTC_Ctype);

    }
 
}


template <class ctype>
void NormalProjector<ctype>::addCornerNodeBundle(PSurface<2,ctype>* cS, int v, int nN)
{
    std::vector<int> neighbors = cS->getTrianglesPerVertex(v);

    for (int i=0; i<neighbors.size(); i++) {

        int corner = cS->triangles(neighbors[i]).getCorner(v);
        cS->addCornerNode(neighbors[i], corner, nN);

    }
        
}


template <class ctype>
bool NormalProjector<ctype>::computeInverseNormalProjection(const StaticVector<ctype,3>& p0_f, const StaticVector<ctype,3>& p1_f, const StaticVector<ctype,3>& p2_f,
                                                     const StaticVector<double,3>& n0, const StaticVector<double,3>& n1, const StaticVector<double,3>& n2,
                                                     const StaticVector<ctype,3>& target, StaticVector<double,3>& x)
{
    int i;
    const double eps = 1e-6;
    // Fix some initial value
    x.assign(1.0);

    // transform to double
    StaticVector<double,3> p0(p0_f[0], p0_f[1], p0_f[2]);
    StaticVector<double,3> p1(p1_f[0], p1_f[1], p1_f[2]);
    StaticVector<double,3> p2(p2_f[0], p2_f[1], p2_f[2]);

    for (i=0; i<10; i++) {

        // compute Newton correction
        StaticVector<double,3> Fxk = x[0]*(p0-p2) + x[1]*(p1-p2) + x[2]*x[0]*(n0-n2) + x[2]*x[1]*(n1-n2) + x[2]*n2 + p2;// - target;
        Fxk[0] -= target[0];
        Fxk[1] -= target[1];
        Fxk[2] -= target[2];

        StaticMatrix<double,3> FPrimexk(p0 - p2 + x[2]*(n0-n2),
                         p1 - p2 + x[2]*(n1-n2),
                         x[0]*(n0-n2) + x[1]*(n1-n2) + n2);

        StaticMatrix<double,3> FPrimexkInv = FPrimexk.inverse();

        StaticVector<double,3> newtonCorrection; // = (-1) * FPrimexk.inverse() * Fxk;
        
        FPrimexkInv.multMatrixVec(-Fxk, newtonCorrection);

        x += newtonCorrection;

    }

    if (x[0]>=-eps && x[1]>=-eps && (x[0]+x[1] <=1+eps)){
        return true;
    } 
        
    return false;
}


template <class ctype>
bool NormalProjector<ctype>::edgeIntersectsNormalFan(const StaticVector<ctype,3>& q0_f, const StaticVector<ctype,3>& q1_f,
                                              const StaticVector<ctype,3>& p0_f, const StaticVector<ctype,3>& p1_f,
                                              const StaticVector<double,3>& n0, const StaticVector<double,3>& n1,
                                              StaticVector<double,3>& x)
{
    int i;
    // transform to double values
    StaticVector<double,3> q0(q0_f[0], q0_f[1], q0_f[2]);
    StaticVector<double,3> q1(q1_f[0], q1_f[1], q1_f[2]);
    StaticVector<double,3> p0(p0_f[0], p0_f[1], p0_f[2]);
    StaticVector<double,3> p1(p1_f[0], p1_f[1], p1_f[2]);

    // Fix some initial value
    // sometimes it only works when the initial value is an intersection...
    x[0] = x[1] = 0.5;
    x[2] = 1;
    StaticVector<double,3> newtonCorrection;

    for (i=0; i<30; i++) {

        // compute Newton correction

        StaticVector<double,3> Fxk = p0-q0 + x[0]*(p1-p0) + x[2]*n0 + x[2]*x[0]*(n1-n0) - x[1]*(q1-q0);

        StaticMatrix<double,3> FPrimexk(p1-p0 + x[2]*(n1-n0),
                         q0-q1,
                         n0 + x[0]*(n1-n0));

        StaticMatrix<double,3> FPrimexkInv = FPrimexk.inverse();
        
        FPrimexkInv.multMatrixVec(-Fxk, newtonCorrection);

        x += newtonCorrection;

    }

    if (x[0]>=0 && x[0]<=1 && x[1]>=0 && x[1]<=1 && newtonCorrection.length()<1e-4){
        return true;
    } 
           
    return false;
}


template <class ctype>
bool NormalProjector<ctype>::rayIntersectsTriangle(const StaticVector<double,3>& basePoint, const StaticVector<double,3>& direction,
                                            const StaticVector<ctype,3>& a_, const StaticVector<ctype,3>& b_, const StaticVector<ctype,3>& c_,
                                            StaticVector<double,2>& localCoords, double& normalDist, double eps)
{
    const StaticVector<double,3> &p = basePoint;

    StaticVector<double,3> a(a_[0], a_[1], a_[2]);
    StaticVector<double,3> b(b_[0], b_[1], b_[2]);
    StaticVector<double,3> c(c_[0], c_[1], c_[2]);

    StaticVector<double,3> e1 = b-a;
    StaticVector<double,3> e2 = c-a;
    e1.normalize();
    e2.normalize();
    bool parallel = fabs(StaticMatrix<double,3>(e1, e2, direction).det()) <eps;
        
        // Cramer's rule
        
        if (!parallel){

            double det = StaticMatrix<double,3>(b-a, c-a, direction).det();
            
            // triangle and edge are not parallel
            double nu = StaticMatrix<double,3>(b-a, c-a, p-a).det() / det;

            double lambda = StaticMatrix<double,3>(p-a, c-a, direction).det() / det;
            if (lambda<-eps) return false;

            double mu = StaticMatrix<double,3>(b-a, p-a, direction).det() / det;
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
            double alpha = StaticMatrix<double,3>(b-a, c-a, p-a).det();
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

    for (int i=0; i<cT.nodes.size(); i++)
        if ((cT.nodes[i].isCORNER_NODE() || cT.nodes[i].isGHOST_NODE()) &&
            cT.nodes[i].getCorner()==corner)
            return i;

    return -1;
}


template <class ctype>
int NormalProjector<ctype>::getCommonTri(const NodeBundle& a, const NodeBundle& b)
{
    for (int i=0; i<a.size(); i++)
        for (int j=0; j<b.size(); j++)
            if (a[i].tri==b[j].tri)
                return a[i].tri;

    return -1;
}


template <class ctype>
std::vector<int> NormalProjector<ctype>::getCommonTris(const NodeBundle& a, const NodeBundle& b)
{
    std::vector<int> result;

    for (int i=0; i<a.size(); i++)
        for (int j=0; j<b.size(); j++)
            if (a[i].tri==b[j].tri)
                result.push_back(a[i].tri);

    return result;
}



template <class ctype>
void NormalProjector<ctype>::setupEdgePointArrays(PSurface<2,ctype>* par)
{
    int i, j;

    for (i=0; i<par->getNumTriangles(); i++) {

        DomainTriangle<ctype>& cT = par->triangles(i);

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
                
            double lambda = cN.getDomainEdgeCoord();
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


// ////////////////////////////////////////////////////////
//   Explicit template instantiations.
//   If you need more, you can add them here.
// ////////////////////////////////////////////////////////

template class NormalProjector<float>;
template class NormalProjector<double>;

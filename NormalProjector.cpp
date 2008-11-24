#include "ContactBoundary.h"
#include "NormalProjector.h"

#include <mclib/McVec2d.h>
#include <mclib/McMat3d.h>
#include <mclib/McVec3d.h>

#include "NodeBundle.h"

#include <stdexcept>
#include <vector>


void NormalProjector::handleSide(Parametrization* par, const ContactBoundary& contactPatch,
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
        
    std::vector<McVec3d> normals(nPoints);
    normals.assign(nPoints, McVec3d(0.0,0.0,0.0));
    
    McDArray<unsigned char> nTriPerVertex(nPoints);

    nTriPerVertex.fill(0);

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
            
            McVec3f a_ = par->vertices(p1) - par->vertices(p0);
            McVec3f b_ = par->vertices(p2) - par->vertices(p0);
            
            McVec3d a(a_[0], a_[1], a_[2]);
            McVec3d b(b_[0], b_[1], b_[2]);
            McVec3d triNormal = a.cross(b);
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
    
//     for (i=0; i<normals.size(); i++)
//         printf("normal %d:  %g %g %g\n", i, normals[i][0], normals[i][1], normals[i][2]);
    // /////////////////////////////////////////////////////////////
    //   Compute the vertex normals of the target side
    // /////////////////////////////////////////////////////////////
    int nTargetPoints = contactPatch.surf->points.size();
    int nTargetTriangles = contactPatch.triIdx.size();
        
    targetNormals.assign(nTargetPoints, McVec3d(0.0,0.0,0.0));
    std::vector<bool> hasTargetNormal;
    hasTargetNormal.assign(nTargetPoints, false);

    for (i=0; i<nTargetTriangles; i++) {
        
        int p0 = contactPatch.triangles(i).points[0];
        int p1 = contactPatch.triangles(i).points[1];
        int p2 = contactPatch.triangles(i).points[2];
        
        McVec3f a_ = contactPatch.surf->points[p1] - contactPatch.surf->points[p0];
        McVec3f b_ = contactPatch.surf->points[p2] - contactPatch.surf->points[p0];
        
        McVec3d a(a_[0], a_[1], a_[2]);
        McVec3d b(b_[0], b_[1], b_[2]);
        McVec3d triNormal = a.cross(b);
        triNormal.normalize();
        
        targetNormals[p0] += triNormal;
        targetNormals[p1] += triNormal;
        targetNormals[p2] += triNormal;
             
        hasTargetNormal[p0] = true;
        hasTargetNormal[p1] = true;
        hasTargetNormal[p2] = true;

    }
    
    for (i=0; i<contactPatch.vertices.size(); i++) {
        if (hasTargetNormal[contactPatch.vertices[i]])
            targetNormals[contactPatch.vertices[i]].normalize();
//             printf("targetNormal %d:  %g %g %g\n", i, 
//                    targetNormals[contactPatch.vertices[i]][0],
//                    targetNormals[contactPatch.vertices[i]][1],
//                    targetNormals[contactPatch.vertices[i]][2]);
    }

    // //////////////////////////////////////////
    // Insert the vertices of the contact boundary as nodes on the intermediate manifold

    // This array stores the preimages of each vertex in the target surface
    std::vector<NodeBundle> projectedTo(surf->points.size());

    // This bitfield marks whether base grid vertices already have a
    // corresponding image
    std::vector<bool> vertexHasBeenHandled(par->getNumVertices(), false);
    
    // Loop over the vertices of the target surface
    for (i=0; i<contactPatch.vertices.size(); i++) {

        McVec2d bestDPos;
        TriangleIdx bestTri = -1;
        double bestDist = std::numeric_limits<double>::max();

        for (j=0; j<par->getNumTriangles(); j++) {

            const McVec3f& p0 = par->vertices(par->triangles(j).vertices[0]);
            const McVec3f& p1 = par->vertices(par->triangles(j).vertices[1]);
            const McVec3f& p2 = par->vertices(par->triangles(j).vertices[2]);

            const McVec3d& n0 = normals[par->triangles(j).vertices[0]];
            const McVec3d& n1 = normals[par->triangles(j).vertices[1]];
            const McVec3d& n2 = normals[par->triangles(j).vertices[2]];

            McVec3d x; // the unknown...

            if (computeInverseNormalProjection(p0, p1, p2, n0, n1, n2, 
                                               surf->points[contactPatch.vertices[i]], x)) {

                // We want that the line from the domain surface to its projection
                // approaches the target surface from the front side, i.e., it should
                // not pass through the body represented by the target surface.
                // We do a simplified test by comparing the connecting segment
                // with the normal at the target surface and the normal at the
                // domain surface
                McVec3f base       = p0*x[0] + p1*x[1] + (1-x[0]-x[1])*p2;
                McVec3d baseNormal = n0*x[0] + n1*x[1] + (1-x[0]-x[1])*n2;
                McVec3d segment(surf->points[contactPatch.vertices[i]][0] - base[0],
                                surf->points[contactPatch.vertices[i]][1] - base[1],
                                surf->points[contactPatch.vertices[i]][2] - base[2]);
                
//                 printf("------------------\n");
//                 printf("base:     %g %g %g\n", base[0], base[1], base[2]);
//                 printf("target:   %g %g %g\n", surf->points[contactPatch.vertices[i]][0],
//                        surf->points[contactPatch.vertices[i]][1],
//                        surf->points[contactPatch.vertices[i]][2]);
//                 printf("normal: %g %g %g\n", targetNormals[contactPatch.vertices[i]][0],
//                        targetNormals[contactPatch.vertices[i]][1],
//                        targetNormals[contactPatch.vertices[i]][2]);
//                 printf("scalar product %g\n", segment.dot(targetNormals[contactPatch.vertices[i]]));
                double distance = segment.length2();

                if (segment.dot(targetNormals[contactPatch.vertices[i]]) > -eps
                    && segment.dot(baseNormal) > 0
                    && distance > 1e-10) {
                    //printf("aborting %g %g %g\n", segment.dot(targetNormals[contactPatch.vertices[i]]), 
                    //       segment.dot(baseNormal), distance);
                    continue;
                }

                // There may be several inverse orthogonal projections.
                // We want the shortest one.

                if (distance < bestDist) {

                    bestDist = distance;
                    bestDPos = McVec2d(x[0], x[1]);
                    bestTri  = j;

                }
                
            }

        }
        
        //printf("vertex: %d, bestTri: %d,  bestDPos %g %g\n", i, bestTri, bestDPos.x, bestDPos.y);

        if (bestTri != -1) {

            // determine which type of node to add
            Node::NodeType newType = Node::INTERIOR_NODE;
            int dir = -1;
            double mu;

            // if the normal projection hits a base grid vertex, this is the vertex
            VertexIdx v = -1;

            if (bestDPos.x < eps) {
                dir = 1;
                mu = 1-bestDPos.y;
                if (bestDPos.y < eps) {
                    newType = Node::CORNER_NODE;
                    v       = par->triangles(bestTri).vertices[2];
                } else if (bestDPos.y > 1-eps) {
                    newType = Node::CORNER_NODE;
                    v       = par->triangles(bestTri).vertices[1];
                } else {
                    newType = Node::TOUCHING_NODE;
                }
            } else if (bestDPos.y < eps) {
                dir = 2;
                mu = bestDPos.x;
                if (bestDPos.x < eps) {
                    newType = Node::CORNER_NODE;
                    v       = par->triangles(bestTri).vertices[2];
                } else if (bestDPos.x > 1-eps) {
                    newType = Node::CORNER_NODE;
                    v       = par->triangles(bestTri).vertices[0];
                } else {
                    newType = Node::TOUCHING_NODE;
                }
            } else if (1-bestDPos.x-bestDPos.y < eps) {
                dir = 0;
                mu = 1-bestDPos.x;
                if (bestDPos.y < eps) {
                    newType = Node::CORNER_NODE;
                    v       = par->triangles(bestTri).vertices[0];
                } else if (bestDPos.y > 1-eps) {
                    newType = Node::CORNER_NODE;
                    v       = par->triangles(bestTri).vertices[1];
                } else {
                    newType = Node::TOUCHING_NODE;
                }
            }
                    
            //printf("newType: %d\n", newType);

            McVec2f bestDPosFloat(bestDPos.x, bestDPos.y);

            if (newType==Node::TOUCHING_NODE) {

                // find the other triangle, if there is one
                TriangleIdx neighboringTri = par->getNeighboringTriangle(bestTri, dir);
                //printf("dir: %d, bestTri: %d,  neighboringTri: %d\n", dir, bestTri, neighboringTri);
                if (neighboringTri == -1) {
                    NodeIdx newNodeNumber = par->addTouchingNode(bestTri, bestDPosFloat, dir, contactPatch.vertices[i]);
                    projectedTo[contactPatch.vertices[i]].resize(1);
                    projectedTo[contactPatch.vertices[i]][0].setValue(bestTri, newNodeNumber);
                } else {
                    // find domain pos on other triangle
                    EdgeIdx commonEdge = par->triangles(bestTri).getCommonEdgeWith(par->triangles(neighboringTri));
                    int dir2 = par->triangles(neighboringTri).getEdge(commonEdge);
                    McVec2f dP2((dir2==0)*(mu) + (dir2==2)*(1-mu), (dir2==0)*(1-mu) + (dir2==1)*(mu));

                    // insert touching node pair
                    NodeIdx newNodeNumber = par->addTouchingNodePair(bestTri, neighboringTri, bestDPosFloat, dP2, 
                                                                     dir, dir2, contactPatch.vertices[i]);
                    projectedTo[contactPatch.vertices[i]].resize(2);
                    projectedTo[contactPatch.vertices[i]][0].setValue(bestTri, newNodeNumber);
                    projectedTo[contactPatch.vertices[i]][1].setValue(neighboringTri, par->triangles(neighboringTri).nodes.size()-1);
                }
                
            } else if (newType==Node::CORNER_NODE) {

                addCornerNodeBundle(par, v, contactPatch.vertices[i]);
                vertexHasBeenHandled[v] = true;
                McSmallArray<TriangleIdx, 12> neighboringTris = par->getTrianglesPerVertex(v);
                projectedTo[contactPatch.vertices[i]].resize(neighboringTris.size());
                for (j=0; j<neighboringTris.size(); j++) {
                    projectedTo[contactPatch.vertices[i]][j].setValue(neighboringTris[j],
                                                                      par->triangles(neighboringTris[j]).nodes.size()-1);
                }

            } else {

                NodeIdx newNodeNumber = par->addInteriorNode(bestTri, bestDPosFloat, contactPatch.vertices[i]);
            
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

        McVec2d bestDPos;
        TriangleIdx bestTri = -1;
        double bestDist = std::numeric_limits<double>::max();

        const McVec3f& basePointFloat = par->vertices(i);
        const McVec3d basePoint(basePointFloat.x, basePointFloat.y, basePointFloat.z);
        const McVec3d& normal    = normals[i];

        //for (j=0; j<reducedContactPatch.triIdx.size(); j++) {
        for (j=0; j<contactPatch.triIdx.size(); j++) {

            McVec2d domainPos;
            double dist;

//             const McVec3f& p0 = surf->points[reducedContactPatch.triangles(j).points[0]];
//             const McVec3f& p1 = surf->points[reducedContactPatch.triangles(j).points[1]];
//             const McVec3f& p2 = surf->points[reducedContactPatch.triangles(j).points[2]];
            const McVec3f& p0 = surf->points[contactPatch.triangles(j).points[0]];
            const McVec3f& p1 = surf->points[contactPatch.triangles(j).points[1]];
            const McVec3f& p2 = surf->points[contactPatch.triangles(j).points[2]];

            if (rayIntersectsTriangle(basePoint, normal, p0, p1, p2, domainPos, dist, eps)) {

                if (dist<bestDist) {
                    bestTri = j;
                    bestDPos = domainPos;
                    bestDist = dist;
                }

            }
            
            // ...
            if (bestTri != -1) {
                //printf("Insert ghostNode at (%f %f)\n", bestDPos.x, bestDPos.y);
                //insertGhostNodeAtVertex(par, i, reducedContactPatch.triIdx[bestTri], bestDPos);
                insertGhostNodeAtVertex(par, i, contactPatch.triIdx[bestTri], bestDPos);
                break;
            }

        }

    }
    //return;
    // ////////////////////////////////////////////////////////////
    // Insert the edges
    // ////////////////////////////////////////////////////////////
    //for (i=0; i<reducedContactPatch.triIdx.size(); i++) {
    for (i=0; i<contactPatch.triIdx.size(); i++) {

        for (j=0; j<3; j++) {
            
//             int from = reducedContactPatch.triangles(i).points[j];
//             int to   = reducedContactPatch.triangles(i).points[(j+1)%3];
            int from = contactPatch.triangles(i).points[j];
            int to   = contactPatch.triangles(i).points[(j+1)%3];

            //if (from < to || reducedContactPatch.containsEdge(from, to)==1) {
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


void NormalProjector::insertEdge(Parametrization* par, 
                                 const std::vector<McVec3d>& normals,
                                 int from, int to, 
                                 const std::vector<NodeBundle>& projectedTo)
{
    //printf("------------------------------\n");
    //printf("Entering insertEdge, from = %d, to = %d\n", from, to);

//     Node::NodeType fromType = par->nodes(projectedTo[from][0]).type;
//     Node::NodeType toType   = par->nodes(projectedTo[to][0]).type;
//     printf("FromType: %d,   toType: %d\n", fromType, toType);

    int enteringEdge=-1;
    NodeBundle curr = projectedTo[from];

    // parameter value for the edge to be inserted
    double lambda = 0;

    while (curr!=projectedTo[to]) {

        // Connect two nodes that are on the same triangle
        if (onSameTriangle(curr, projectedTo[to])) {
            
            // Get the common triangle
            McSmallArray<TriangleIdx, 2> commonTris = getCommonTris(curr, projectedTo[to]);
            for (int i=0; i<commonTris.size(); i++) {
                par->triangles(commonTris[i]).addEdge(curr.triToIdx(commonTris[i]), 
                                                      projectedTo[to].triToIdx(commonTris[i]));
            }
            
            break;
        }

        try {
            
            Node::NodeType currType = par->nodes(curr[0]).type;
            switch (currType) {
            case Node::TOUCHING_NODE:
                
                insertEdgeFromTouchingNode(par, normals, from, to, lambda, projectedTo, curr, enteringEdge);
                break;
                
            case Node::GHOST_NODE:
            case Node::CORNER_NODE:
                //printf("** corner\n");
                insertEdgeFromCornerNode(par, normals, from, to, lambda, projectedTo, curr, enteringEdge);
                //printf("-- corner\n");
                break;
                
            case Node::INTERSECTION_NODE:
                insertEdgeFromIntersectionNode(par, normals, from, to, lambda, projectedTo, curr, enteringEdge);
                break;
                
            case Node::INTERIOR_NODE:
                //printf("** interior\n");
                insertEdgeFromInteriorNode(par, normals, from, to, lambda, projectedTo, curr, enteringEdge);
                //printf("-- interior\n");
                break;
                
            }
            
        } catch (std::runtime_error e) {
            std::cout << "Exception caught!" << std::endl;
            return;
        }
    }
    
}

void NormalProjector::insertEdgeFromInteriorNode(Parametrization* par, 
                                                 const std::vector<McVec3d>& normals,
                                                 int from, int to, double &lambda,
                                                 const std::vector<NodeBundle>& projectedTo,
                                                 NodeBundle& curr, int& enteringEdge)
{
    int i;
    TriangleIdx cT = curr[0].tri;

    //curr.print();

    // loop over the three edges of the current triangle (except for the entering edge) and
    // check whether the paramPolyEdge leaves the triangle via this edge
    for (i=0; i<3; i++) {
            
        if (i==enteringEdge)
            continue;
            
        //printf("i = %d, enteringEdge = %d \n", i, enteringEdge);
        McVec3d x;
        VertexIdx p = par->triangles(curr[0].tri).vertices[i];
        VertexIdx q = par->triangles(curr[0].tri).vertices[(i+1)%3];
        //printf("p: %d   q: %d\n", p, q);
        const Surface* surf = par->surface;


        if (edgeIntersectsNormalFan(surf->points[from], surf->points[to],
                                    par->vertices(p), par->vertices(q),
                                    normals[p], normals[q], x)) {

//             printf("Intersection found (%f %f)\n", x[0], x[1]);
//             printf("p: %d, q: %d\n", p, q);
//             printf("i: %d  entEdge: %d\n", i, enteringEdge);
            const double& newLambda = x[1];
            const double& mu        = x[0];
                
            if (newLambda < lambda) {
                //printf("newLambda: %f   lambda: %f\n", newLambda, lambda);
                throw(std::runtime_error("[FromInteriorNode] Error: the normal projection is not continuous!"));
            }

            int corner = -1;
            if (mu<0.00001) 
                corner = i;
            else if (mu>9.9999)
                corner = (i+1)%3;

//             printf("corner: %d\n", corner);
            if (corner==-1) {
                // get neighboring triangle
                TriangleIdx neighboringTri = par->getNeighboringTriangle(curr[0].tri, i);
                    
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
                McVec2d dom1((i==0)*(1-mu) + (i==2)*mu, (i==0)*mu + (i==1)*(1-mu));
                McVec2d dom2((e==0)*mu + (e==2)*(1-mu), (e==0)*(1-mu) + (e==1)*mu);
                    
                McVec3f image = surf->points[from] + newLambda*(surf->points[to]-surf->points[from]);

                McVec2f dom1Float(dom1.x, dom1.y);
                McVec2f dom2Float(dom2.x, dom2.y);
                NodeIdx newNodeIn  = par->addIntersectionNodePair(curr[0].tri, neighboringTri,
                                                                  dom1Float, dom2Float, i, e, image);
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
                
                //printf("the vertex in question %d\n", cSurf->triangles(cT).vertices[corner]);
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

void NormalProjector::insertEdgeFromIntersectionNode(Parametrization* par, 
                                                 const std::vector<McVec3d>& normals,
                                                 int from, int to, double &lambda,
                                                 const std::vector<NodeBundle>& projectedTo,
                                                 NodeBundle& curr, int& enteringEdge)
{
    int i;
    TriangleIdx cTIdx = curr[0].tri;

    //curr.print();

    // loop over the three edges of the current triangle (except for the entering edge) and
    // check whether the paramPolyEdge leaves the triangle via this edge
    //int currentEdge = cT.nodes[curr[0].idx].getDomainEdge();
    //printf("currentEdge: %d\n", currentEdge);

    for (i=0; i<3; i++) {
            
        if (i==enteringEdge)
            continue;
            
        //   printf("i = %d, enteringEdge = %d \n", i, enteringEdge);
        McVec3d x;
        VertexIdx p = par->triangles(curr[0].tri).vertices[i];
        VertexIdx q = par->triangles(curr[0].tri).vertices[(i+1)%3];
        //printf("p: %d   q: %d\n", p, q);
        const Surface* surf = par->surface;
        if (edgeIntersectsNormalFan(surf->points[from], surf->points[to],
                                    par->vertices(p), par->vertices(q),
                                    normals[p], normals[q], x)) {
//             printf("Intersection found (%f %f)\n", x[0], x[1]);
//             printf("p: %d, q: %d\n", p, q);
//             printf("i: %d  entEdge: %d\n", i, enteringEdge);
            const double& newLambda = x[1];
            const double& mu        = x[0];
                
            if (newLambda < lambda) {
                //printf("newLambda: %f   lambda: %f\n", newLambda, lambda);
                throw(std::runtime_error("[FromIntersectionNode] Error: the normal projection is not continuous!"));
                    
            }

            int corner = -1;
            if (mu<0.00001) 
                corner = i;
            else if (mu>9.9999)
                corner = (i+1)%3;

//             printf("corner: %d\n", corner);
            if (corner==-1) {
                // get neighboring triangle
                TriangleIdx neighboringTri = par->getNeighboringTriangle(curr[0].tri, i);
                    
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
                McVec2f dom1((i==0)*(1-mu) + (i==2)*mu, (i==0)*mu + (i==1)*(1-mu));
                McVec2f dom2((e==0)*mu + (e==2)*(1-mu), (e==0)*(1-mu) + (e==1)*mu);
                    
                McVec3f image = surf->points[from] + newLambda*(surf->points[to]-surf->points[from]);
                    
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
                
                //printf("the vertex in question %d\n", cSurf->triangles(cT).vertices[corner]);
            }
            break;
        }
            
    }
    if (i==3) {
#ifndef NDEBUG
        printf("No intersection found!\n");
#endif
        par->triangles(curr[0].tri).nodes.remove(curr[0].idx);
        curr = projectedTo[to];
        return;
        assert(false);
    }
        
}

void NormalProjector::insertEdgeFromTouchingNode(Parametrization* par,
                                                 const std::vector<McVec3d>& normals,
                                                 int from, int to, double &lambda,
                                                 const std::vector<NodeBundle>& projectedTo,
                                                 NodeBundle& curr, int& enteringTri)
{
    //printf("Edge from TouchingNode\n");
    const Surface* surf = par->surface;

//     printf("curr:\n");
//     curr.print();
//     printf("enteringEdge: %d\n", enteringEdge);
    // The other end of the edge is *not* on this triangle
    for (int i=0; i<curr.size(); i++) {

        // I don't think the following if-clause is ever true
        assert(enteringTri==-1);
        if (curr[i].tri == enteringTri)
            continue;

        DomainTriangle& cT = par->triangles(curr[i].tri);
        int currentEdge = cT.nodes[curr[i].idx].getDomainEdge();
        
        for (int j=0; j<3; j++) {
            
            if (j==currentEdge)
                continue;
            
            McVec3d x;
            VertexIdx p = cT.vertices[j];
            VertexIdx q = cT.vertices[(j+1)%3];

            if (edgeIntersectsNormalFan(surf->points[from], surf->points[to],
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
                    TriangleIdx neighboringTri = par->getNeighboringTriangle(curr[i].tri, j);
                    //printf("neighboringTri = %d\n", neighboringTri);
                    // if no neighboring triangle --> error
                    if (neighboringTri==-1) {
                        
                        printf("[FromTouchingNode] Error: Normal images leaves intermediate surface!\n");
                        abort();
                        
                    }
                    // add intersection nodes on both sides
                    
                    // better: using getEdge()
                    int e = par->triangles(neighboringTri).getCorner(q);

                    // the domain position of the new intersection node on this triangle
                    McVec2f dom1((j==0)*(1-mu) + (j==2)*mu, (j==0)*mu + (j==1)*(1-mu));
                    McVec2f dom2((e==0)*mu + (e==2)*(1-mu), (e==0)*(1-mu) + (e==1)*mu);
                    
                    McVec3f image = surf->points[from] + lambda*(surf->points[to]-surf->points[from]);
                    
                    NodeIdx newNodeIn  = par->addIntersectionNodePair(curr[i].tri, neighboringTri,
                                                                        dom1, dom2, j, e, image);
                    NodeIdx newNodeOut = par->triangles(neighboringTri).nodes.size()-1;
                    
                    // insert new parameter edge
                    par->triangles(curr[i].tri).addEdge(curr[i].idx, newNodeIn);

                    curr.resize(1);
                    curr[0].setValue(neighboringTri, newNodeOut);
                    enteringTri = e;
                    
                } else if (corner== ((currentEdge+2)%3)) {
                    //printf("opp Vertex\n");
                    // parameter polyedge leaves BG triangle through the opposite vertex
                    //cSurf->addParEdge(side, curr[i].tri, curr[i].idx, getCornerNode(cT, corner));
                    par->triangles(curr[i].tri).addEdge(curr[i].idx, getCornerNode(cT, corner));
                    enteringTri = curr[i].tri;
                    curr = par->getNodeBundleAtVertex(cT.vertices[corner]);

                } else {
                    // parameter polyedge leaves BG triangle through an adjacent vertex
                    //printf("adj vertex\n");
                    NodeBundle target = par->getNodeBundleAtVertex(cT.vertices[corner]);
                    //target.print();
                    for (int k=0; k<curr.size(); k++)
                        //cSurf->addParEdge(side, curr[k].tri, curr[k].idx, target.triToIdx(curr[k].tri));
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


void NormalProjector::insertEdgeFromCornerNode(Parametrization* par,
                                               const std::vector<McVec3d>& normals, int from, int to, double &lambda,
                                                 const std::vector<NodeBundle>& projectedTo,
                                                 NodeBundle& curr, int& enteringEdge)
{
    int i;

//     printf("Edge from CornerNode\n");
//      printf("enteringTri: %d\n", enteringEdge);
//      printf("curr:\n");
//      curr.print();
//      printf("++++++++++\n");

    const Surface* surf = par->surface;

//     printf("from:\n");
//     projectedTo[from].print();
//     printf("to:\n");
//     projectedTo[to].print();
    
    // The other end of the edge is *not* on this triangle
    for (i=0; i<curr.size(); i++) {
        
        TriangleIdx cT = curr[i].tri;

        // it is called enteringEdge, but the incoming value is the entering*Tri*!
        if (cT == enteringEdge)
            continue;

        int thisCorner = par->triangles(cT).nodes[curr[i].idx].getCorner();
        int oppEdge = (thisCorner+1)%3;
        //printf("Testing Triangle %d\n", cT);
        McVec3d x;
        VertexIdx p = par->triangles(cT).vertices[(thisCorner+1)%3];
        VertexIdx q = par->triangles(cT).vertices[(thisCorner+2)%3];
        //printf("p: %d   q: %d\n", p, q);
        if (edgeIntersectsNormalFan(surf->points[from], surf->points[to],
                                    par->vertices(p), par->vertices(q),
                                    //normals[from], normals[to], x)) {
                                    normals[p], normals[q], x)) {
            
            //printf("newLambda %f   lambda %f\n", x[0], lambda);
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
            
            //printf("Corner! %d\n", corner);
            if (corner==-1) {
                // parameter polyedge is leaving basegrid triangle 
                // through the opposite edge
                
                // get neighboring triangle
                TriangleIdx neighboringTri = par->getNeighboringTriangle(cT, oppEdge);
                //printf("neighboringTri = %d\n", neighboringTri);
                // if no neighboring triangle --> error
                if (neighboringTri==-1) {
                    
                    throw(std::runtime_error("[FromCornerNode] Error: Normal images leaves intermediate surface!"));
                    
                }
                // add intersection nodes on both sides
                
                // better: using getEdge()
                int e = par->triangles(neighboringTri).getCorner(q);

                // the domain position of the new intersection node on this triangle
                McVec2f dom1((oppEdge==0)*(1-mu) + (oppEdge==2)*mu, (oppEdge==0)*mu + (oppEdge==1)*(1-mu));
                McVec2f dom2((e==0)*mu + (e==2)*(1-mu), (e==0)*(1-mu) + (e==1)*mu);
                
                McVec3f image = surf->points[from] + lambda*(surf->points[to]-surf->points[from]);
                
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
                McSmallArray<TriangleIdx, 2> commonTris = getCommonTris(curr, target);
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










bool NormalProjector::edgeCanBeInserted(const Parametrization* par, 
                                        const std::vector<McVec3d>& normals,
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

    Node::NodeType currType = par->nodes(curr[0]).type;
    TriangleIdx currTri     = curr[0].tri;

    // parameter value for the edge to be inserted
    double lambda = 0;

    while (true) {

        // If the two nodes are on the same triangle it is surely possible to enter the edge
        if (onSameTriangle(currTri, projectedTo[to]))
            return true;

        switch (currType) {
        case Node::TOUCHING_NODE:

            if (!testInsertEdgeFromTouchingNode(par, normals, 
                                                from, to, 
                                                lambda, projectedTo, 
                                                curr, currType, currTri, enteringEdge))
                return false;

            break;
            
        case Node::GHOST_NODE:
        case Node::CORNER_NODE:

            if (!testInsertEdgeFromCornerNode(par, normals, 
                                              from, to, 
                                              lambda, projectedTo, 
                                              curr, currType, currTri, enteringEdge))
                return false;
            break;

        case Node::INTERSECTION_NODE:
            if (!testInsertEdgeFromIntersectionNode(par, normals, from, to, lambda, projectedTo, currType, currTri, enteringEdge))
                return false;
            break;
            
        case Node::INTERIOR_NODE:
            if (!testInsertEdgeFromInteriorNode(par, normals, from, to, lambda, projectedTo, currType, currTri, enteringEdge))
                return false;
            break;
            
        default:
            assert(0);
        }
            
    }
    
    std::cout << "should not occur" << std::endl;
    return true;
}

bool NormalProjector::testInsertEdgeFromInteriorNode(const Parametrization* par, 
                                                     const std::vector<McVec3d>& normals,
                                                     int from, int to, double &lambda,
                                                     const std::vector<NodeBundle>& projectedTo,
                                                     Node::NodeType& currType, TriangleIdx& currTri,
                                                     int& enteringEdge)
{
    // loop over the three edges of the current triangle (except for the entering edge) and
    // check whether the paramPolyEdge leaves the triangle via this edge
    int i;
    for (i=0; i<3; i++) {
            
        if (i==enteringEdge)
            continue;
            
        //printf("i = %d, enteringEdge = %d \n", i, enteringEdge);
        McVec3d x;
        VertexIdx p = par->triangles(currTri).vertices[i];
        VertexIdx q = par->triangles(currTri).vertices[(i+1)%3];
        //printf("p: %d   q: %d\n", p, q);
        const Surface* surf = par->surface;


        if (edgeIntersectsNormalFan(surf->points[from], surf->points[to],
                                    par->vertices(p), par->vertices(q),
                                    normals[p], normals[q], x)) {

//             printf("Intersection found (%f %f)\n", x[0], x[1]);
//             printf("p: %d, q: %d\n", p, q);
//             printf("i: %d  entEdge: %d\n", i, enteringEdge);
            const double& newLambda = x[1];
            const double& mu        = x[0];
                
            if (newLambda < lambda) {
                //printf("newLambda: %f   lambda: %f\n", newLambda, lambda);
                // Error: the normal projection is not continuous!
                return false;
            }

            int corner = -1;
            if (mu<0.00001) 
                corner = i;
            else if (mu>9.9999)
                corner = (i+1)%3;

//             printf("corner: %d\n", corner);
            if (corner==-1) {
                // get neighboring triangle
                TriangleIdx neighboringTri = par->getNeighboringTriangle(currTri, i);
                    
                // if no neighboring triangle --> error
                if (neighboringTri==-1) {
                    // Error: Normal images leaves domain surface!
                    return false;
                }
                    
                // add intersection nodes on both sides
                    
                // the domain position of the new intersection node on the next triangle
                // better: using getEdge()
                int e = par->triangles(neighboringTri).getCorner(q);
                    
                currType = Node::INTERSECTION_NODE;
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

bool NormalProjector::testInsertEdgeFromIntersectionNode(const Parametrization* par, 
                                                         const std::vector<McVec3d>& normals,
                                                         int from, int to, double &lambda,
                                                         const std::vector<NodeBundle>& projectedTo,
                                                         Node::NodeType& currType, TriangleIdx& currTri,
                                                         int& enteringEdge)
{
    // loop over the three edges of the current triangle (except for the entering edge) and
    // check whether the paramPolyEdge leaves the triangle via this edge
    //int currentEdge = cT.nodes[curr[0].idx].getDomainEdge();
    //printf("currentEdge: %d\n", currentEdge);

    int i;
    for (i=0; i<3; i++) {
            
        if (i==enteringEdge)
            continue;
            
        McVec3d x;
        VertexIdx p = par->triangles(currTri).vertices[i];
        VertexIdx q = par->triangles(currTri).vertices[(i+1)%3];
        
        const Surface* surf = par->surface;
        if (edgeIntersectsNormalFan(surf->points[from], surf->points[to],
                                    par->vertices(p), par->vertices(q),
                                    normals[p], normals[q], x)) {
//             printf("Intersection found (%f %f)\n", x[0], x[1]);
//             printf("p: %d, q: %d\n", p, q);
//             printf("i: %d  entEdge: %d\n", i, enteringEdge);
            const double& newLambda = x[1];
            const double& mu        = x[0];
                
            if (newLambda < lambda) {
                //printf("newLambda: %f   lambda: %f\n", newLambda, lambda);
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
                TriangleIdx neighboringTri = par->getNeighboringTriangle(currTri, i);
                    
                // if no neighboring triangle --> error
                if (neighboringTri==-1) {
                    // Error: Normal images leaves domain surface!
                    return false;
                }
                    
                // add intersection nodes on both sides
                    
                // the domain position of the new intersection node on the next triangle
                // better: using getEdge()
                int e = par->triangles(neighboringTri).getCorner(q);

                currType = Node::INTERSECTION_NODE;
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


bool NormalProjector::testInsertEdgeFromTouchingNode(const Parametrization* par,
                                                     const std::vector<McVec3d>& normals,
                                                     int from, int to, double &lambda,
                                                     const std::vector<NodeBundle>& projectedTo,
                                                     const NodeBundle& curr,
                                                     Node::NodeType& currType, TriangleIdx& currTri,
                                                     int& enteringEdge)
{
    //printf("Edge from TouchingNode\n");
    const Surface* surf = par->surface;

    // The other end of the edge is *not* on this triangle
    for (int i=0; i<curr.size(); i++) {

        const DomainTriangle& cT = par->triangles(curr[i].tri);
        int currentEdge = cT.nodes[curr[i].idx].getDomainEdge();
        
        for (int j=0; j<3; j++) {
            
            if (j==currentEdge)
                continue;
        
            McVec3d x;
            VertexIdx p = cT.vertices[j];
            VertexIdx q = cT.vertices[(j+1)%3];

            if (edgeIntersectsNormalFan(surf->points[from], surf->points[to],
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
                    TriangleIdx neighboringTri = par->getNeighboringTriangle(curr[i].tri, j);

                    // if no neighboring triangle --> error
                    if (neighboringTri==-1)
                        return false;
                    
                    // add intersection nodes on both sides
                
                    // better: using getEdge()
                    int e = par->triangles(neighboringTri).getCorner(q);
                
                    currType = Node::INTERSECTION_NODE;
                    currTri  = neighboringTri;
                    lambda   = newLambda;
                    enteringEdge = e;
                    
                    return true;
                
                } else {
                    // Shouldn't happen: we have to leave the triangle through an
                    // edge.  If the parameter edge should go through a domain vertex, there would
                    // be a node on that vertex.
                    abort();
                }
                
            }
            
        }

    }
    
    printf("No intersection found!\n");
    return false;

}


bool NormalProjector::testInsertEdgeFromCornerNode(const Parametrization* par,
                                                   const std::vector<McVec3d>& normals, 
                                                   int from, int to, double &lambda,
                                                   const std::vector<NodeBundle>& projectedTo,
                                                   const NodeBundle& curr, 
                                                   Node::NodeType& currType, TriangleIdx& currTri,
                                                   int& enteringEdge)
{
    int i;

    const Surface* surf = par->surface;

    // The other end of the edge is *not* on this triangle
    for (i=0; i<curr.size(); i++) {
        
        TriangleIdx cT = curr[i].tri;

#if 0   // Pointless: there is no entering anything!
        // it is called enteringEdge, but the incoming value is the entering*Tri*!
        if (cT == enteringEdge)
            continue;
#endif

        int thisCorner = par->triangles(cT).nodes[curr[i].idx].getCorner();
        int oppEdge = (thisCorner+1)%3;
        //printf("Testing Triangle %d\n", cT);
        McVec3d x;
        VertexIdx p = par->triangles(cT).vertices[(thisCorner+1)%3];
        VertexIdx q = par->triangles(cT).vertices[(thisCorner+2)%3];
        //printf("p: %d   q: %d\n", p, q);
        if (edgeIntersectsNormalFan(surf->points[from], surf->points[to],
                                    par->vertices(p), par->vertices(q),
                                    //normals[from], normals[to], x)) {
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
            
            //printf("Corner! %d\n", corner);
            if (corner==-1) {
                // parameter polyedge is leaving basegrid triangle 
                // through the opposite edge
                
                // get neighboring triangle
                TriangleIdx neighboringTri = par->getNeighboringTriangle(cT, oppEdge);
                //printf("neighboringTri = %d\n", neighboringTri);
                // if no neighboring triangle --> error
                if (neighboringTri==-1)
                    return false;

                // add intersection nodes on both sides
                
                // better: using getEdge()
                int e = par->triangles(neighboringTri).getCorner(q);

                currType = Node::INTERSECTION_NODE;
                currTri  = neighboringTri;
                lambda   = newLambda;
                enteringEdge = e;

                return true;

            } else {

                // Shouldn't happen: we have to leave the triangle through an
                // edge.  If the parameter edge should go through a domain vertex, there would
                // be a node on that vertex.
                abort();
                
            }
            
        }
        
    }

    printf("no intersection found!\n");
    return false;
}


bool NormalProjector::onSameTriangle(const NodeBundle& a, const NodeBundle& b) const
{
    for (int i=0; i<a.size(); i++)
        for (int j=0; j<b.size(); j++)
            if (a[i].tri==b[j].tri)
                return true;

    return false;
}

bool NormalProjector::onSameTriangle(const TriangleIdx& tri, const NodeBundle& b) const
{
    for (int j=0; j<b.size(); j++)
        if (tri==b[j].tri)
            return true;

    return false;
}

void NormalProjector::insertGhostNodeAtVertex(Parametrization* par, int v, 
                                              int targetTri, const McVec2d& localTargetCoords)
{
    McSmallArray<TriangleIdx, 12> neighbors = par->getTrianglesPerVertex(v);
    //printf("localtargetCoords (%f %f)\n", localTargetCoords.x, localTargetCoords.y);
    for (int i=0; i<neighbors.size(); i++) {

        const DomainTriangle& cT = par->triangles(neighbors[i]);
        
        int corner = cT.getCorner(v);
        McVec2f lTC_Float(localTargetCoords.x, localTargetCoords.y);
        par->addGhostNode(neighbors[i], corner, targetTri, lTC_Float);

    }
 
}

void NormalProjector::addCornerNodeBundle(Parametrization* cS, int v, int nN)
{
    McSmallArray<TriangleIdx, 12> neighbors = cS->getTrianglesPerVertex(v);
    //printf("localtargetCoords (%f %f)\n", localTargetCoords.x, localTargetCoords.y);
    for (int i=0; i<neighbors.size(); i++) {

        const DomainTriangle& cT = cS->triangles(neighbors[i]);
        
        int corner = cT.getCorner(v);
        cS->addCornerNode(neighbors[i], corner, nN);

    }
        
}

bool NormalProjector::computeInverseNormalProjection(const McVec3f& p0_f, const McVec3f& p1_f, const McVec3f& p2_f,
                                                     const McVec3d& n0, const McVec3d& n1, const McVec3d& n2,
                                                     const McVec3f& target, McVec3d& x)
{
    int i;
    const double eps = 1e-6;
    // Fix some initial value
    x.setValue(1.0, 1.0, 1.0);

    // transform to double
    McVec3d p0(p0_f[0], p0_f[1], p0_f[2]);
    McVec3d p1(p1_f[0], p1_f[1], p1_f[2]);
    McVec3d p2(p2_f[0], p2_f[1], p2_f[2]);

    for (i=0; i<10; i++) {

        // compute Newton correction
        McVec3d Fxk = x[0]*(p0-p2) + x[1]*(p1-p2) + x[2]*x[0]*(n0-n2) + x[2]*x[1]*(n1-n2) + x[2]*n2 + p2;// - target;
        Fxk[0] -= target[0];
        Fxk[1] -= target[1];
        Fxk[2] -= target[2];

        //printf("Fxk = (%f %f %f)\n", Fxk[0], Fxk[1], Fxk[2]);

        McMat3d FPrimexk(p0 - p2 + x[2]*(n0-n2),
                         p1 - p2 + x[2]*(n1-n2),
                         x[0]*(n0-n2) + x[1]*(n1-n2) + n2);

        McMat3d FPrimexkInv = FPrimexk.inverse();

        McVec3d newtonCorrection; // = (-1) * FPrimexk.inverse() * Fxk;
        
        FPrimexkInv.multMatrixVec(-Fxk, newtonCorrection);

        x += newtonCorrection;

        //if (debug) {
            //printf("corr = (%f %f %f)\n", newtonCorrection[0], newtonCorrection[1], newtonCorrection[2]);
        //  printf("x = (%f %f %f)\n", x[0], x[1], x[2]);
        //}
    }

    if (x[0]>=-eps && x[1]>=-eps && (x[0]+x[1] <=1+eps)){
        //printf("x = (%f %f %f)\n", x[0], x[1], x[2]);
        return true;
    } 
        
    return false;
}

bool NormalProjector::edgeIntersectsNormalFan(const McVec3f& q0_f, const McVec3f& q1_f,
                                              const McVec3f& p0_f, const McVec3f& p1_f,
                                              const McVec3d& n0, const McVec3d& n1,
                                              McVec3d& x)
{
    int i;
    // transform to double values
    McVec3d q0(q0_f[0], q0_f[1], q0_f[2]);
    McVec3d q1(q1_f[0], q1_f[1], q1_f[2]);
    McVec3d p0(p0_f[0], p0_f[1], p0_f[2]);
    McVec3d p1(p1_f[0], p1_f[1], p1_f[2]);

    // Fix some initial value
    // sometimes it only works when the initial value is an intersection...
    x.setValue(0.5, 0.5, 1);
    McVec3d newtonCorrection;

    //printf("--------- in Newton solver -----------\n");
    for (i=0; i<30; i++) {

        // compute Newton correction

        McVec3d Fxk = p0-q0 + x[0]*(p1-p0) + x[2]*n0 + x[2]*x[0]*(n1-n0) - x[1]*(q1-q0);

        //printf("Fxk: (%f %f %f)\n", Fxk.x, Fxk.y, Fxk.z);

        McMat3d FPrimexk(p1-p0 + x[2]*(n1-n0),
                         q0-q1,
                         n0 + x[0]*(n1-n0));

//         printf("FPrimexk: %f \t %f \t %f\n", FPrimexk[0][0],FPrimexk[0][1],FPrimexk[0][2]);
//         printf("          %f \t %f \t %f\n", FPrimexk[1][0],FPrimexk[1][1],FPrimexk[1][2]);
//         printf("          %f \t %f \t %f\n\n", FPrimexk[2][0],FPrimexk[2][1],FPrimexk[2][2]);

//         printf("det: %f\n", FPrimexk.det());
        McMat3d FPrimexkInv = FPrimexk.inverse();
        
//         printf("FPrimexkInv: %f \t %f \t %f\n", FPrimexkInv[0][0],FPrimexkInv[0][1],FPrimexkInv[0][2]);
//         printf("             %f \t %f \t %f\n", FPrimexkInv[1][0],FPrimexkInv[1][1],FPrimexkInv[1][2]);
//         printf("             %f \t %f \t %f\n\n", FPrimexkInv[2][0],FPrimexkInv[2][1],FPrimexkInv[2][2]);

//         McMat3f Id = FPrimexk*FPrimexkInv;
//         printf("Id: %f \t %f \t %f\n", Id[0][0],Id[0][1],Id[0][2]);
//         printf("    %f \t %f \t %f\n", Id[1][0],Id[1][1],Id[1][2]);
//         printf("    %f \t %f \t %f\n\n", Id[2][0],Id[2][1],Id[2][2]);

        FPrimexkInv.multMatrixVec(-Fxk, newtonCorrection);

        x += newtonCorrection;

//         printf("corr = (%f %f %f)\n", newtonCorrection[0], newtonCorrection[1], newtonCorrection[2]);
//          printf("x = (%f %f %f)\n", x[0], x[1], x[2]);
    }
    //assert(newtonCorrection.x > -50);
    if (x[0]>=0 && x[0]<=1 && x[1]>=0 && x[1]<=1 && newtonCorrection.length()<1e-4){
        //printf("x = (%f %f %f)\n", x[0], x[1], x[2]);
        return true;
    } 
           
    return false;
}

bool NormalProjector::rayIntersectsTriangle(const McVec3d& basePoint, const McVec3d& direction,
                                            const McVec3f& a_, const McVec3f& b_, const McVec3f& c_,
                                            McVec2d& localCoords, double& normalDist, double eps)
{
    const McVec3d &p = basePoint;

    McVec3d a(a_[0], a_[1], a_[2]);
    McVec3d b(b_[0], b_[1], b_[2]);
    McVec3d c(c_[0], c_[1], c_[2]);

    McVec3d e1 = b-a;
    McVec3d e2 = c-a;
    e1.normalize();
    e2.normalize();
    bool parallel = fabs(McMat3d(e1, e2, direction).det()) <eps;
        
        // Cramer's rule
        
        if (!parallel){

            double det = McMat3d(b-a, c-a, direction).det();
            
            // triangle and edge are not parallel
            double nu = McMat3d(b-a, c-a, p-a).det() / det;

            double lambda = McMat3d(p-a, c-a, direction).det() / det;
            if (lambda<-eps) return false;

            double mu = McMat3d(b-a, p-a, direction).det() / det;
            if (mu<-eps) return false;

            if (lambda + mu > 1+eps) 
                return false;
            else {
                localCoords[0] = 1-lambda-mu;
                localCoords[1] = lambda;
                normalDist     = -nu;
//                 printf("---------------------\n");
//                 printf("det: %g\n", det);
//                 printf("lambda: %g,   mu %g\n", lambda, mu);
//                 printf("a: (%g %g %g),  b (%g %g %g)  c (%g %g %g)\n",
//                        a[0], a[1], a[2], b[0], b[1], b[2], c[0], c[1], c[2]);
                //McVec3d where = p + nu*(-direction);
                //printf("where (%f %f %f)\n", where.x, where.y, where.z);
                //McVec3d w6 = a*(1-localCoords[0]-localCoords[1]) + b*localCoords[1] + c*localCoords[0];
                //McVec3d w6 = a*localCoords[0] + b*localCoords[1] + c*(1-localCoords[0]-localCoords[1]);

                //printf("w6: (%f %f %f)\n", w6.x, w6.y, w6.z);
                return true;
            }

        } else {

            // triangle and edge are parallel
            double alpha = McMat3d(b-a, c-a, p-a).det();
            if (alpha<-eps || alpha>eps)
                return false;
            else {
                printf("ray and triangle are parallel!\n");
                return false;

            }
                
        }


}


void NormalProjector::computeVertexNormals()
{
}

NodeIdx NormalProjector::getCornerNode(const DomainTriangle& cT, int corner)
{
    assert(corner>=0 && corner<3);

    for (int i=0; i<cT.nodes.size(); i++)
        if ((cT.nodes[i].isCORNER_NODE() || cT.nodes[i].isGHOST_NODE()) &&
            cT.nodes[i].getCorner()==corner)
            return i;

    return -1;
}

TriangleIdx NormalProjector::getCommonTri(const NodeBundle& a, const NodeBundle& b)
{
    for (int i=0; i<a.size(); i++)
        for (int j=0; j<b.size(); j++)
            if (a[i].tri==b[j].tri)
                return a[i].tri;

    return -1;
}

McSmallArray<TriangleIdx, 2> NormalProjector::getCommonTris(const NodeBundle& a, const NodeBundle& b)
{
    McSmallArray<TriangleIdx, 2> result;

    for (int i=0; i<a.size(); i++)
        for (int j=0; j<b.size(); j++)
            if (a[i].tri==b[j].tri)
                result.append(a[i].tri);

    return result;
}


void NormalProjector::setupEdgePointArrays(Parametrization* par)
{
    int i, j;

    for (i=0; i<par->getNumTriangles(); i++) {

        DomainTriangle& cT = par->triangles(i);

        cT.edgePoints[0].clear();
        cT.edgePoints[1].clear();
        cT.edgePoints[2].clear();
            
        for (j=0; j<cT.nodes.size(); j++) {
                
            Node& cN = cT.nodes[j];
                
            if (cN.isINTERIOR_NODE())
                continue;
                
            if (cN.isCORNER_NODE() || cN.isGHOST_NODE()) {
                int corner = cN.getCorner();
                cT.edgePoints[corner].insert(0, j);
                cT.edgePoints[(corner+2)%3].append(j);
                continue;
            } 
                
            double lambda = cN.getDomainEdgeCoord();
            int domainEdge = cN.getDomainEdge();
            McSmallArray<int, 2>& cEP = cT.edgePoints[domainEdge];
            
            int idx = 0;
            while (idx<cEP.size() && cT.nodes[cEP[idx]].getDomainEdgeCoord(domainEdge)<lambda) {
                idx++;
            }                
            
            cEP.insert(idx, j);
            
        }
        
    }   
    
}

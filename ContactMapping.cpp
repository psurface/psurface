#include <limits>
#include <stdexcept>

#ifdef PSURFACE_STANDALONE
#include "TargetSurface.h"
#else
#include "hxsurface/Surface.h"
#endif

#include <psurface/ContactMapping.h>
#include <psurface/NormalProjector.h>
#include <psurface/StaticMatrix.h>
#include <psurface/DirectionFunction.h>

template <class ctype>
void ContactMapping<2,ctype>::build(const std::vector<std::tr1::array<ctype,2> >& coords1,  ///< The vertex coordinates of the first surface
               const std::vector<std::tr1::array<int,2> >& tri1,       ///< The triangles of the first surface
               const std::vector<std::tr1::array<ctype,2> >& coords2,  ///< The vertices of the second surface
               const std::vector<std::tr1::array<int,2> >& tri2,
                                    const DirectionFunction<2,ctype>* domainDirection,
                                    const DirectionFunction<2,ctype>* targetDirection
               )
{
    int numVertices1 = coords1.size();
    int numVertices2 = coords2.size();
    int nTri1  = tri1.size();
    int nTri2  = tri2.size();

#if 0
    printf("----- 1 -----\n");

    for (int i=0; i<nTri1; i++)
        printf("-- %d %d\n", tri1[2*i], tri1[2*i+1]);

    printf("----- 2 -----\n");

    for (int i=0; i<nTri2; i++)
        printf("-- %d %d\n", tri2[2*i], tri2[2*i+1]);
#endif

    // //////////////////////////////////////////////////
    //   Build domain surface and its normal field
    // //////////////////////////////////////////////////

    psurface_.domainVertices.resize(numVertices1);
    for (int i=0; i<numVertices1; i++) 
        for (int j=0; j<2; j++)
            psurface_.domainVertices[i][j] = coords1[i][j];

    // Build the domain segments
    psurface_.domainSegments.resize(nTri1);
    for (int i=0; i<nTri1; i++) {
        psurface_.domainSegments[i].points[0] = tri1[i][0];
        psurface_.domainSegments[i].points[1] = tri1[i][1];
    }

    // ///////////////////////////////
    //   Build the domain normal field
    // ///////////////////////////////
    std::vector<StaticVector<ctype, 2> > domainNormals;
    std::vector<StaticVector<ctype, 2> > targetNormals;

    computeDiscreteDomainDirections(domainDirection, domainNormals);

    // //////////////////////////////////////////////////
    //   Build range surface and its normal field
    // //////////////////////////////////////////////////

    // first mark the vertices that are actually used
    psurface_.targetVertices.resize(numVertices2);
    for (int i=0; i<numVertices2; i++) 
        for (int j=0; j<2; j++)
            psurface_.targetVertices[i][j] = coords2[i][j];

    // /////////////////////////////////////////////////////
    //   Build the segments-per-vertex arrays
    // /////////////////////////////////////////////////////

    std::vector<std::tr1::array<int, 2> > segPerVertex1(psurface_.domainVertices.size());
    for (size_t i=0; i<segPerVertex1.size(); i++)
        segPerVertex1[i][0] = segPerVertex1[i][1] = -1;

    for (int i=0; i<nTri1; i++) {

        //printf("segment %d:  %d %d  --  %d %d\n", i, tri2[2*i], tri2[2*i+1], used2[tri2[2*i]],used2[tri2[2*i+1]]);
        for (int j=0; j<2; j++) {

            int p = tri1[i][j];
            if (segPerVertex1[p][0]==-1)
                segPerVertex1[p][0] = i;
            else
                segPerVertex1[p][1] = i;

        }

    }

    // use this to construct the neighbor relationships between segments
    for (size_t i=0; i<psurface_.domainSegments.size(); i++) {

        int vertex0 = psurface_.domainSegments[i].points[0];
        int other0 = (segPerVertex1[vertex0][0] == i) ? segPerVertex1[vertex0][1] : segPerVertex1[vertex0][0];
        psurface_.domainSegments[i].neighbor[0] = other0;

        int vertex1 = psurface_.domainSegments[i].points[1];
        int other1 = (segPerVertex1[vertex1][0] == i) ? segPerVertex1[vertex1][1] : segPerVertex1[vertex1][0];
        psurface_.domainSegments[i].neighbor[1] = other1;

        //printf("Segment %d neighbors:  %d  %d\n", i, other0, other1);
    }

    

    // Build the segments-per-vertex arrays for the target vertices
    std::vector<std::tr1::array<int, 2> > segPerVertex2(psurface_.targetVertices.size());
    for (size_t i=0; i<segPerVertex2.size(); i++)
        segPerVertex2[i][0] = segPerVertex2[i][1] = -1;

    for (int i=0; i<nTri2; i++) {

        //printf("segment %d:  %d %d  --  %d %d\n", i, tri2[2*i], tri2[2*i+1], tri2[2*i], tri2[2*i+1]);
        for (int j=0; j<2; j++) {

            int p = tri2[i][j];
            if (segPerVertex2[p][0]==-1)
                segPerVertex2[p][0] = i;
            else
                segPerVertex2[p][1] = i;

        }

    }

    computeDiscreteTargetDirections(tri2, targetDirection, targetNormals);

    // ///////////////////////////////////////////////////////////////////////
    //   Project the vertices of the target surface onto the domain surface
    // ///////////////////////////////////////////////////////////////////////
    const ctype eps = 1e-10;

    for (size_t i=0; i<psurface_.targetVertices.size(); i++) {

        ctype bestLocalPos = std::numeric_limits<ctype>::max();  // init to something
        int bestSegment = -1;
        ctype bestDist = std::numeric_limits<ctype>::max();
        
        for (int j=0; j<psurface_.domainSegments.size(); j++) {

            const StaticVector<ctype,2>& p0 = psurface_.domainVertices[psurface_.domainSegments[j].points[0]];
            const StaticVector<ctype,2>& p1 = psurface_.domainVertices[psurface_.domainSegments[j].points[1]];

            const StaticVector<ctype,2>& n0 = domainNormals[psurface_.domainSegments[j].points[0]];
            const StaticVector<ctype,2>& n1 = domainNormals[psurface_.domainSegments[j].points[1]];

            ctype local; // the unknown...

            if (NormalProjector<ctype>::computeInverseNormalProjection(p0, p1, n0, n1, 
                                                                        psurface_.targetVertices[i], local)) {

                // We want that the line from the domain surface to its projection
                // approaches the target surface from the front side, i.e., it should
                // not pass through the body represented by the target surface.
                // We do a simplified test by comparing the connecting segment
                // with the normal at the target surface and the normal at the
                // domain surface
                /** \todo Rewrite this once we have expression templates */
                StaticVector<ctype,2> base;
                StaticVector<ctype, 2> baseNormal;
                StaticVector<ctype, 2> segment;

                for (int k=0; k<2; k++) {
                    base[k]       = (1-local)*p0[k] + local*p1[k];
                    baseNormal[k] = (1-local)*n0[k] + local*n1[k];
                    segment[k]    = psurface_.targetVertices[i][k] - base[k];
                }

                ctype distance = segment.length2();

                if (segment.dot(targetNormals[i]) > -0.0001
                    && segment.dot(baseNormal) > -0.0001
                    && distance > 1e-8) {
                    //printf("aborting %g %g %g\n", segment * targetNormals[i], segment * baseNormal, distance);
                    continue;
                }

                // There may be several inverse orthogonal projections.
                // We want the shortest one.

                if (distance < bestDist) {

                    bestDist = distance;
                    bestLocalPos = local;
                    bestSegment  = j;

                }
                
            }

        }

        // /////////////////////////////////////////////
        //   We have found a valid projection
        // /////////////////////////////////////////////
        if (bestSegment != -1) {

            typename PSurface<1,ctype>::DomainSegment& bS = psurface_.domainSegments[bestSegment];

            if (bestLocalPos < eps) {

                // Insert as new first element
                bS.nodes.insert(bS.nodes.begin(), typename PSurface<1,ctype>::Node(0, 1, true, true, segPerVertex2[i][0], segPerVertex2[i][1]));

                // Look for left neighbor segment
                if (psurface_.domainSegments[bestSegment].neighbor[0] != -1) {

                    psurface_.domainSegments[psurface_.domainSegments[bestSegment].neighbor[0]].nodes.push_back( typename PSurface<1,ctype>::Node(1, 0, true, true, 
                                                                                                  segPerVertex2[i][0], segPerVertex2[i][1]) );
                    
                }

            } else if (bestLocalPos > 1-eps) {
                    
                typename PSurface<1,ctype>::Node newNode(1, 0, true, true, segPerVertex2[i][0], segPerVertex2[i][1]);
                bS.nodes.push_back(newNode);

                // Look for right neighbor segment
                if (psurface_.domainSegments[bestSegment].neighbor[1] != -1) {

                    typename PSurface<1,ctype>::DomainSegment& rightNeighborSegment = psurface_.domainSegments[psurface_.domainSegments[bestSegment].neighbor[1]];
                    rightNeighborSegment.nodes.insert(rightNeighborSegment.nodes.begin(), 
                                                      typename PSurface<1,ctype>::Node(0, 1, true, true, segPerVertex2[i][0], segPerVertex2[i][1]));

                }
            } else {
                int nNodes = bS.nodes.size();
                
                bS.nodes.resize(nNodes+1);
                int j=nNodes-1;
                for (; j>=0; j--) {
                    if (bS.nodes[j].domainLocalPosition > bestLocalPos)
                        bS.nodes[j+1] = bS.nodes[j];
                    else
                        break;
                }             
                
                bS.nodes[j+1] = typename PSurface<1,ctype>::Node(bestLocalPos, 0, false, true, segPerVertex2[i][0], segPerVertex2[i][1]);
                
            }
            
        }

    }

    // //////////////////////////////////////////////////////////////////////
    //   Insert missing nodes that belong to vertices of the domain segment
    // //////////////////////////////////////////////////////////////////////

    for (int i=0; i<psurface_.domainSegments.size(); i++) {

        typename PSurface<1,ctype>::DomainSegment& cS = psurface_.domainSegments[i];

        // Insert node belonging to domain vertex to the segment to the left of the vertex
        if (cS.nodes.size()==0
            || !cS.nodes[0].isNodeOnVertex
            || (cS.nodes.size()==1 && cS.nodes[0].isNodeOnVertex && cS.nodes[0].domainLocalPosition > 1-eps)) {

            ctype rangeLocalPosition;
            int rangeSegment;

            if (NormalProjector<ctype>::normalProjection(psurface_.domainVertices[cS.points[0]], domainNormals[cS.points[0]],
                                                          rangeSegment, rangeLocalPosition,
                                                          tri2, coords2)) {

                typename PSurface<1,ctype>::Node newNode(0, rangeLocalPosition, true, false, rangeSegment, rangeSegment);
                
                cS.nodes.insert(cS.nodes.begin(), newNode);

            }

        }

        // Insert node belonging to domain vertex to the segment to the right of the vertex
        if (cS.nodes.size()==0
            || !cS.nodes.back().isNodeOnVertex
            || (cS.nodes.size()==1 && cS.nodes[0].isNodeOnVertex && cS.nodes[0].domainLocalPosition < eps)) {

            ctype rangeLocalPosition;
            int rangeSegment;

            if (NormalProjector<ctype>::normalProjection(psurface_.domainVertices[cS.points[1]], domainNormals[cS.points[1]],
                                                          rangeSegment, rangeLocalPosition,
                                                          tri2, coords2)) {

                typename PSurface<1,ctype>::Node newNode(1, rangeLocalPosition, true, false, rangeSegment, rangeSegment);
       
                cS.nodes.push_back(newNode);
            }

        }

    }

#if 0
    for (int i=0; i<psurface_.domainSegments.size(); i++) {
        printf(" --- segment %d ---   (%d  -->  %d)\n", i, 
               psurface_.domainSegments[i].points[0],psurface_.domainSegments[i].points[1]);
        for (int j=0; j<psurface_.domainSegments[i].nodes.size(); j++)
            std::cout << psurface_.domainSegments[i].nodes[j];

        std::cout << std::endl;
    }
#endif    

    // /////////////////////////////////////////////////////
    //   Insert edges
    // /////////////////////////////////////////////////////

    /** \todo Only works if the relevant domain is a single connected component */
    for (int i=0; i<psurface_.domainSegments.size(); i++) {

        std::vector<typename PSurface<1,ctype>::Node>& nodes = psurface_.domainSegments[i].nodes;

        ////////////////////////////////
        for (int j=0; j<int(nodes.size())-1; j++) {

            if (nodes[j].rangeSegments[0] == nodes[j+1].rangeSegments[0])
                nodes[j].rightRangeSegment = nodes[j].rangeSegments[0];
            else if (nodes[j].rangeSegments[0] == nodes[j+1].rangeSegments[1])
                nodes[j].rightRangeSegment = nodes[j].rangeSegments[0];
            else if (nodes[j].rangeSegments[1] == nodes[j+1].rangeSegments[0])
                nodes[j].rightRangeSegment = nodes[j].rangeSegments[1];
            else if (nodes[j].rangeSegments[1] == nodes[j+1].rangeSegments[1])
                nodes[j].rightRangeSegment = nodes[j].rangeSegments[1];
            else
                throw(std::runtime_error("Segment of the PSurface<1> data structure is inconsistent!"));
            
            if (nodes[j].rightRangeSegment == -1)
                throw(std::runtime_error("Segment of the PSurface<1> data structure is inconsistent!"));

        }

    }

}

template <class ctype>
void ContactMapping<2,ctype>::computeDiscreteDomainDirections(const DirectionFunction<2,ctype>* direction,
                                                         std::vector<StaticVector<ctype,2> >& normals)
{
    size_t nElements = psurface_.domainSegments.size();

    if (!direction) {

        // //////////////////////////////////////////////////////////
        //  The contact directions are given as the vertex normals
        // //////////////////////////////////////////////////////////
        normals.resize(psurface_.domainVertices.size());
        for (size_t i=0; i<psurface_.domainVertices.size(); i++) 
            normals[i] = StaticVector<ctype,2>(0);

        for (size_t i=0; i<nElements; i++) {

            // Compute segment normal
            int v0 = psurface_.domainSegments[i].points[0];
            int v1 = psurface_.domainSegments[i].points[1];

            StaticVector<ctype,2> segment;
            segment[0] = psurface_.domainVertices[v1][0] - psurface_.domainVertices[v0][0];
            segment[1] = psurface_.domainVertices[v1][1] - psurface_.domainVertices[v0][1];

            StaticVector<ctype,2> segmentNormal;
            segmentNormal[0] =  segment[1];
            segmentNormal[1] = -segment[0];

            segmentNormal /= segmentNormal.length();

            normals[psurface_.domainSegments[i].points[0]] += segmentNormal;
            normals[psurface_.domainSegments[i].points[1]] += segmentNormal;

        }

        for (size_t i=0; i<normals.size(); i++)
            normals[i] /= normals[i].length();
        
    } else {

        // Sample the provided analytical contact direction field
        normals.resize(psurface_.domainVertices.size());
        for (size_t i=0; i<psurface_.domainVertices.size(); i++) {

            if (dynamic_cast<const AnalyticDirectionFunction<2,ctype>*>(direction))
                normals[i] = (*dynamic_cast<const AnalyticDirectionFunction<2,ctype>*>(direction))(psurface_.domainVertices[i]);
            else if (dynamic_cast<const DiscreteDirectionFunction<2,ctype>*>(direction))
                normals[i] = (*dynamic_cast<const DiscreteDirectionFunction<2,ctype>*>(direction))(i);
            else
                throw(std::runtime_error("Domain direction function not properly set!"));
                
        }

    }

}

template <class ctype>
void ContactMapping<2,ctype>::computeDiscreteTargetDirections(const std::vector<std::tr1::array<int,2> >& elements,
                                                              const DirectionFunction<2,ctype>* direction,
                                                              std::vector<StaticVector<ctype,2> >& normals)
{
    // Build the normal field
    normals.resize(psurface_.targetVertices.size());
    for (size_t i=0; i<psurface_.targetVertices.size(); i++) 
        normals[i] = StaticVector<ctype,2>(0);
    
    if (!direction) {
        
        for (int i=0; i<elements.size(); i++) {
            
            // Compute segment normal
            int v0 = elements[i][0];
            int v1 = elements[i][1];
            StaticVector<ctype,2> segment;
            segment[0] = psurface_.targetVertices[v1][0] - psurface_.targetVertices[v0][0];
            segment[1] = psurface_.targetVertices[v1][1] - psurface_.targetVertices[v0][1];
            
            StaticVector<ctype,2> segmentNormal;
            segmentNormal[0] =  segment[1];
            segmentNormal[1] = -segment[0];
            
            segmentNormal /= segmentNormal.length();
            
            normals[elements[i][0]] += segmentNormal;
            normals[elements[i][1]] += segmentNormal;
            
        }
        
        for (size_t i=0; i<normals.size(); i++)
            normals[i] /= normals[i].length();
        
    } else {
        
        // Sample the provided analytical contact direction field
        normals.resize(psurface_.targetVertices.size());
        for (size_t i=0; i<psurface_.targetVertices.size(); i++) {
            
            if (dynamic_cast<const AnalyticDirectionFunction<2,ctype>*>(direction))
                normals[i] = (*dynamic_cast<const AnalyticDirectionFunction<2,ctype>*>(direction))(psurface_.targetVertices[i]);
            else if (dynamic_cast<const DiscreteDirectionFunction<2,ctype>*>(direction))
                normals[i] = (*dynamic_cast<const DiscreteDirectionFunction<2,ctype>*>(direction))(i);
            else
                throw(std::runtime_error("Target direction function not properly set!"));
            
        }
        
    }

}

template <class ctype>
void ContactMapping<3,ctype>::build(const std::vector<std::tr1::array<ctype,3> >& coords1,  ///< The vertices of the first surface
                         const std::vector<std::tr1::array<int,3> >& tri1,       ///< The triangles of the first surface
                         const std::vector<std::tr1::array<ctype,3> >& coords2,  ///< The vertices of the second surface
                         const std::vector<std::tr1::array<int,3> >& tri2,
                                    const DirectionFunction<3,ctype>* domainDirection,
                                    const DirectionFunction<3,ctype>* targetDirection)
{
    int nVert1 = coords1.size();
    int nVert2 = coords2.size();
    int nTri1  = tri1.size();
    int nTri2  = tri2.size();

    // Create a first Surface object
    Surface* surface1_ = new Surface;

#ifndef PSURFACE_STANDALONE
    // Amira Surface class needs a 'patch' structure
    surface1_->patches.resize(1);
    surface1_->patches[0] = new Surface::Patch;
    surface1_->patches[0]->innerRegion = 0;
    surface1_->patches[0]->outerRegion = 1;
    surface1_->patches[0]->boundaryId  = 0;
#endif

    surface1_->points.resize(nVert1);
    for (int i=0; i<nVert1; i++) 
        for (int j=0; j<3; j++)
            surface1_->points[i][j] = coords1[i][j];

    surface1_->triangles.resize(nTri1);
#ifndef PSURFACE_STANDALONE
    surface1_->patches[0]->triangles.resize(nTri1);
#endif
    for (int i=0; i<nTri1; i++) {

        surface1_->triangles[i].points[0] = tri1[i][0];
        surface1_->triangles[i].points[1] = tri1[i][1];
        surface1_->triangles[i].points[2] = tri1[i][2];

#ifndef PSURFACE_STANDALONE
        surface1_->triangles[i].patch = 0;
        surface1_->patches[0]->triangles[i] = i;
#endif
    }

    // Create a second Surface object
    Surface* surface2_ = new Surface;
    
#ifndef PSURFACE_STANDALONE
    surface2_->patches.resize(1);
    surface2_->patches[0] = new Surface::Patch;
    surface2_->patches[0]->innerRegion = 0;
    surface2_->patches[0]->outerRegion = 1;
    surface2_->patches[0]->boundaryId  = 0;
#endif

    surface2_->points.resize(nVert2);
    for (int i=0; i<nVert2; i++)
        for (int j=0; j<3; j++)
            surface2_->points[i][j] = coords2[i][j];

    surface2_->triangles.resize(nTri2);
#ifndef PSURFACE_STANDALONE
    surface2_->patches[0]->triangles.resize(nTri2);
#endif
    for (int i=0; i<nTri2; i++) {

        surface2_->triangles[i].points[0] = tri2[i][0];
        surface2_->triangles[i].points[1] = tri2[i][1];
        surface2_->triangles[i].points[2] = tri2[i][2];

#ifndef PSURFACE_STANDALONE
        surface2_->triangles[i].patch = 0;
        surface2_->patches[0]->triangles[i] = i;
#endif
    }

//     ContactToolBox<ctype>::buildContactSurface(&psurface_, surface1_, surface2_,
//                                                domainDirection, targetDirection);
   // set up parametrization
    psurface_.surface = const_cast<Surface*>(surface2_);
    psurface_.patches.resize(1);
    psurface_.patches[0].innerRegion = 0;
    psurface_.patches[0].outerRegion = 1;
    psurface_.patches[0].boundaryId  = 0;
            
    // ///////
    const_cast<Surface*>(surface1_)->removeUnusedPoints();
    const_cast<Surface*>(surface2_)->removeUnusedPoints();
    
    std::cout << surface1_->points.size() << " resp. "
              << surface2_->points.size() << " contact nodes found!" << std::endl;

    std::cout << "Contact patches contain " << surface1_->triangles.size() 
              << " (resp. " << surface2_->triangles.size() << ") triangles." << std::endl;
    
    // the nonmortar side becomes the base grid of the parametrization
    for (size_t i=0; i<surface1_->points.size(); i++) {
        StaticVector<ctype,3> newVertex;
        for (int j=0; j<3; j++)
            newVertex[j] = surface1_->points[i][j];
        psurface_.newVertex(newVertex);
    }
    
    for (size_t i=0; i<surface1_->triangles.size(); i++) {
        
        int newTri = psurface_.createSpaceForTriangle(surface1_->triangles[i].points[0],
                                                  surface1_->triangles[i].points[1],
                                                  surface1_->triangles[i].points[2]);
        psurface_.integrateTriangle(newTri);
        psurface_.triangles(newTri).patch = 0;
        
    }
    
    // compute projection
    NormalProjector<ctype> projector(&psurface_);
    
    projector.project(surface2_, domainDirection, targetDirection);

}


// ///////////////////////////////////////////////////////////////////////
//   Explicitly instantiate 'float' and 'double' versions of this code
// ///////////////////////////////////////////////////////////////////////

template class ContactMapping<2,float>;
template class ContactMapping<2,double>;

template class ContactMapping<3,float>;
template class ContactMapping<3,double>;


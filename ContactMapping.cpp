#include <limits>
#include <stdexcept>

#include <psurface/ContactMapping.h>
#include <psurface/StaticMatrix.h>

template <class ctype>
void ContactMapping<2,ctype>::build(const std::vector<std::tr1::array<double,2> >& coords1,  ///< The vertex coordinates of the first surface
               const std::vector<std::tr1::array<int,2> >& tri1,       ///< The triangles of the first surface
               const std::vector<std::tr1::array<double,2> >& coords2,  ///< The vertices of the second surface
               const std::vector<std::tr1::array<int,2> >& tri2,
               float epsilon,
               void (*obsDirections)(const double* pos, double* dir)
               )
{
    int nVert1 = coords1.size();
    int nVert2 = coords2.size();
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

    // first mark the vertices that are actually used
    std::vector<int> used1(nVert1);
    for (size_t i=0; i<used1.size(); i++)
        used1[i] = -1;

    for (int i=0; i<nTri1; i++) {
        used1[tri1[i][0]] = 1;
        used1[tri1[i][1]] = 1;
    }

    int numVertices1 = 0;
    for (int i=0; i<nVert1; i++)
        if (used1[i]==1)
            used1[i] = numVertices1++;

    vertices.resize(numVertices1);
    for (int i=0; i<nVert1; i++) 
        if (used1[i]!=-1)
            for (int j=0; j<2; j++)
                vertices[used1[i]][j] = coords1[i][j];

    // Build the domain segments
    domainSegments.resize(nTri1);
    for (int i=0; i<nTri1; i++) {
        domainSegments[i].points[0] = used1[tri1[i][0]];
        domainSegments[i].points[1] = used1[tri1[i][1]];
    }

    // ///////////////////////////////
    //   Build the normal field
    // ///////////////////////////////
    if (!obsDirections) {

        // //////////////////////////////////////////////////////////
        //  The contact directions are given as the vertex normals
        // //////////////////////////////////////////////////////////
        domainNormals.resize(vertices.size());
        for (size_t i=0; i<vertices.size(); i++) 
            domainNormals[i] = 0;

        for (int i=0; i<nTri1; i++) {

            // Compute segment normal
            int v0 = tri1[i][0];
            int v1 = tri1[i][1];

            StaticVector<double,2> segment;
            segment[0] = coords1[v1][0] - coords1[v0][0];
            segment[1] = coords1[v1][1] - coords1[v0][1];

            StaticVector<double,2> segmentNormal;
            segmentNormal[0] =  segment[1];
            segmentNormal[1] = -segment[0];

            segmentNormal /= segmentNormal.length();

            domainNormals[used1[tri1[i][0]]]   += segmentNormal;
            domainNormals[used1[tri1[i][1]]] += segmentNormal;

//             std::cout << "Normal: " << segmentNormal << "   ";
//             printf("added to %d %d   ---   %d %d\n", tri1[2*i], tri1[2*i+1],
//                    used1[tri1[2*i]], used1[tri1[2*i+1]]);
        }

        for (size_t i=0; i<domainNormals.size(); i++) {
            domainNormals[i] /= domainNormals[i].length();
            //std::cout << "Normal " << i << ":   " << domainNormals[i] << std::endl;
        }
        

    } else {

        // Sample the provided analytical contact direction field
        domainNormals.resize(vertices.size());
        for (size_t i=0; i<vertices.size(); i++) 
            obsDirections(&vertices[i][0], &domainNormals[i][0]);

    }


    // //////////////////////////////////////////////////
    //   Build range surface and its normal field
    // //////////////////////////////////////////////////

    // first mark the vertices that are actually used
    std::vector<int> used2(nVert2);
    for (size_t i=0; i<used2.size(); i++)
        used2[i] = -1;

    for (int i=0; i<nTri2; i++) {
        used2[tri2[i][0]] = 1;
        used2[tri2[i][1]] = 1;
    }

    int numVertices2 = 0;
    for (int i=0; i<nVert2; i++)
        if (used2[i]==1)
            used2[i] = numVertices2++;

    targetVertices.resize(numVertices2);
    for (int i=0; i<nVert2; i++) 
        if (used2[i]!=-1)
            for (int j=0; j<2; j++)
                targetVertices[used2[i]][j] = coords2[i][j];

    // /////////////////////////////////////////////////////
    //   Build the segments-per-vertex arrays
    // /////////////////////////////////////////////////////

    std::vector<std::tr1::array<int, 2> > segPerVertex1(vertices.size());
    for (size_t i=0; i<segPerVertex1.size(); i++)
        segPerVertex1[i][0] = segPerVertex1[i][1] = -1;

    for (int i=0; i<nTri1; i++) {

        //printf("segment %d:  %d %d  --  %d %d\n", i, tri2[2*i], tri2[2*i+1], used2[tri2[2*i]],used2[tri2[2*i+1]]);
        for (int j=0; j<2; j++) {

            int p = used1[tri1[i][j]];
            if (segPerVertex1[p][0]==-1)
                segPerVertex1[p][0] = i;
            else
                segPerVertex1[p][1] = i;

        }

    }

    // use this to construct the neighbor relationships between segments
    for (size_t i=0; i<domainSegments.size(); i++) {

        int vertex0 = domainSegments[i].points[0];
        int other0 = (segPerVertex1[vertex0][0] == i) ? segPerVertex1[vertex0][1] : segPerVertex1[vertex0][0];
        domainSegments[i].neighbor[0] = other0;

        int vertex1 = domainSegments[i].points[1];
        int other1 = (segPerVertex1[vertex1][0] == i) ? segPerVertex1[vertex1][1] : segPerVertex1[vertex1][0];
        domainSegments[i].neighbor[1] = other1;

        //printf("Segment %d neighbors:  %d  %d\n", i, other0, other1);
    }

    

    // Build the segments-per-vertex arrays for the target vertices
    std::vector<std::tr1::array<int, 2> > segPerVertex2(targetVertices.size());
    for (size_t i=0; i<segPerVertex2.size(); i++)
        segPerVertex2[i][0] = segPerVertex2[i][1] = -1;

    for (int i=0; i<nTri2; i++) {

        //printf("segment %d:  %d %d  --  %d %d\n", i, tri2[2*i], tri2[2*i+1], used2[tri2[2*i]],used2[tri2[2*i+1]]);
        for (int j=0; j<2; j++) {

            int p = used2[tri2[i][j]];
            if (segPerVertex2[p][0]==-1)
                segPerVertex2[p][0] = i;
            else
                segPerVertex2[p][1] = i;

        }

    }

#if 0
    for (int i=0; i<segPerVertex2.size(); i++)
        printf("i %d:  %d  %d\n", i, segPerVertex2[i][0], segPerVertex2[i][1]);

    exit(0);
#endif

    // Build the normal field
    targetNormals.resize(targetVertices.size());
    for (size_t i=0; i<targetVertices.size(); i++) 
        targetNormals[i] = 0;

    for (int i=0; i<nTri2; i++) {

        // Compute segment normal
        int v0 = tri2[i][0];
        int v1 = tri2[i][1];
        StaticVector<double,2> segment;
        segment[0] = coords2[v1][0] - coords2[v0][0];
        segment[1] = coords2[v1][1] - coords2[v0][1];

        StaticVector<double,2> segmentNormal;
        segmentNormal[0] =  segment[1];
        segmentNormal[1] = -segment[0];

        segmentNormal /= segmentNormal.length();

        targetNormals[used2[tri2[i][0]]] += segmentNormal;
        targetNormals[used2[tri2[i][1]]] += segmentNormal;

    }

    for (size_t i=0; i<targetNormals.size(); i++) {
        targetNormals[i] /= targetNormals[i].length();
        //std::cout << "Normal " << i << ":   " << targetNormals[i] << std::endl;
    }
    //exit(0);

    // ///////////////////////////////////////////////////////////////////////
    //   Project the vertices of the target surface onto the domain surface
    // ///////////////////////////////////////////////////////////////////////
    const double eps = 1e-10;

    for (size_t i=0; i<targetVertices.size(); i++) {

        double bestLocalPos = std::numeric_limits<double>::max();  // init to something
        int bestSegment = -1;
        double bestDist = std::numeric_limits<double>::max();
        
        for (int j=0; j<domainSegments.size(); j++) {

            const StaticVector<double,2>& p0 = vertices[domainSegments[j].points[0]];
            const StaticVector<double,2>& p1 = vertices[domainSegments[j].points[1]];

            const StaticVector<double,2>& n0 = domainNormals[domainSegments[j].points[0]];
            const StaticVector<double,2>& n1 = domainNormals[domainSegments[j].points[1]];

            double local; // the unknown...

            if (computeInverseNormalProjection(p0, p1, n0, n1, 
                                               targetVertices[i], local)) {

                // We want that the line from the domain surface to its projection
                // approaches the target surface from the front side, i.e., it should
                // not pass through the body represented by the target surface.
                // We do a simplified test by comparing the connecting segment
                // with the normal at the target surface and the normal at the
                // domain surface
                /** \todo Rewrite this once we have expression templates */
                StaticVector<double,2> base;
                StaticVector<double, 2> baseNormal;
                StaticVector<double, 2> segment;

                for (int k=0; k<2; k++) {
                    base[k]       = (1-local)*p0[k] + local*p1[k];
                    baseNormal[k] = (1-local)*n0[k] + local*n1[k];
                    segment[k]    = targetVertices[i][k] - base[k];
                }

                double distance = segment.length2();

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

            DomainSegment& bS = domainSegments[bestSegment];

            if (bestLocalPos < eps) {

                // Insert as new first element
                bS.nodes.insert(bS.nodes.begin(), Node(0, 1, true, true, segPerVertex2[i][0], segPerVertex2[i][1]));

                // Look for left neighbor segment
                if (domainSegments[bestSegment].neighbor[0] != -1) {

                    domainSegments[domainSegments[bestSegment].neighbor[0]].nodes.push_back( Node(1, 0, true, true, 
                                                                                                  segPerVertex2[i][0], segPerVertex2[i][1]) );
                    
                }

            } else if (bestLocalPos > 1-eps) {
                    
                Node newNode(1, 0, true, true, segPerVertex2[i][0], segPerVertex2[i][1]);
                bS.nodes.push_back(newNode);

                // Look for right neighbor segment
                if (domainSegments[bestSegment].neighbor[1] != -1) {

                    DomainSegment& rightNeighborSegment = domainSegments[domainSegments[bestSegment].neighbor[1]];
                    rightNeighborSegment.nodes.insert(rightNeighborSegment.nodes.begin(), 
                                                      Node(0, 1, true, true, segPerVertex2[i][0], segPerVertex2[i][1]));

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
                
                bS.nodes[j+1] = Node(bestLocalPos, 0, false, true, segPerVertex2[i][0], segPerVertex2[i][1]);
                
            }
            
        }

    }

    // //////////////////////////////////////////////////////////////////////
    //   Insert missing nodes that belong to vertices of the domain segment
    // //////////////////////////////////////////////////////////////////////

    for (int i=0; i<domainSegments.size(); i++) {

        DomainSegment& cS = domainSegments[i];

        // Insert node belonging to domain vertex to the segment to the left of the vertex
        if (cS.nodes.size()==0
            || !cS.nodes[0].isNodeOnVertex
            || (cS.nodes.size()==1 && cS.nodes[0].isNodeOnVertex && cS.nodes[0].domainLocalPosition > 1-eps)) {

            double rangeLocalPosition;
            int rangeSegment;

            if (normalProjection(vertices[cS.points[0]], domainNormals[cS.points[0]],
                                 rangeSegment, rangeLocalPosition,
                                 tri2, coords2)) {

                Node newNode(0, rangeLocalPosition, true, false, rangeSegment, rangeSegment);
                
                cS.nodes.insert(cS.nodes.begin(), newNode);

            }

        }

        // Insert node belonging to domain vertex to the segment to the right of the vertex
        if (cS.nodes.size()==0
            || !cS.nodes.back().isNodeOnVertex
            || (cS.nodes.size()==1 && cS.nodes[0].isNodeOnVertex && cS.nodes[0].domainLocalPosition < eps)) {

            double rangeLocalPosition;
            int rangeSegment;

            if (normalProjection(vertices[cS.points[1]], domainNormals[cS.points[1]],
                                 rangeSegment, rangeLocalPosition,
                                 tri2, coords2)) {

                Node newNode(1, rangeLocalPosition, true, false, rangeSegment, rangeSegment);
       
                cS.nodes.push_back(newNode);
            }

        }

    }

#if 0
    for (int i=0; i<domainSegments.size(); i++) {
        printf(" --- segment %d ---   (%d  -->  %d)\n", i, 
               domainSegments[i].points[0],domainSegments[i].points[1]);
        for (int j=0; j<domainSegments[i].nodes.size(); j++)
            std::cout << domainSegments[i].nodes[j];

        std::cout << std::endl;
    }
#endif    

    // /////////////////////////////////////////////////////
    //   Insert edges
    // /////////////////////////////////////////////////////

    /** \todo Only works if the relevant domain is a single connected component */
    for (int i=0; i<domainSegments.size(); i++) {

        std::vector<Node>& nodes = domainSegments[i].nodes;

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
                throw(std::runtime_error("Segmentx of the Contact Mapping Data structure is inconsistent!"));
            
            if (nodes[j].rightRangeSegment == -1)
                throw(std::runtime_error("Segmentx of the Contact Mapping Data structure is inconsistent!"));

        }

    }

}

template <class ctype>
void ContactMapping<2,ctype>::getOverlaps(std::vector<IntersectionPrimitive<1,float> >& overlaps)
{
    for (int i=0; i<domainSegments.size(); i++) {

        const DomainSegment&    cS = domainSegments[i];
        const std::vector<Node>& nodes = domainSegments[i].nodes;

        if (!isCompletelyCovered(i))
            continue;

        ////////////////////////////////
        for (int j=0; j<int(nodes.size())-1; j++) {
            
            /** \todo Should be in here for true edge handling */
            // Don't do anything if the current pair of points is not connected by an edge
            if (nodes[j].rightRangeSegment == -1)
                continue;

            // //////////////////////////////////////////////
            // Assemble new overlap
            // //////////////////////////////////////////////
            IntersectionPrimitive<1,float> newOverlap;
            newOverlap.tris[0] = i;
            newOverlap.tris[1] = nodes[j].rightRangeSegment;
            
            newOverlap.localCoords[0][0][0] = nodes[j].domainLocalPosition;
            newOverlap.localCoords[0][1][0] = nodes[j+1].domainLocalPosition;
            
            // if the target of a node is a vertex on the target surface, its
            // rangeLocalPosition is always 0.  But its equivalent coordinate
            // in the two overlaps that contain it has to be once 1 and once zero.
            // That explains the following conditional clause
            newOverlap.localCoords[1][0][0] = (nodes[j].isNodeOnTargetVertex) ? 1 : nodes[j].rangeLocalPosition;

            newOverlap.localCoords[1][1][0] = nodes[j+1].rangeLocalPosition;
            
            // Compute the world position of the overlap on the domain side
            /** \todo Rewrite this once we have expression templates */
            newOverlap.points[0][0] = vertices[cS.points[0]][0] * (1-nodes[j].domainLocalPosition)
                + vertices[cS.points[1]][0] * nodes[j].domainLocalPosition;
            newOverlap.points[0][1] = vertices[cS.points[0]][1] * (1-nodes[j].domainLocalPosition)
                + vertices[cS.points[1]][1] * nodes[j].domainLocalPosition;

            newOverlap.points[1][0] = vertices[cS.points[0]][0] * (1-cS.nodes[j+1].domainLocalPosition)
                + vertices[cS.points[1]][0] * cS.nodes[j+1].domainLocalPosition;
            newOverlap.points[1][1] = vertices[cS.points[0]][1] * (1-cS.nodes[j+1].domainLocalPosition)
                + vertices[cS.points[1]][1] * cS.nodes[j+1].domainLocalPosition;
            
            overlaps.push_back(newOverlap);
        }
        
    }

#if 0
    for (int i=0; i<overlaps.size(); i++)
        printf("overlap %d,   nonmortar (%g  -->  %g),    mortar (%g  -->  %g)\n", i,
               overlaps[i].localCoords[0][0][0], overlaps[i].localCoords[0][1][0],
               overlaps[i].localCoords[1][0][0], overlaps[i].localCoords[1][1][0]);
#endif
//     exit(0);
}

template <class ctype>
bool ContactMapping<2,ctype>::isCompletelyCovered(int i) const
{
    // A domain segment is completely covered by the mapping if the first
    // and last node coincide with the ends of the domain boundary segment and if
    // every pair of adjacent nodes on the boundary segment is connected by an
    // edge.

    const DomainSegment& cS = domainSegments[i];

    if (cS.nodes.size()<2)
        return false;

    if (!cS.nodes[0].isNodeOnVertex || !cS.nodes[cS.nodes.size()-1].isNodeOnVertex)
        return false;

#if 0
    /** \todo Should be in here for true edge handling */
    for (int i=0; i<cS.nodes.size()-1; i++)
        if (cS.nodes[i].rightRangeSegment == -1)
            return false;
#endif    

    return true;
}

template <class ctype>
bool ContactMapping<2,ctype>::computeInverseNormalProjection(const StaticVector<double,2>& p0,
                                                       const StaticVector<double,2>& p1,
                                                       const StaticVector<double,2>& n0,
                                                       const StaticVector<double,2>& n1,
                                                       const StaticVector<double,2>& q,
                                                       double& local)
{
    double a = (p1[1]-p0[1])*(n1[0]-n0[0]) - (p1[0]-p0[0])*(n1[1]-n0[1]);
    double b = -(q[1]-p0[1])*(n1[0]-n0[0]) + (p1[1]-p0[1])*n0[0] + (q[0]-p0[0])*(n1[1]-n0[1]) - (p1[0]-p0[0])*n0[1];
    double c = -(q[1]-p0[1])*n0[0] + (q[0]-p0[0])*n0[1];

    // Is the quadratic formula degenerated to a linear one?
    if (std::abs(a) < 1e-10) {
        local = -c/b;
        //printf("mu:  %g,  old local %g\n", mu, ((q[0]-p0[0]) / (p1[0]-p0[0])));
        
        return local >= 0 && local <= 1;
    }

    // The abc-formula
    double mu_0 = (-b + std::sqrt(b*b - 4*a*c))/(2*a);
    double mu_1 = (-b - std::sqrt(b*b - 4*a*c))/(2*a);

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
bool ContactMapping<2,ctype>::normalProjection(const StaticVector<double,2>& base,
                                         const StaticVector<double,2>& direction,
                                         int& bestSegment,
                                         double& rangeLocalPosition,
                                         const std::vector<std::tr1::array<int,2> >& targetSegments,
                                         const std::vector<std::tr1::array<double, 2> >& coords) const
{
    bestSegment = -1;
    int nTargetSegments = targetSegments.size();
    double bestDistance = std::numeric_limits<double>::max();

    for (int i=0; i<nTargetSegments; i++) {

        StaticVector<double,2> p0, p1;
        p0[0] = coords[targetSegments[i][0]][0];
        p0[1] = coords[targetSegments[i][0]][1];

        p1[0] = coords[targetSegments[i][1]][0];
        p1[1] = coords[targetSegments[i][1]][1];

        double distance, targetLocal;
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
bool ContactMapping<2,ctype>::
rayIntersectsLine(const StaticVector<double, 2>& basePoint, 
                  const StaticVector<double, 2>& direction,
                  const StaticVector<double, 2>& a, 
                  const StaticVector<double, 2>& b, 
                  double& distance, double& targetLocal) const
{
    // we solve the equation basePoint + x_0 * normal = a + x_1 * (b-a)

    StaticMatrix<double,2> mat;
    mat[0][0] = direction[0];
    mat[1][0] = direction[1];
    mat[0][1] = a[0]-b[0];
    mat[1][1] = a[1]-b[1];

    /** \todo Easier with expression templates */
    StaticVector<double,2> rhs;
    rhs[0] = a[0]-basePoint[0];
    rhs[1] = a[1]-basePoint[1];

    StaticVector<double,2> x;

    // Solve the system.  If it is singular the normal and the segment
    // are parallel and there is no intersection

    double detinv = mat[0][0]*mat[1][1]-mat[0][1]*mat[1][0];
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

// ///////////////////////////////////////////////////////////////////////
//   Explicitly instantiate 'float' and 'double' versions of this code
// ///////////////////////////////////////////////////////////////////////

template class ContactMapping<2,float>;
template class ContactMapping<2,double>;

#include <psurface/PSurface.h>
#include <psurface/ContactToolBox.h>
#include <vector>


void ContactToolBox::extractMergedGrid(PSurface<2,float>* psurface,
                                       std::vector<IntersectionPrimitive<2,float> >& mergedGrid)
{
    if (psurface->getNumTriangles()==0)
        return;


    // ///////////////////////////////////////////////////
    // Set up point location structure
    // Can we use the routine in PSurface<2,float> ???
    // ///////////////////////////////////////////////////
    for (int i=0; i<psurface->getNumTriangles(); i++) {

        DomainTriangle<float>& cT = psurface->triangles(i);

        cT.insertExtraEdges();

        for (int k=0; k<cT.nodes.size(); k++)
            cT.makeCyclicGeometrically(cT.nodes[k]);

        // the standard insertExtraEdges can miss edges if ghost nodes are present.  
        // We insert them now
        for (int k=0; k<3; k++){
            // size() returns an unsigned type, which underflows if edgePoints[k] is empty
            for (int l=0; l<((int)cT.edgePoints[k].size())-1; l++) {

                PlaneParam<float>::DirectedEdgeIterator cE = cT.getDirectedEdgeIterator(cT.edgePoints[k][l], cT.edgePoints[k][l+1]);

                if (cE.isValid() && cE.getDPrev().from() != cE.getONext().to())
                    cT.addEdge(cE.getONext().to(), cE.to(), true);
                
            }
            
        }

        // and since we have added more edges, we have to redo the cyclic ordering
        for (int k=0; k<cT.nodes.size(); k++)
            cT.makeCyclicGeometrically(cT.nodes[k]);
        
        //
        for (int k=0; k<3; k++){
            for (int l=0; l<cT.edgePoints[k].size(); l++) {
                
                if (!cT.nodes[cT.edgePoints[k][l]].isOnCorner()) {
                    cT.nodes[cT.edgePoints[k][l]].setDomainEdge(k);
                    cT.nodes[cT.edgePoints[k][l]].setDomainEdgePosition(l);
                }
            }
            
        }
        
    }

    psurface->surface->computeTrianglesPerPoint();

    // get the array that relates the base grid triangles with the whole nonmortar surface
    std::vector<int> nonMortarTargetTris = psurface->domainSurfaceTriangleNumbers;

    //
    for (int i=0; i<psurface->getNumTriangles(); i++) {

        const DomainTriangle<float>& cT = psurface->triangles(i);

        if (cT.nodes.size()<3)
            continue;

        ////////////////////////////////
        PlaneParam<float>::TriangleIterator cPT;
        for (cPT = cT.firstTriangle(); cPT.isValid(); ++cPT) {

            
            int targetTri = -1;
            try {
                targetTri = psurface->getImageSurfaceTriangle(i, cPT.vertices());
            } catch (PSurface<2,float>::ParamError){
                printf("exception caught!\n");
                targetTri = -1;
            }

            if (targetTri==-1)
                continue;
            
            assert(targetTri>=0 && targetTri<psurface->surface->triangles.size());


            // //////////////////////////////////////////////
            // Assemble the triangles
            // //////////////////////////////////////////////
            mergedGrid.push_back(IntersectionPrimitive<2,float>());
            mergedGrid.back().tris[0] = nonMortarTargetTris[i];
            mergedGrid.back().tris[1] = targetTri;

            for (int j=0; j<3; j++) {

                // Local coordinates in the domain triangle
                mergedGrid.back().localCoords[0][j] = cT.nodes[cPT.vertices(j)].domainPos();
            
                // Local coordinates in the target triangle
                mergedGrid.back().localCoords[1][j] = psurface->getLocalTargetCoords(GlobalNodeIdx(i, cPT.vertices(j)), targetTri);
            
                // world coordinates in the domain triangle
                mergedGrid.back().points[j] = 
                    PlaneParam<float>::linearInterpol<StaticVector<float,3> >(cT.nodes[cPT.vertices(j)].domainPos(),
                                                                              psurface->vertices(cT.vertices[0]), 
                                                                              psurface->vertices(cT.vertices[1]), 
                                                                              psurface->vertices(cT.vertices[2]));

            }            
            
        }
        
    }

}


void ContactToolBox::extractMergedGrid(PSurface<1,double>* psurface,
                                       std::vector<IntersectionPrimitive<1,float> >& mergedGrid)
{
    for (size_t i=0; i<psurface->domainSegments.size(); i++) {

        const PSurface<1,double>::DomainSegment&    cS = psurface->domainSegments[i];
        const std::vector<PSurface<1,double>::Node>& nodes = psurface->domainSegments[i].nodes;

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
            
            // Compute the world position of the overlap on the domain side */
            newOverlap.points[0][0] = psurface->vertices[cS.points[0]][0] * (1-nodes[j].domainLocalPosition)
                + psurface->vertices[cS.points[1]][0] * nodes[j].domainLocalPosition;
            newOverlap.points[0][1] = psurface->vertices[cS.points[0]][1] * (1-nodes[j].domainLocalPosition)
                + psurface->vertices[cS.points[1]][1] * nodes[j].domainLocalPosition;

            newOverlap.points[1][0] = psurface->vertices[cS.points[0]][0] * (1-cS.nodes[j+1].domainLocalPosition)
                + psurface->vertices[cS.points[1]][0] * cS.nodes[j+1].domainLocalPosition;
            newOverlap.points[1][1] = psurface->vertices[cS.points[0]][1] * (1-cS.nodes[j+1].domainLocalPosition)
                + psurface->vertices[cS.points[1]][1] * cS.nodes[j+1].domainLocalPosition;
            
            mergedGrid.push_back(newOverlap);
        }
        
    }

#if 0
    for (int i=0; i<mergedGrid.size(); i++)
        printf("overlap %d,   nonmortar (%g  -->  %g),    mortar (%g  -->  %g)\n", i,
               mergedGrid[i].localCoords[0][0][0], mergedGrid[i].localCoords[0][1][0],
               mergedGrid[i].localCoords[1][0][0], mergedGrid[i].localCoords[1][1][0]);
#endif
//     exit(0);
}

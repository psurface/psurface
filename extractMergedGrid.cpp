#include <psurface/Parametrization.h>
#include <psurface/ContactToolBox.h>
#include <vector>

//#define MY_DB

void ContactToolBox::extractMergedGrid(Parametrization* cPar,
                                       std::vector<IntersectionPrimitive<float> >& mergedGrid)
{
    if (cPar->getNumTriangles()==0)
        return;


    // ///////////////////////////////////////////////////
    // Set up point location structure
    // Can we use the routine in Parametrization ???
    // ///////////////////////////////////////////////////
    for (int i=0; i<cPar->getNumTriangles(); i++) {

        DomainTriangle& cT = cPar->triangles(i);

        cT.insertExtraEdges();

        for (int k=0; k<cT.nodes.size(); k++)
            cT.makeCyclicGeometrically(cT.nodes[k]);

        // the standard insertExtraEdges can miss edges if ghost nodes are present.  
        // We insert them now
        for (int k=0; k<3; k++){
            // size() returns an unsigned type, which underflows if edgePoints[k] is empty
            for (int l=0; l<((int)cT.edgePoints[k].size())-1; l++) {

                PlaneParam::DirectedEdgeIterator cE = cT.getDirectedEdgeIterator(cT.edgePoints[k][l], cT.edgePoints[k][l+1]);

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

    cPar->surface->computeTrianglesPerPoint();

    // get the array that relates the base grid triangles with the whole nonmortar surface
#if 0
    HxParameter* nonMortarTargetTriParam = cPar->params->find("targetTris");
    assert(nonMortarTargetTriParam);
    assert(nonMortarTargetTriParam->dim() == cPar->getNumTriangles());
    std::vector<int> nonMortarTargetTris(nonMortarTargetTriParam->dim());
    nonMortarTargetTriParam->getNum(&nonMortarTargetTris[0]);
#else
    std::vector<int> nonMortarTargetTris = cPar->domainSurfaceTriangleNumbers;
#endif

    //
    for (int i=0; i<cPar->getNumTriangles(); i++) {

        const DomainTriangle& cT = cPar->triangles(i);

        if (cT.nodes.size()<3)
            continue;

#ifdef PSURFACE_EXTRACT_ONLY_COMPLETELY_COVERED_FACES
        if (!isCompletelyCovered(cPar, i, &cT))
            continue;
#endif

        ////////////////////////////////
        PlaneParam::TriangleIterator cPT;
        for (cPT = cT.firstTriangle(); cPT.isValid(); ++cPT) {

            
            int targetTri = -1;
            try {
                targetTri = cPar->getImageSurfaceTriangle(i, cPT.vertices());
            } catch (Parametrization::ParamError){
                printf("exception caught!\n");
                targetTri = -1;
            }

            if (targetTri==-1)
                continue;
            
            assert(targetTri>=0 && targetTri<cPar->surface->triangles.size());

            std::tr1::array<IntersectionAlt, 3> intersections;
            for (int j=0; j<3; j++) {
                intersections[j].pos = cT.nodes[cPT.vertices(j)].domainPos();
                intersections[j].localTargetCoords = cPar->getLocalTargetCoords(GlobalNodeIdx(i, cPT.vertices(j)), 
                                                                                targetTri);
            }

            // //////////////////////////////////////////////
            // Assemble the triangles
            mergedGrid.push_back(IntersectionPrimitive<float>());
            mergedGrid.back().tris[0] = nonMortarTargetTris[i];
            mergedGrid.back().tris[1] = targetTri;
            
            mergedGrid.back().localCoords[0][0] = intersections[0].pos;
            mergedGrid.back().localCoords[0][1] = intersections[1].pos;
            mergedGrid.back().localCoords[0][2] = intersections[2].pos;
            
            mergedGrid.back().localCoords[1][0] = intersections[0].localTargetCoords;
            mergedGrid.back().localCoords[1][1] = intersections[1].localTargetCoords;
            mergedGrid.back().localCoords[1][2] = intersections[2].localTargetCoords;
            
            // world coordinate function
            mergedGrid.back().points[0] = 
                PlaneParam::linearInterpol<StaticVector<float,3> >(intersections[0].pos,
                                                    cPar->vertices(cT.vertices[0]), 
                                                    cPar->vertices(cT.vertices[1]), 
                                                    cPar->vertices(cT.vertices[2]));
            mergedGrid.back().points[1] = 
                PlaneParam::linearInterpol<StaticVector<float,3> >(intersections[1].pos, 
                                                    cPar->vertices(cT.vertices[0]), 
                                                    cPar->vertices(cT.vertices[1]),
                                                    cPar->vertices(cT.vertices[2]));
            mergedGrid.back().points[2] = 
                PlaneParam::linearInterpol<StaticVector<float,3> >(intersections[2].pos, 
                                                    cPar->vertices(cT.vertices[0]), 
                                                    cPar->vertices(cT.vertices[1]), 
                                                    cPar->vertices(cT.vertices[2]));           
            
            
        }
        
    }
}

#ifdef PSURFACE_EXTRACT_ONLY_COMPLETELY_COVERED_FACES
bool ContactToolBox::isCompletelyCovered(Parametrization* cPar, int tri, const DomainTriangle* cT)
{
    // Count number of CORNER/GHOST-nodes
    int nCorners = 0;
    for (int i=0; i<cT->nodes.size(); i++)
        if (cT->nodes[i].isCORNER_NODE() ||
            cT->nodes[i].isGHOST_NODE())
            nCorners++;

    // Count number of parametrization triangles
    PlaneParam::TriangleIterator cPT;
    int nTris = 0;

    for (cPT = cT->firstTriangle(); cPT.isValid(); ++cPT) {
        int targetTri = cPar->getImageSurfaceTriangle(tri, cPT.vertices());

        if (targetTri!=-1)
            nTris++;
    }

    // compare with expected number of triangles (by Euler's formula)
    return nCorners==3 && (nTris == 1 - cT->nodes.size() + cT->getNumEdges());
}
#endif

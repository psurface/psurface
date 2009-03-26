#include <psurface/Parametrization.h>
#include <psurface/ContactToolBox.h>
#include <vector>

//#define MY_DB

void ContactToolBox::extractMergedGrid(Parametrization* cPar,
                                       std::vector<IntersectionPrimitive>& mergedGrid)
{
    int i, j, k, l;

    if (cPar->getNumTriangles()==0)
        return;


    ////////////////////////////////////////////
    // Set up point location structure
    // Can we use the routine in Parametrization ???
    for (i=0; i<cPar->getNumTriangles(); i++) {

        DomainTriangle& cT = cPar->triangles(i);
        //             if (i==206)
        //                 ParamToolBox::display(cT, i);
        cT.insertExtraEdges();
        //             if (i==206)
        //                 ParamToolBox::display(cT, i);
        for (k=0; k<cT.nodes.size(); k++)
            cT.makeCyclicGeometrically(cT.nodes[k]);

        // the standard insertExtraEdges can miss edges if ghost nodes are present.  
        // We insert them now
        for (k=0; k<3; k++){
            for (l=0; l<cT.edgePoints[k].size()-1; l++) {

                PlaneParam::DirectedEdgeIterator cE = cT.getDirectedEdgeIterator(cT.edgePoints[k][l], cT.edgePoints[k][l+1]);

                if (cE.isValid() && cE.getDPrev().from() != cE.getONext().to())
                    cT.addEdge(cE.getONext().to(), cE.to(), true);
                
            }
            
        }

        // and since we have added more edges, we have to redo the cyclic ordering
        for (k=0; k<cT.nodes.size(); k++)
            cT.makeCyclicGeometrically(cT.nodes[k]);
        
        //            if (i==206)
        //                 ParamToolBox::display(cT, i);
        
        //
        for (k=0; k<3; k++){
            for (l=0; l<cT.edgePoints[k].size(); l++) {
                
                if (!cT.nodes[cT.edgePoints[k][l]].isOnCorner()) {
                    cT.nodes[cT.edgePoints[k][l]].setDomainEdge(k);
                    cT.nodes[cT.edgePoints[k][l]].setDomainEdgePosition(l);
                }
            }
            
        }
        
    }
//     return;
    cPar->surface->computeTrianglesPerPoint();

    // get the array that relates the base grid triangles with the whole nonmortar surface
    HxParameter* nonMortarTargetTriParam = cPar->params->find("targetTris");
    assert(nonMortarTargetTriParam);
    assert(nonMortarTargetTriParam->dim() == cPar->getNumTriangles());
    std::vector<int> nonMortarTargetTris(nonMortarTargetTriParam->dim());
    nonMortarTargetTriParam->getNum(&nonMortarTargetTris[0]);

    //
#ifndef MY_DB
    for (i=0; i<cPar->getNumTriangles(); i++) {
#else
    for (i=0; i<1; i++) {
#endif   
        const DomainTriangle& cT = cPar->triangles(i);
#ifdef MY_DB
        printf("merging triangle %d,  %d nodes\n", i, cT.nodes.size());
#endif
        if (cT.nodes.size()<3)
            continue;

#ifdef MY_DB
        ParamToolBox::display(cT, i);
        cT.print(true, true, true);
#endif

        if (!isCompletelyCovered(cPar, i, &cT))
            continue;

        //break;
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

            McSArray<IntersectionAlt, 3> intersections;
            for (j=0; j<3; j++) {
                intersections[j].pos = cT.nodes[cPT.vertices(j)].domainPos();
                intersections[j].localTargetCoords = cPar->getLocalTargetCoords(GlobalNodeIdx(i, cPT.vertices(j)), 
                                                                                targetTri);
            }

            // //////////////////////////////////////////////
            // Assemble the triangles
            mergedGrid.push_back(IntersectionPrimitive());
            mergedGrid.back().tris[0] = nonMortarTargetTris[i];
            mergedGrid.back().tris[1] = targetTri;
            
            mergedGrid.back().localCoords[0][0] = intersections[0].pos;
            mergedGrid.back().localCoords[0][1] = intersections[1].pos;
            mergedGrid.back().localCoords[0][2] = intersections[2].pos;
            
            mergedGrid.back().localCoords[1][0] = intersections[0].localTargetCoords;
            mergedGrid.back().localCoords[1][1] = intersections[1].localTargetCoords;
            mergedGrid.back().localCoords[1][2] = intersections[2].localTargetCoords;
            
#ifndef MY_DB
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
            
            
#else
            // output in barycentric coordinates (for debugging)
            mergedGrid[newTri].points[0] = StaticVector<float,3>(intersections[0].pos.x, intersections[0].pos.y, 0);
            mergedGrid[newTri].points[1] = StaticVector<float,3>(intersections[1].pos.x, intersections[1].pos.y, 0);
            mergedGrid[newTri].points[2] = StaticVector<float,3>(intersections[2].pos.x, intersections[2].pos.y, 0);
#endif
            
            
        }
        
    }
}

bool ContactToolBox::isCompletelyCovered(Parametrization* cPar, int tri, const DomainTriangle* cT)
{
    // Count number of CORNER/GHOST-nodes
    int nCorners = 0;
    for (int i=0; i<cT->nodes.size(); i++)
        if (cT->nodes[i].isCORNER_NODE() ||
            cT->nodes[i].isGHOST_NODE())
            nCorners++;

    // debug count
    static int count = 0;
    // Count number of parametrization triangles
    PlaneParam::TriangleIterator cPT;
    int nTris = 0;

    for (cPT = cT->firstTriangle(); cPT.isValid(); ++cPT) {
        int targetTri = cPar->getImageSurfaceTriangle(tri, cPT.vertices());

        if (targetTri!=-1)
            nTris++;
    }

    count++;

    //     printf("nCorners; %d   True number of triangles: %d,    expected: %d\n",
    //            nCorners, nTris, 1 - cT->nodes.size() + cT->getNumEdges());
    
    // compare with expected number of triangles (by Euler's formula)
    return nCorners==3 && (nTris == 1 - cT->nodes.size() + cT->getNumEdges());
}

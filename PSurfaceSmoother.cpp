#include <psurface/PSurface.h>
#include "PSurfaceSmoother.h"

#include <psurface/CircularPatch.h>
#include <psurface/DomainPolygon.h>
#include <psurface/SparseMatrix.h>


template <class ctype>
void PSurfaceSmoother<ctype>::applyEdgeRelaxation(PSurface<2,ctype>* psurface, int edge, 
                                                  bool keepPatches, std::vector<unsigned int>& nodeStack)
{
        
    McEdge& cE = psurface->edges(edge);

    if (cE.triangles.size()!=2)
        return;
        
    StaticVector<float,2> quadCoords[4];
    DomainPolygon quadri(psurface);
    bool flipped;

    psurface->triangles(cE.triangles[0]).checkConsistency("PreRelax0");
    psurface->triangles(cE.triangles[1]).checkConsistency("PreRelax1");

    ParamToolBox::mergeTwoTrianglesIntoQuadrangle(cE.triangles[0], cE.triangles[1], 
                                                  quadri, flipped, quadCoords, nodeStack, psurface);
        
    // apply the desired parametrization

    if (psurface->triangles(cE.triangles[0]).patch != psurface->triangles(cE.triangles[1]).patch && 
        keepPatches) {
            applyHorizontalRelaxation(quadri);
        } else
            quadri.applyParametrization();
        
    // undo the merge
    CircularPatch<float> cutter(2, psurface);
        
    // all this true copying is rather inefficient...
    std::tr1::array<DomainTriangle<float>, 2> backupTriangles;
    backupTriangles[0] = psurface->triangles(cE.triangles[0]);
    backupTriangles[1] = psurface->triangles(cE.triangles[1]);

    cutter[0] = cE.triangles[0];
    cutter[1] = cE.triangles[1];

    // this is necessary because if quadri.triangulate fails the nodeStack is corrupted
    std::vector<unsigned int> tempNodeStack(nodeStack);

    if (!quadri.triangulate(cutter, tempNodeStack)){
            
        std::cerr << "Couldn't cut quadrangle -- aborting" << std::endl;

        psurface->triangles(cE.triangles[0]) = backupTriangles[0];
        psurface->triangles(cE.triangles[1]) = backupTriangles[1];

        return;
    }

    nodeStack = tempNodeStack;

    if (flipped)
        psurface->triangles(cutter[1]).flip();

    if (psurface->triangles(cE.triangles[0]).patch != psurface->triangles(cE.triangles[1]).patch && 
        keepPatches) {
        psurface->triangles(cE.triangles[0]).applyParametrization(psurface->iPos);
        psurface->triangles(cE.triangles[1]).applyParametrization(psurface->iPos);
    }

    psurface->triangles(cE.triangles[0]).checkConsistency("PostRelax0");
    psurface->triangles(cE.triangles[1]).checkConsistency("PostRelax1");

    psurface->integrateTriangle(cE.triangles[0]);
    psurface->integrateTriangle(cE.triangles[1]);

}

template <class ctype>
void PSurfaceSmoother<ctype>::applyHorizontalRelaxation(DomainPolygon& quadri, PSurface<2,ctype>* psurface)
{
    // compute lambdas
    SparseMatrix<float> lambda_ij(quadri.nodes.size());

    quadri.computeFloaterLambdas(lambda_ij, psurface->iPos);

    // build matrix 
    lambda_ij *= -1;

    for (int i=0; i<lambda_ij.nRows(); i++)
        lambda_ij.setEntry(i, i, 1);

    // compute the right side, only x-coordinates are interesting
    std::vector<float> b_x(quadri.nodes.size());

    std::fill(b_x.begin(), b_x.end(), 0);

    for (int i=0; i<quadri.nodes.size(); i++)
        if (!quadri.nodes[i].isINTERIOR_NODE())
            b_x[i] = quadri.nodes[i].domainPos()[0];

    // solve the system
    int maxIter=3000;
    std::vector<float> residual;
    std::vector<float> result = b_x;

    // xCoords
    for (int i=0; i<quadri.nodes.size(); i++)
        result[i] = quadri.nodes[i].domainPos()[0];

#warning Smoothing linear system is not actually solved!
    //lambda_ij.BiCGSTAB(b_x, result, residual, &maxIter, 0.000001);

    for (size_t i=0; i<quadri.nodes.size(); i++)
        if (quadri.nodes[i].isINTERIOR_NODE())
            //quadri.nodes[i].domainPos.x = result[i];
            quadri.nodes[i].setDomainPos(StaticVector<float,2>(result[i], quadri.nodes[i].domainPos()[1]));

}

template <class ctype>
void PSurfaceSmoother<ctype>::applyVertexRelaxation()
{
#if 0
    const int numVertices = par->getNumVertices();

    theWorkArea->startWorking(numVertices, "");

    int vCounter = 0;

    DomainVertex* centerPoint;
    // loop over all regular vertices
    for (centerPoint=par->vertices.first(); centerPoint; 
         centerPoint=par->vertices.succ(centerPoint), vCounter++) {

        theWorkArea->progressStep();

        if (!(vCounter%100) && theWorkArea->wasInterrupted())
            break;
        
        int i, j, k;
        
        for (DomainTriangle* cT=par->triangles.first(); cT; cT=par->triangles.succ(cT))
            cT->checkConsistency("Beginning of relaxation Loop\n");

        ////////////////////////////////////////
        // merge the star of the current node into one DomainPolygon

        DomainPolygon fullStar;
        int newCenterNode = -1;
        McDArray<DomainTriangle*> fullStarTris(0);

        if (!ParamToolBox::mergeStarIntoPolygon(centerPoint, fullStar, fullStarTris, newCenterNode))
            continue;


        ////////////////////////////////////////
        // apply the requested parametrization

        //ParamToolBox::display(fullStar, vCounter);

        if (portRadius.getValue() > 0.995) {
            
            fullStar.applyParametrization();

        } else {

            //////////////////////////////////////////////////////////////////////
            // apply parametrization to only a part of the graph

            // find shortest edge length
            float minLength = FLT_MAX;
            for (i=0; i<fullStar.boundaryPoints.size(); i++)
                if (fullStar.nodes[fullStar.cornerNode(i)].domainPos.length() < minLength)
                    minLength = centerPoint->edges[i]->length();

            float freeRadius = minLength*portRadius.getValue();

            McDArray<int> interiorNodes(0);
            McDArray<int> boundaryNodes(0);
            int cN;
            
            for (cN=0; cN<fullStar.nodes.size(); cN++)
                if (fullStar.nodes[cN].domainPos.length() < freeRadius)         
                    interiorNodes.append(cN);
                else
                    boundaryNodes.append(cN);
            
            assert(boundaryNodes.size()>=3);
            
            // compute lambdas
            McSparseMatrix<float, false> lambda_ij(interiorNodes.size()+boundaryNodes.size());
            
            fullStar.computeFloaterLambdas(lambda_ij, interiorNodes, boundaryNodes);
            
            
            // build matrix 
            lambda_ij *= -1;
            
            for (i=0; i<lambda_ij.nRows(); i++)
                lambda_ij.setEntry(i, i, 1);
            
            
            // compute the right side, split in x and y coordinates
            // this doesn't really use the sparse matrix well
            
            McDArray<float> b_x(interiorNodes.size()+boundaryNodes.size());
            McDArray<float> b_y(interiorNodes.size()+boundaryNodes.size());
            
            b_x.fill(0);
            b_y.fill(0);
            
            for (j=0; j<boundaryNodes.size(); j++)
                b_x[interiorNodes.size()+j] = fullStar.nodes[boundaryNodes[j]].domainPos.x;
            
            for (j=0; j<boundaryNodes.size(); j++)
                b_y[interiorNodes.size()+j] = fullStar.nodes[boundaryNodes[j]].domainPos.y;
            
            // solve the system
            int maxIter=3000;
            McDArray<float> residue;
            McDArray<float> result = b_x;
            
            // xCoords
            for (i=0; i<interiorNodes.size(); i++)
                result[i] = fullStar.nodes[interiorNodes[i]].domainPos.x;
            
            lambda_ij.SOR(b_x, result, residue, &maxIter, 0.000001, 0.9);
            
            for (i=0; i<interiorNodes.size(); i++)
                fullStar.nodes[interiorNodes[i]].domainPos.x = result[i];
            
            // yCoords
            maxIter = 3000;
            result = b_y;
            for (i=0; i<interiorNodes.size(); i++)
                result[i] = fullStar.nodes[interiorNodes[i]].domainPos.y;
            
            lambda_ij.SOR(b_y, result, residue, &maxIter, 0.000001, 0.9);
            
            for (i=0; i<interiorNodes.size(); i++)
                fullStar.nodes[interiorNodes[i]].domainPos.y = result[i];
            
        }


        //ParamToolBox::display(fullStar, 20+vCounter);

        //////////////////////////////////////////////////////////////////////
        // look for the nodes that is closest to (0, 0).  It will become
        // the new centerPoint

        if (!portKeepPatches.getValue(1)) {
            int cN;
            float minDist = FLT_MAX;
            
            for (cN=0; cN<fullStar.nodes.size(); cN++){
                
                if (fullStar.nodes[cN].domainPos.length2() < minDist){
                    newCenterNode = cN;
                    minDist = fullStar.nodes[cN].domainPos.length2();
                }
            }
        }

        if (fullStar.nodes[newCenterNode].type != PlaneParam::Node::INTERIOR_NODE)
            printf("Warning:  New centernode is not INTERIOR_NODE!\n");

        //////////////////////////////////////////////////////////////////////
        // recut the star

        // we first do a cut the polygon into 'pizza slices'

        *centerPoint = fullStar.nodes[newCenterNode].imagePos;

        for (i=0; i<fullStarTris.size(); i++){
            fullStar.slice(newCenterNode, centerPoint, i*3);
            fullStar.checkConsistency("Slicing");
            //display(fullStar, i+2);
        }

        //ParamToolBox::display(fullStar, 40+vCounter);
        
        // the polygon has been cut.  Move each slice to its
        // original triangle

        int cN;

        for (i=0; i<fullStarTris.size(); i++) {

            for (cN=0; cN<fullStar.nodes.size(); cN++)
                fullStar.nodes[cN].location = PlaneParam::Node::IN_POLYGON;

            DomainTriangle* cT = fullStarTris[i];
            
            int offset=0;
            if (cT->points[0]==centerPoint)
                offset = 0;
            else if (cT->points[1]==centerPoint)
                offset = 1;
            else
                offset = 2;
        

            // copy edgePoint arrays
            for (j=0; j<3; j++){
                cT->edgePoints[(j+offset)%3] = fullStar.edgePoints[(3*i+1+j)%fullStar.edgePoints.size()];
                fullStar.edgePoints[(3*i+1+j)%fullStar.edgePoints.size()].clear();
            }

            // copy nodes using a graph-search algorithm
            for (j=0; j<3; j++){
                for (k=0; k<cT->edgePoints[j].size(); k++)
                    if (fullStar.nodes[cT->edgePoints[j][k]].location != PlaneParam::Node::IN_TRIANGLE)
                        moveSubGraph(cT->edgePoints[j][k], fullStar, newCenterNode);
            }

            // make a copy of the centerNode
            fullStar.nodes.appendSpace(1);
            int localCenterNode = fullStar.nodes.size()-1;
            fullStar.nodes[localCenterNode].setValue(fullStar.nodes[newCenterNode].domainPos, 
                                                     fullStar.nodes[newCenterNode].imagePos, 
                                                     PlaneParam::Node::CORNER_NODE);

            fullStar.nodes[localCenterNode].location = PlaneParam::Node::IN_TRIANGLE;

            for (j=fullStar.nodes[newCenterNode].degree()-1; j>=0; j--)
                if (fullStar.nodes[fullStar.nodes[newCenterNode].neighbors[j]].location == PlaneParam::Node::IN_TRIANGLE) {
                    fullStar.nodes[localCenterNode].neighbors.append(fullStar.nodes[newCenterNode].neighbors[j]);
                    fullStar.nodes[fullStar.nodes[newCenterNode].neighbors[j]].replaceReferenceTo(newCenterNode, localCenterNode);
                    fullStar.nodes[newCenterNode].neighbors.remove(j);
                }

            cT->edgePoints[0+offset][0] = localCenterNode;
            cT->edgePoints[(2+offset)%3].last() = localCenterNode;

            //////////////////////////////////////
            // sort out the nodes that belong onto the triangle
            int numTriNodes = 0;
            int triNode;

            for (triNode=0; triNode<fullStar.nodes.size(); triNode++) 
                if (fullStar.nodes[triNode].location == PlaneParam::Node::IN_TRIANGLE)
                    numTriNodes++;
            
            int triCount = 0;
            cT->nodes.resize(numTriNodes);
            McDArray<int> offArr(fullStar.nodes.size());
            
            for (triNode=0; triNode<fullStar.nodes.size(); triNode++)
                if (fullStar.nodes[triNode].location == PlaneParam::Node::IN_TRIANGLE) {
                    cT->nodes[triCount] = fullStar.nodes[triNode];
                    fullStar.invalidate(triNode);
                    offArr[triNode] = triCount;
                    triCount++;
                }
            
            for (j=0; j<numTriNodes; j++)
                for (k=0; k<cT->nodes[j].neighbors.size(); k++)
                    cT->nodes[j].neighbors[k] = offArr[cT->nodes[j].neighbors[k]];
            
            for (j=0; j<3; j++)
                for (k=0; k<cT->edgePoints[j].size(); k++)
                    cT->edgePoints[j][k] = offArr[cT->edgePoints[j][k]];



            /////////////////////////////////////
            cT->checkConsistency("After Vertex Relaxation\n");

            cT->installBarycentricCoordinates();


            fullStar.checkConsistency("before garbage collection\n");

            // reuse offArr
            fullStar.garbageCollection(offArr);
            newCenterNode -= offArr[newCenterNode];

            fullStar.checkConsistency("after garbage collection\n");


        }

    }

    theWorkArea->stopWorking();
 
#endif   
}

template <class ctype>
void PSurfaceSmoother<ctype>::moveSubGraph(int startingNode, DomainPolygon& from, int centerNode)
{
#if 0
    if (startingNode==centerNode)
        return;

    from.nodes[startingNode].location = Node::IN_TRIANGLE;

    for (int i=0; i<from.nodes[startingNode].degree(); i++)
        if (from.nodes[from.nodes[startingNode].neighbors(i)].location!=Node::IN_TRIANGLE)
            moveSubGraph(from.nodes[startingNode].neighbors(i), from, centerNode);
#endif
}

#include <psurface/ContactToolBox.h>
#include <psurface/ContactBoundary.h>
#include <psurface/NormalProjector.h>

#include <mclib/McSArray.h>

#include <mclib/McOctree.h>
#include "MyMcVec3f.h"

// Debug
// #include <hxcluster/HxCluster.h>
// #include <Amira/HxObjectPool.h>

void ContactToolBox::buildContactSurface(Parametrization* cPar, 
                                         const Surface* surf1,  const Surface* surf2,
                                         float epsilon, 
                                         void (*obsDirections)(const double* pos, double* dir))
{
    int i;
    
    // set up parametrization
    cPar->surface = const_cast<Surface*>(surf2);
    cPar->patches.resize(1);
    cPar->patches[0].innerRegion = 0;
    cPar->patches[0].outerRegion = 1;
    cPar->patches[0].boundaryId  = 0;
            
    // ///////
    McSArray<ContactBoundary, 2> contactBoundary;
    
    const_cast<Surface*>(surf1)->removeUnusedPoints();
    const_cast<Surface*>(surf2)->removeUnusedPoints();
    
    contactBoundary[0].init(surf1);
    contactBoundary[1].init(surf2);
    
    contactOracle(surf1, surf2, contactBoundary[0].vertices, contactBoundary[1].vertices, epsilon);
    
    if (!contactBoundary[0].vertices.size()) {
        
        printf("No contact surface found!\n");
        return;
        
    }
    
    printf("%ld resp. %ld contact nodes found!\n", contactBoundary[0].vertices.size(),
           contactBoundary[1].vertices.size());


    // //////////////////
    // BEGIN DEBUG
    // //////////////////
#if 0
    char nodesName[50];
    // display nodes
    sprintf(nodesName, "Poly.0.nodes");
    HxCluster* theNodes = new HxCluster();
    theNodes->setLabel(nodesName);

    theNodes->points.resize(contactBoundary[0].vertices.size());
    theNodes->ids.resize(contactBoundary[0].vertices.size());

    theNodes->dataColumns.resize(0);

    for (i=0; i<contactBoundary[0].vertices.size(); i++) {
        theNodes->points[i] = McVec3f(contactBoundary[0].surf->points[contactBoundary[0].vertices[i]].x,
                                      contactBoundary[0].surf->points[contactBoundary[0].vertices[i]].y,
                                      contactBoundary[0].surf->points[contactBoundary[0].vertices[i]].z);

        theNodes->ids[i] = 0;
    }

    theNodes->writePSI(nodesName);
    theObjectPool->addObject(theNodes);

    // display nodes
    sprintf(nodesName, "Poly.1.nodes");
    theNodes = new HxCluster();
    theNodes->setLabel(nodesName);

    theNodes->points.resize(contactBoundary[1].vertices.size());
    theNodes->ids.resize(contactBoundary[1].vertices.size());

    theNodes->dataColumns.resize(0);

    for (i=0; i<contactBoundary[1].vertices.size(); i++) {
        theNodes->points[i] = McVec3f(contactBoundary[1].surf->points[contactBoundary[1].vertices[i]].x,
                                      contactBoundary[1].surf->points[contactBoundary[1].vertices[i]].y,
                                      contactBoundary[1].surf->points[contactBoundary[1].vertices[i]].z);

        theNodes->ids[i] = 0;
    }

    theNodes->writePSI(nodesName);
    theObjectPool->addObject(theNodes);
#endif
    // ///////////////////////////////////
    // END DEBUG
    // //////////////////////////////////

    // create the two contact patches
    
    computeContactPatch(surf1, contactBoundary[0]);
    computeContactPatch(surf2, contactBoundary[1]);
    
    std::cout << "Contact patches contain " << contactBoundary[0].triIdx.size() 
              << " (resp. " << contactBoundary[1].triIdx.size() << ") triangles." << std::endl;
    
    // the nonmortar side becomes the base grid of the parametrization
    for (i=0; i<contactBoundary[0].vertices.size(); i++)
        cPar->newVertex(surf1->points[contactBoundary[0].vertices[i]]);
    
    cPar->params->insert(new HxParameter("targetTris",
                                         contactBoundary[0].triIdx.size(),
                                         contactBoundary[0].triIdx.dataPtr()));
    
    McDArray<int> vertexOffsets = contactBoundary[0].getVertexOffsets();
    for (i=0; i<contactBoundary[0].triIdx.size(); i++) {
        
        TriangleIdx newTri = cPar->createSpaceForTriangle(vertexOffsets[contactBoundary[0].triangles(i).points[0]],
                                                          vertexOffsets[contactBoundary[0].triangles(i).points[1]],
                                                          vertexOffsets[contactBoundary[0].triangles(i).points[2]]);
        cPar->integrateTriangle(newTri);
        cPar->triangles(newTri).patch = 0;
        
    }
    
    // compute projection
    NormalProjector* projector = new NormalProjector;
    
    projector->handleSide(cPar, contactBoundary[1], obsDirections);
    
    delete projector;

}


void ContactToolBox::contactOracle(const Surface* surf1, const Surface* surf2,
                                   McDArray<int>& contactNodes1, McDArray<int>& contactNodes2,
                                   float epsilon)
{
    int i, j;
    const float epsSquared = epsilon*epsilon;
    
    McBox3f bbox1, bbox2;
    surf1->getBoundingBox(&bbox1[0]);
    surf2->getBoundingBox(&bbox2[0]);
    
    bbox1.extendByEps(epsilon);
    bbox2.extendByEps(epsilon);

    // The possible contact patches must be in intersectBox
    McBox3f intersectBox = bbox1.intersectWith(bbox2);

    // We first put the vertices of surface1 into an octree
    McOctree<MyMcVec3f> octree1(bbox1);
    McDArray<MyMcVec3f> points1(surf1->points.size());        

    for (int i=0; i<surf1->points.size(); i++) {
        points1[i] = surf1->points[i];
        octree1.insert(points1[i]);
    }

    octree1.enableUniqueLookup(points1.size(), points1.dataPtr());
    
    // We first put the vertices of surface2 into an octree
    McOctree<MyMcVec3f> octree2(intersectBox);
    McDArray<MyMcVec3f> points2(surf2->points.size());        

    for (int i=0; i<surf2->points.size(); i++){
        
        points2[i] = surf2->points[i];
        if (intersectBox.contains(surf2->points[i]))
            octree2.insert(points2[i]);
        
    }
    
    octree2.enableUniqueLookup(points2.size(), points2.dataPtr());
    
    // Two bitfields to mark the contact nodes
    std::vector<bool> contactField2(surf2->points.size(), false);
    
    // ///////////////////////////////////////////////////////////////////////
    //   Completely keep the domain (nonmortar) side, and make
    //   each vertex on the mortar side be part of the contact surface if it
    //   is less than epsilon away from any nonmortar triangle.
    // ///////////////////////////////////////////////////////////////////////


    // ///////////////////////////////////////////////////////////////////////
    //   Loop over all vertices in surface2 and check whether it is no more
    //   than epsilon away from surface1
    // ///////////////////////////////////////////////////////////////////////

    //  Loop over all triangles in surface1
    for (int i=0; i<surf1->triangles.size(); i++) {

        const McVec3f& p0 = surf1->points[surf1->triangles[i].points[0]];
        const McVec3f& p1 = surf1->points[surf1->triangles[i].points[1]];
        const McVec3f& p2 = surf1->points[surf1->triangles[i].points[2]];

        //  Look up the octree for points in a conversative neighborhood
        //  of the triangle.  The triangle's boundingbox + epsilon will do
        McBox3f queryBox(p0, p1);
        queryBox.extendBy(p2);
        queryBox.extendByEps(epsilon);

        McDArray<int> result;
        octree2.lookupIndex(queryBox, result);
                                       
        for (j=0; j<result.size(); j++) {

            // Don't recompute everything is the vertex is already marked
            if (contactField2[result[j]])
                continue;

            const McVec3f& candidatePoint = surf2->points[result[j]];

            McVec3f q = getClosestPointOnTriangle(p0, p1, p2, candidatePoint);

            if ( (q-candidatePoint).length2() < epsSquared)
                contactField2[result[j]] = true;

        }
        
    }

    //  Loop over all triangles in surface2
    for (int i=0; i<surf2->triangles.size(); i++) {

        const McVec3f& p0 = surf2->points[surf2->triangles[i].points[0]];
        const McVec3f& p1 = surf2->points[surf2->triangles[i].points[1]];
        const McVec3f& p2 = surf2->points[surf2->triangles[i].points[2]];

        //  Look up the octree for points in a conversative neighborhood
        //  of the triangle.  The triangle's boundingbox + epsilon will do
        McBox3f queryBox(p0, p1);
        queryBox.extendBy(p2);
        queryBox.extendByEps(epsilon);

        McDArray<int> result;
        octree1.lookupIndex(queryBox, result);

        // If the bounding box contains any vertices from surface one we keep
        // the whole triangle
        if (result.size() > 0) 
            for (int j=0; j<3; j++)
                contactField2[surf2->triangles[i].points[j]] = true;

    }



    //  All vertices of surface1 belong to the result contact surface

    int c = 0;

    contactNodes1.resize(surf1->points.size());
    for (i=0; i<surf1->points.size(); i++)
        contactNodes1[c++] = i;

    c = 0;
    int nSetBits2 = 0;
    for (i=0; i<contactField2.size(); i++)
        if (contactField2[i])
            nSetBits2++;

    contactNodes2.resize(nSetBits2);
    for (i=0; i<contactField2.size(); i++)
        if (contactField2[i])
            contactNodes2[c++] = i;

}

McVec3f ContactToolBox::getClosestPointOnTriangle(const McVec3f& p0,
                                                  const McVec3f& p1,
                                                  const McVec3f& p2,          
                                                  const McVec3f& candidate)
{
        
    // local base
    McVec3f a = p1 - p0;
    McVec3f b = p2 - p0;
    McVec3f c = a.cross(b);
    c.normalize();
        
    McVec3f x = candidate - p0;
        
    // write x in the new base  (Cramer's rule)
    //McMat3f numerator(a, b, c);
    double denominatorDet = McMat3f(a, b, c).det();
    McMat3f alphaMat(x, b, c);
    McMat3f betaMat(a, x, c);
    McMat3f gammaMat(a, b, x);
    
    float alpha = alphaMat.det()/denominatorDet;
    float beta  = betaMat.det()/denominatorDet;
    //float gamma = gammaMat.det()/denominatorDet;
    
    // check whether orthogonal projection onto the ab plane is in triangle
    if (alpha>=0 && beta>=0 && (1-alpha-beta)>=0) {
        // The orthogonal projection of the candidate point on the plane
        // spanned by the triangle
        return p0 + alpha*a + beta*b;
    }
    
    // ////////////////////////////////////////////////////////////////////////////////
    //   The candidate point is not 'over' the triangle.  Then the close point to it
    //   _on_ the triangle is the closest point on the boundary of the triangle.
    // ////////////////////////////////////////////////////////////////////////////////

    float bestDist = std::numeric_limits<float>::max();
    McVec3f bestPoint;

    // I need the points in an array
    McVec3f points[3];
    points[0] = p0;
    points[1] = p1;
    points[2] = p2;
    
    // check point against edges
    for (int i=0; i<3; i++){
        
        McVec3f from = points[i];
        McVec3f to   = points[(i+1)%3];
        
        McVec3f edge = to - from;
        
        float projectLength = edge.dot(candidate - from)/edge.length();
        McVec3f projection = edge/edge.length() * projectLength;
        
        float orthoDist = ((candidate-from) - projection).length();
        
        if (projectLength>=0 && projectLength<=edge.length() && orthoDist<bestDist) {
            bestDist = orthoDist;
            bestPoint = projection + from;
        }
    }
    
    // check point against vertices
    for (int i=0; i<3; i++){
        float dist = (candidate - points[i]).length();
        if (dist < bestDist){
            bestDist = dist;
            bestPoint = points[i];
        }
    }
    
    
    return bestPoint;
}


void ContactToolBox::computeContactPatch(const Surface* surf, ContactBoundary& cBound)
{
    int i;
    
    //////////////////////////////////////////////////
    // create triangles
    cBound.vertices.sort(&mcStandardCompare);

    for (i=0; i<surf->triangles.size(); i++) {
        
        const McVec3i& p = surf->triangles[i].points;
        
        if (cBound.vertices.findSorted(p[0], mcStandardCompare)>=0 &&
            cBound.vertices.findSorted(p[1], mcStandardCompare)>=0 &&
            cBound.vertices.findSorted(p[2], mcStandardCompare)>=0) {
            
            cBound.triIdx.append(i);
        }
    }
}

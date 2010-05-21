#include <tr1/array>
#include <vector>
#include <set>

#include <psurface/ContactToolBox.h>
#include <psurface/ContactBoundary.h>
#include <psurface/NormalProjector.h>
#include <psurface/MultiDimOctree.h>
#include <psurface/PointIntersectionFunctor.h>

template <class ctype>
void ContactToolBox<ctype>::buildContactSurface(PSurface<2,ctype>* cPar, 
                                                const Surface* surf1,  const Surface* surf2,
                                                const DirectionFunction<3,ctype>* domainDirection,
                                                const DirectionFunction<3,ctype>* targetDirection
                                                )
{
    // set up parametrization
    cPar->surface = const_cast<Surface*>(surf2);
    cPar->patches.resize(1);
    cPar->patches[0].innerRegion = 0;
    cPar->patches[0].outerRegion = 1;
    cPar->patches[0].boundaryId  = 0;
            
    // ///////
    std::tr1::array<ContactBoundary, 2> contactBoundary;
    
    const_cast<Surface*>(surf1)->removeUnusedPoints();
    const_cast<Surface*>(surf2)->removeUnusedPoints();
    
    contactBoundary[0].init(surf1);
    contactBoundary[1].init(surf2);
    
    contactBoundary[0].vertices.resize(surf1->points.size());
    for (int i=0; i<surf1->points.size(); i++)
        contactBoundary[0].vertices[i] = i;

    contactBoundary[1].vertices.resize(surf2->points.size());
    for (int i=0; i<surf2->points.size(); i++)
        contactBoundary[1].vertices[i] = i;

    if (!contactBoundary[0].vertices.size()) {
        
        printf("No contact surface found!\n");
        return;
        
    }
    
    std::cout << contactBoundary[0].vertices.size() << " resp. "
              << contactBoundary[1].vertices.size() << " contact nodes found!" << std::endl;

    // create the two contact patches
    
    computeContactPatch(surf1, contactBoundary[0]);
    computeContactPatch(surf2, contactBoundary[1]);
    
    std::cout << "Contact patches contain " << contactBoundary[0].triIdx.size() 
              << " (resp. " << contactBoundary[1].triIdx.size() << ") triangles." << std::endl;
    
    // the nonmortar side becomes the base grid of the parametrization
    for (size_t i=0; i<contactBoundary[0].vertices.size(); i++) {
        StaticVector<ctype,3> newVertex;
        for (int j=0; j<3; j++)
            newVertex[j] = surf1->points[contactBoundary[0].vertices[i]][j];
        cPar->newVertex(newVertex);
    }
    
    cPar->domainSurfaceTriangleNumbers = contactBoundary[0].triIdx;

    std::vector<int> vertexOffsets = contactBoundary[0].getVertexOffsets();
    for (size_t i=0; i<contactBoundary[0].triIdx.size(); i++) {
        
        int newTri = cPar->createSpaceForTriangle(vertexOffsets[contactBoundary[0].triangles(i).points[0]],
                                                  vertexOffsets[contactBoundary[0].triangles(i).points[1]],
                                                  vertexOffsets[contactBoundary[0].triangles(i).points[2]]);
        cPar->integrateTriangle(newTri);
        cPar->triangles(newTri).patch = 0;
        
    }
    
    // compute projection
    NormalProjector<ctype> projector(cPar);
    
    projector.project(contactBoundary[1], domainDirection, targetDirection);

}


template <class ctype>
void ContactToolBox<ctype>::computeContactPatch(const Surface* surf, ContactBoundary& cBound)
{
    //////////////////////////////////////////////////
    // create triangles
    //////////////////////////////////////////////////

    /** \todo Maybe cbound.vertices can be a std::set?
        Then we wouldn't have to copy
    */
    std::set<int> vertexSet(cBound.vertices.begin(), cBound.vertices.end());

    for (int i=0; i<surf->triangles.size(); i++) {
        
        if (vertexSet.find(surf->triangles[i].points[0]) != vertexSet.end() &&
            vertexSet.find(surf->triangles[i].points[1]) != vertexSet.end() &&
            vertexSet.find(surf->triangles[i].points[2]) != vertexSet.end()) {
            
            cBound.triIdx.push_back(i);
        }
    }
}

// ////////////////////////////////////////////////////////
//   Explicit template instantiations.
//   If you need more, you can add them here.
// ////////////////////////////////////////////////////////

template class ContactToolBox<float>;
template class ContactToolBox<double>;

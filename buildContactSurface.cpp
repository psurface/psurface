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
    const_cast<Surface*>(surf1)->removeUnusedPoints();
    const_cast<Surface*>(surf2)->removeUnusedPoints();
    
    std::cout << surf1->points.size() << " resp. "
              << surf2->points.size() << " contact nodes found!" << std::endl;

    std::cout << "Contact patches contain " << surf1->triangles.size() 
              << " (resp. " << surf2->triangles.size() << ") triangles." << std::endl;
    
    // the nonmortar side becomes the base grid of the parametrization
    for (size_t i=0; i<surf1->points.size(); i++) {
        StaticVector<ctype,3> newVertex;
        for (int j=0; j<3; j++)
            newVertex[j] = surf1->points[i][j];
        cPar->newVertex(newVertex);
    }
    
    for (size_t i=0; i<surf1->triangles.size(); i++) {
        
        int newTri = cPar->createSpaceForTriangle(surf1->triangles[i].points[0],
                                                  surf1->triangles[i].points[1],
                                                  surf1->triangles[i].points[2]);
        cPar->integrateTriangle(newTri);
        cPar->triangles(newTri).patch = 0;
        
    }
    
    // compute projection
    NormalProjector<ctype> projector(cPar);
    
    projector.project(surf2, domainDirection, targetDirection);

}


// ////////////////////////////////////////////////////////
//   Explicit template instantiations.
//   If you need more, you can add them here.
// ////////////////////////////////////////////////////////

template class ContactToolBox<float>;
template class ContactToolBox<double>;

#include <psurface/PSurface.h>
#include <psurface/contact.h>
#include <psurface/ContactToolBox.h>
#include <psurface/IntersectionPrimitive.h>

#include <vector>
#include <tr1/array>

static Parametrization* cPar;
static Surface* surf1;
static Surface* surf2;

void buildContactMapping(const std::vector<std::tr1::array<double,3> >& coords1,  ///< The vertices of the first surface as \f$x_0 ,y_0 ,z_0, x_1, y_1, z_1 ...\f$
                         const std::vector<std::tr1::array<int,3> >& tri1,       ///< The triangles of the first surface
                         const std::vector<std::tr1::array<double,3> >& coords2,  ///< The vertices of the second surface
                         const std::vector<std::tr1::array<int,3> >& tri2,
                         float epsilon, void (*obsDirections)(const double* pos, double* dir))
{
    int nVert1 = coords1.size();
    int nVert2 = coords2.size();
    int nTri1  = tri1.size();
    int nTri2  = tri2.size();

    // Create a first Surface object
    surf1 = new Surface;

#ifndef PSURFACE_STANDALONE
    // Amira Surface class needs a 'patch' structure
    surf1->patches.resize(1);
    surf1->patches[0] = new Surface::Patch;
    surf1->patches[0]->innerRegion = 0;
    surf1->patches[0]->outerRegion = 1;
    surf1->patches[0]->boundaryId  = 0;
#endif

    surf1->points.resize(nVert1);
    for (int i=0; i<nVert1; i++) 
        for (int j=0; j<3; j++)
            surf1->points[i][j] = coords1[i][j];

    surf1->triangles.resize(nTri1);
#ifndef PSURFACE_STANDALONE
    surf1->patches[0]->triangles.resize(nTri1);
#endif
    for (int i=0; i<nTri1; i++) {

        surf1->triangles[i].points[0] = tri1[i][0];
        surf1->triangles[i].points[1] = tri1[i][1];
        surf1->triangles[i].points[2] = tri1[i][2];

#ifndef PSURFACE_STANDALONE
        surf1->triangles[i].patch = 0;
        surf1->patches[0]->triangles[i] = i;
#endif
    }

    // Create a second Surface object
    surf2 = new Surface;
    
#ifndef PSURFACE_STANDALONE
    surf2->patches.resize(1);
    surf2->patches[0] = new Surface::Patch;
    surf2->patches[0]->innerRegion = 0;
    surf2->patches[0]->outerRegion = 1;
    surf2->patches[0]->boundaryId  = 0;
#endif

    surf2->points.resize(nVert2);
    for (int i=0; i<nVert2; i++)
        for (int j=0; j<3; j++)
            surf2->points[i][j] = coords2[i][j];

    surf2->triangles.resize(nTri2);
#ifndef PSURFACE_STANDALONE
    surf2->patches[0]->triangles.resize(nTri2);
#endif
    for (int i=0; i<nTri2; i++) {

        surf2->triangles[i].points[0] = tri2[i][0];
        surf2->triangles[i].points[1] = tri2[i][1];
        surf2->triangles[i].points[2] = tri2[i][2];

#ifndef PSURFACE_STANDALONE
        surf2->triangles[i].patch = 0;
        surf2->patches[0]->triangles[i] = i;
#endif
    }

#if defined HAVE_AMIRAMESH || !defined PSURFACE_STANDALONE
    // For debugging
    surf1->write("testSurf1.surf", 1);
    surf2->write("testSurf2.surf", 1);
#endif

    // Create the parametrization
    if (!cPar)
        delete cPar;

    cPar = new Parametrization;

    ContactToolBox::buildContactSurface(cPar, surf1, surf2, epsilon, obsDirections);

}


void getMergedGrid(std::vector<IntersectionPrimitive<2,float> >& overlaps)
{
    ContactToolBox::extractMergedGrid(cPar, overlaps);
}

void deleteContactSurface()
{
    if (surf1 != NULL) {
        delete surf1;
        surf1 = NULL;
    }

    if (surf2 != NULL) {
        delete surf2;
        surf2 = NULL;
    }

    if (cPar != NULL) {
        delete cPar;
        cPar = NULL;
    }
}



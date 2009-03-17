#include <vector>
#include <string>

#include <psurface/Parametrization.h>
#include <psurface/psurface.h>
#include <amiramesh/AmiraMesh.h>

// Using a map instead of a vector here seems nice, but I need
// the index operator once.
std::vector<std::pair<std::string, Parametrization*> > domains;

Parametrization* currentDomain;

int psurface::LoadMesh(const char* label, const char* filename)
{
    AmiraMesh* am = AmiraMesh::read(filename);
    
    if (!am) {
        printf("An error has occured while reading %s\n", filename);
        return ERROR;
    }

    Parametrization* newDomain = Parametrization::readAmiraMesh(am, filename);

    // The following delete should be here but it sometimes crashes...
    //delete(am);

    if (!newDomain) {
        printf("An error has occured while reading %s\n", filename);
        return ERROR;
    }
    
    domains.push_back(std::pair<std::string, Parametrization*>(label, newDomain));

    return OK;    
}


int psurface::RemoveDomain(const char* label)
{
    std::string domainLabel = label;

    std::vector<std::pair<std::string, Parametrization*> >::iterator i = domains.begin();

    for (; i!= domains.end(); ++i)
        if (i->first == domainLabel) {

            if (i->second == currentDomain)
                currentDomain = NULL;

            delete i->second;
            domains.erase(i);

            return OK;

        }

    return ERROR;
}

void psurface::RemoveAllDomains()
{
    std::vector<std::pair<std::string,Parametrization*> >::iterator i = domains.begin();
    for (; i!=domains.end(); ++i) 
        delete i->second;

    domains.clear();
    currentDomain = NULL;
}




int psurface::StartEditingDomain(const char* label)
{
    std::string domainLabel = label;

    std::vector<std::pair<std::string, Parametrization*> >::iterator i = domains.begin();

    for (; i!=domains.end(); ++i)
        if (i->first == domainLabel) {

            currentDomain = i->second;

            // create the point location data structure, if necessary
            if (!currentDomain->hasUpToDatePointLocationStructure)
                currentDomain->createPointLocationStructure();

            return OK;

        }

    return ERROR;
}

// This is the main access function
int psurface::CallParametrization(int tri, double* p, int* res, double* coords, int seed)
{
    McVec2f input(p[0], p[1]);
    McVec3i result;
    McVec2f localCoords;

    int status = currentDomain->map(tri, input, result, localCoords, seed);

    res[0] = result[0];
    res[1] = result[1];
    res[2] = result[2];

    coords[0] = localCoords[0];
    coords[1] = localCoords[1];
    coords[2] = localCoords[2];

    return (status) ? OK : ERROR;
}

// This is a shortcut function for the surface point position
int psurface::CallPositionParametrization(int tri, double* p, double* res)
{
    McVec2f input(p[0], p[1]);
    McVec3f result; 

    int status = currentDomain->positionMap(tri, input, result);

    res[0] = result[0];
    res[1] = result[1];
    res[2] = result[2];

    return (status) ? OK : ERROR;
}

// This is a shortcut function for the surface point position
// Here we can explicitly choose a parametrization from the list
int psurface::CallPositionParametrizationForDomain(int domain, int tri, double* p, double* res)
{
    McVec2f input(p[0], p[1]);
    McVec3f result; 

    int status = domains[domain].second->positionMap(tri, input, result);

    res[0] = result[0];
    res[1] = result[1];
    res[2] = result[2];

    return (status) ? OK : ERROR;
}

// This is the surface normals access function
int psurface::CallDirectNormalParametrization(int tri, double* p, double* res)
{
    McVec2f input(p[0], p[1]);
    McVec3f result; 

    int status = currentDomain->directNormalMap(tri, input, result);

    res[0] = result[0];
    res[1] = result[1];
    res[2] = result[2];

    return (status) ? OK : ERROR;
}


// Returns the number vertices in the currently chosen parametrization
int psurface::GetNoOfNodes()
{
    return currentDomain->getNumVertices();
}

// Returns the number triangles in the currently chosen parametrization
int psurface::GetNoOfSegments()
{
    return currentDomain->getNumTriangles();
}

int psurface::GetNoOfPatches()
{
    return currentDomain->patches.size();
}

int psurface::GetNoOfDomains()
{
    return domains.size();
}

void psurface::GetNodeNumbersOfSegment(int* points, int tri)
{
    points[0] = currentDomain->triangles(tri).vertices[0];
    points[1] = currentDomain->triangles(tri).vertices[1];
    points[2] = currentDomain->triangles(tri).vertices[2];
}


void psurface::GetLeftAndRightSideOfSegment(int* left, int* right, int tri)
{
    *left  = currentDomain->patches[currentDomain->triangles(tri).patch].innerRegion;
    *right = currentDomain->patches[currentDomain->triangles(tri).patch].outerRegion;
}

int psurface::GetBoundaryIdOfSegment(int tri)
{
    return currentDomain->patches[currentDomain->triangles(tri).patch].boundaryId;
}

int psurface::GetPatchNoOfSegment(int tri)
{
    return currentDomain->triangles(tri).patch;
}


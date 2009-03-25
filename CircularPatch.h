#ifndef CIRCULAR_PATCH
#define CIRCULAR_PATCH

#include <psurface/Box.h>
#include <mclib/McSArray.h>

#include "Domains.h"
#include "Parametrization.h"

#include "psurfaceAPI.h"


class Parametrization;

/** This class represents the retriangulation of a small hole in a surface.
    It's basically a bunch of triangles.  The feature of this class are the
    routines that evaluate the quality of the retriangulation. */
class PSURFACE_API CircularPatch {
public:
    

    /**@name constructors */
    //@{
    ///
    CircularPatch() {
        triangles.resize(0); 
        innerEdges.resize(0);
        par = NULL;
    }

    ///
    CircularPatch(Parametrization* param) {
        triangles.resize(0); 
        innerEdges.resize(0);
        par = param;
    }

    ///
    CircularPatch(int size, Parametrization* param) {
        triangles.resize(size);
        triangles.fill(-1U);

        innerEdges.resize(size-1);
        McSArray<int, 2> emptyArray;
        emptyArray.fill(-1U);
        innerEdges.fill(emptyArray);

        par = param;
    }

    ///
    CircularPatch(const std::vector<int>& array, Parametrization* param) {
        triangles.resize(array.size());
        for (size_t i=0; i<array.size(); i++)
            triangles[i] = array[i];

        par = param;
    }

    //@}

    ///
    void resize(int size) {
        triangles.resize(size);
        triangles.fill(-1);
        
        innerEdges.resize(size-1);
        McSArray<int, 2> emptyArray;
        emptyArray.fill(-1);
        innerEdges.fill(emptyArray);
    }

    ///
    int& last() {
        return triangles.last();
    }
    
    ///
    int& operator[](int i) {
        return triangles[i];
    }

    ///
    const int operator[](int i) const {
        return triangles[i];
    }
    
    ///
    int size() const { return triangles.size(); }

    ///
    void getBoundingBox(Box<float,3> &bbox) const;

    ///
    void killAll(){
        for (int i=0; i<triangles.size(); i++){
            par->removeTriangle(triangles[i]);
        }
    }

    /**@name Evaluatation methods */
    //@{
    /// gets the smallest internal angle found in any of the triangles
    float getMinInteriorAngle() const{
        int i;
        float minAngle = 2*M_PI;
        for (i=0; i<size(); i++){
            float currentMinAngle = par->minInteriorAngle(i);
            if (currentMinAngle<minAngle)
                minAngle = currentMinAngle;
        }

        return minAngle;
    }

    ///
    bool hasSmallDihedralAngles(float threshold, const Parametrization* par, 
                                 const DomainVertex* centerVertex) const;

    /// returns the largest triangle aspect ratio
    float maxAspectRatio() const {
        int i;
        float maxRatio = 0;
        for (i=0; i<size(); i++){
            const float currentAspectRatio = par->aspectRatio(i);
            if (currentAspectRatio>maxRatio)
                maxRatio = currentAspectRatio;
        }

        return maxRatio;
    }

    ///
    bool hasSelfintersections() const;

    /// tests whether the patch intersects the given parametrization, except for the immediate neighbors of center
    bool intersectsParametrization(const std::vector<int> &closeEdges) const;

    /// tests whether insertion of this patch would lead to topology changes
    bool inducesTopologyChange() const;

    //@}

    float distanceTo(const class McVec3f &) const ;

    McSmallArray<McSArray<int, 2>, 50> innerEdges;

private:

    McSmallArray<int, 50> triangles;
public:
    Parametrization* par;

};

#endif

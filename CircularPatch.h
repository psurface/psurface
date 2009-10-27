#ifndef CIRCULAR_PATCH
#define CIRCULAR_PATCH

#include <psurface/Box.h>
#include <psurface/Domains.h>
#include <psurface/PSurface.h>


template <int dim, class ctype>
class PSurface;

/** This class represents the retriangulation of a small hole in a surface.
    It's basically a bunch of triangles.  The feature of this class are the
    routines that evaluate the quality of the retriangulation. */
template <class ctype>
class CircularPatch {
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
    CircularPatch(PSurface<2,ctype>* param) {
        triangles.resize(0); 
        innerEdges.resize(0);
        par = param;
    }

    ///
    CircularPatch(int size, PSurface<2,ctype>* param) {
        triangles.resize(size);
        triangles.assign(size,-1U);

        innerEdges.resize(size-1);
        std::tr1::array<int, 2> emptyArray;
        emptyArray.assign(-1U);
        innerEdges.assign(innerEdges.size(),emptyArray);

        par = param;
    }

    ///
    CircularPatch(const std::vector<int>& array, PSurface<2,ctype>* param) {
        triangles.resize(array.size());
        for (size_t i=0; i<array.size(); i++)
            triangles[i] = array[i];

        par = param;
    }

    //@}

    ///
    void resize(int size) {
        triangles.resize(size);
        triangles.assign(size,-1);
        
        innerEdges.resize(size-1);
        std::tr1::array<int, 2> emptyArray;
        emptyArray.assign(-1);
        innerEdges.assign(innerEdges.size(), emptyArray);
    }

    ///
    int& last() {
        return triangles.back();
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
    void getBoundingBox(Box<ctype,3> &bbox) const;

    ///
    void killAll(){
        for (int i=0; i<triangles.size(); i++){
            par->removeTriangle(triangles[i]);
        }
    }

    /**@name Evaluatation methods */
    //@{
    /// gets the smallest internal angle found in any of the triangles
    ctype getMinInteriorAngle() const{
        int i;
        ctype minAngle = 2*M_PI;
        for (i=0; i<size(); i++){
            ctype currentMinAngle = par->minInteriorAngle(i);
            if (currentMinAngle<minAngle)
                minAngle = currentMinAngle;
        }

        return minAngle;
    }

    ///
    bool hasSmallDihedralAngles(ctype threshold, const PSurface<2,ctype>* par, 
                                 const McVertex<ctype>* centerVertex) const;

    /// returns the largest triangle aspect ratio
    ctype maxAspectRatio() const {
        int i;
        ctype maxRatio = 0;
        for (i=0; i<size(); i++){
            const ctype currentAspectRatio = par->aspectRatio(i);
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

    ctype distanceTo(const class StaticVector<ctype,3> &) const ;

    std::vector<std::tr1::array<int, 2> > innerEdges;

private:

    std::vector<int> triangles;
public:
    PSurface<2,ctype>* par;

};

#endif

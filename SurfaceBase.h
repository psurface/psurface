#ifndef SURFACE_BASE_H
#define SURFACE_BASE_H

#include <vector>

// Check for VC9 / VS2008 with installed feature pack.
#if defined(_MSC_VER) && (_MSC_VER>=1500)
    #if defined(_CPPLIB_VER) && _CPPLIB_VER>=505
        #include <array>
    #else
        #error Please install the Visual Studio 2008 SP1 for TR1 support.
    #endif
#else
    #include <tr1/array>
#endif

#include <limits>

#include "StaticVector.h"
#include "SurfaceParts.h"

#include "psurfaceAPI.h"

/** A simple pointer-based surface.  The @c VertexType, @c EdgeType, and @c TriangleType
    classes have to be derived from McVertex, McEdge, McTriangle,
    respectively.
    
    For an example of the use of this class look at the module HxParametrization,
    which is derived from McPointerSurface, and its three constituent classes
    DomainVertex, DomainEdge, and DomainTriangle.
    @see McVertex, McEdge, McTriangle
*/
template <class VertexType, class EdgeType, class TriangleType>
class PSURFACE_API SurfaceBase {

public:

    /** \brief The type used for coordinates */
    typedef typename VertexType::coordtype ctype;

    ///
    SurfaceBase(){}

    ///
    void clear() {
        triangleArray.resize(0);
        freeTriangleStack.resize(0);
        edgeArray.resize(0);
        freeEdgeStack.resize(0);
        vertexArray.resize(0);
        freeVertexStack.resize(0);
    }

    /**@name Procedural Access to the Elements */
    //@{
    ///
    TriangleType& triangles(int i) {return triangleArray[i];}

    ///
    const TriangleType& triangles(int i) const {return triangleArray[i];}

    ///
    EdgeType& edges(int i) {return edgeArray[i];}

    ///
    const EdgeType& edges(int i) const {return edgeArray[i];}

    ///
    VertexType& vertices(int i) {return vertexArray[i];}

    ///
    const VertexType& vertices(int i) const {return vertexArray[i];}

    ///
    size_t getNumTriangles() const {
        return triangleArray.size();
    }

    ///
    size_t getNumEdges() const {
        return edgeArray.size();
    }

    ///
    size_t getNumVertices() const {
        return vertexArray.size();
    }

    //@}

    /**@name Insertion and Removal of Elements */
    //@{
    /// removes a triangle and maintains a consistent data structure
    void removeTriangle(int tri);

    /// removes an edge
    void removeEdge(int edge){
        vertices(edges(edge).from).removeReferenceTo(edge);
        vertices(edges(edge).to).removeReferenceTo(edge);

        freeEdgeStack.push_back(edge);
    };

    void removeVertex(int vertex) {
        assert(!vertices(vertex).degree());
        freeVertexStack.push_back(vertex);
    }

    int newVertex(const StaticVector<ctype,3>& p);

    int newEdge(int a, int b);
    
    int createSpaceForTriangle(int a, int b, int c);

    void integrateTriangle(int triIdx);
    //@}

    /**@name Topological Queries */
    //@{

    /** Given two vertices, this routine returns the edge that connects them,
        or -1, if no such edge exists. */
    int findEdge(unsigned int a, unsigned int b) const;

    /** Given three vertices, this routine returns the triangle that connects them,
        or -1, if no such triangle exists. */
    int findTriangle(int a, int b, int c) const;

    /// Tests whether the two edges are connected by a common triangle
    int findCommonTriangle(int a, int b) const;

    /// 
    std::vector<int> getTrianglesPerVertex(int v) const;
    ///
    std::vector<int> getNeighbors(int v) const;

    int getNeighboringTriangle(int tri, int side) const;
    
    //@}

    /**@name Geometrical Queries */
    //@{

    /// gives the smallest interior angle of a triangle
    ctype minInteriorAngle(int n) const;

    /// returns the aspect ratio
    ctype aspectRatio(int n) const;

        /// returns the normal vector
    StaticVector<ctype,3> normal(int tri) const {
        const StaticVector<ctype,3> a = vertices(triangles(tri).vertices[1]) - vertices(triangles(tri).vertices[0]);
        const StaticVector<ctype,3> b = vertices(triangles(tri).vertices[2]) - vertices(triangles(tri).vertices[0]);
        StaticVector<ctype,3> n = a.cross(b);
        n.normalize();
        return n;
    }

    ///
    ctype smallestDihedralAngle(int edge) const {
        ctype minAngle = std::numeric_limits<ctype>::max();
        for (int i=0; i<edges(edge).triangles.size(); i++)
            for (int j=i+1; j<edges(edge).triangles.size(); j++)
                minAngle = std::min(minAngle,dihedralAngle(edges(edge).triangles[i], edges(edge).triangles[j]));

        return minAngle;
    }

    /// gives the surface area
    ctype area(int tri) const { 
        StaticVector<ctype,3> a = vertices(triangles(tri).vertices[1]) - vertices(triangles(tri).vertices[0]);
        StaticVector<ctype,3> b = vertices(triangles(tri).vertices[2]) - vertices(triangles(tri).vertices[0]);

        return fabs(0.5 * (a.cross(b)).length());
    }

    /// gives the dihedral angle with a neighboring triangle
    ctype dihedralAngle(int first, int second) const {
        StaticVector<ctype,3> n1 = normal(first);
        StaticVector<ctype,3> n2 = normal(second);

        ctype scalProd = n1.dot(n2);
        if (scalProd < -1) scalProd = -1;
        if (scalProd >  1) scalProd =  1;

        return (triangles(first).isCorrectlyOriented(triangles(second))) ? acos(-scalProd) : acos(scalProd);
    }

    ctype length(int e) const {
        return (vertices(edges(e).from) - vertices(edges(e).to)).length();
    }

    //@}

     /**@name Intersection Tests */
    //@{

    /** Tests whether this triangle intersects the given edge.  This routine is not
        faster than the one that returns the intersection point. */
    bool intersectionTriangleEdge(int tri, 
                                  const McEdge*edge,
                                  ctype eps=0) const;

    /** Tests whether this triangle intersects the given edge, and returns the intersection
        point if there is one. If not, the variable @c where is untouched. */
    bool intersectionTriangleEdge(int tri, 
                                  const McEdge *edge, 
                                  StaticVector<ctype,3>& where, 
                                  bool& parallel, ctype eps=0) const;

    /// Tests whether the point is inside the triangle given by the three argument points.
    static bool pointInTriangle(const StaticVector<ctype,2>& p,
                                const StaticVector<ctype,2>& a, 
                                const StaticVector<ctype,2>& b, 
                                const StaticVector<ctype,2>& c, ctype eps=0);
				
    static bool lineIntersection2D(const StaticVector<ctype,2> &p1, const StaticVector<ctype,2> &p2, const StaticVector<ctype,2> &p3, 
				   const StaticVector<ctype,2> &p4, ctype eps=0);
    //@}

public:
    void garbageCollection();
    
protected:
    ///
    std::vector<TriangleType> triangleArray;

    ///
    std::vector<VertexType> vertexArray;

    ///
    std::vector<EdgeType> edgeArray;

public:
    ///
    std::vector<int> freeTriangleStack;
protected:
    ///
    std::vector<int>     freeEdgeStack;

    ///
    std::vector<int>  freeVertexStack;

};

#endif

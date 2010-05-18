#ifndef SURFACE_BASE_H
#define SURFACE_BASE_H

#include <vector>
#include <set>
#include <tr1/array>
#include <algorithm>
#include <limits>
#include <iostream>

#include <psurface/StaticVector.h>
#include <psurface/StaticMatrix.h>
#include "McPointerSurfaceParts.h"


/** A simple pointer-based surface.  The @c VertexType, @c EdgeType, and @c TriangleType
    classes have to be derived from McVertex, McEdge, McTriangle,
    respectively.
    
    For an example of the use of this class look at the module HxParametrization,
    which is derived from McPointerSurface, and its three constituent classes
    DomainVertex, DomainEdge, and DomainTriangle.
    @see McVertex, McEdge, McTriangle
*/
template <class VertexType, class EdgeType, class TriangleType>
class SurfaceBase {

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
    void removeTriangle(int tri){
        int i;
        
        // update all three neighboring edges
        for (i=0; i<3; i++){
            
            int thisEdge = triangles(tri).edges[i];
            if (thisEdge==-1)
                continue;
            
            if (edges(thisEdge).triangles.size() == 1){
                // remove this edge
                removeEdge(thisEdge);
            }else {
                edges(thisEdge).removeReferenceTo(tri);           
            }
            
            triangles(tri).edges[i] = -1;
        }

        // remove triangle
        //triangles(tri).invalidate();
        freeTriangleStack.push_back(tri);
    }

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

    void addVertex(VertexType* vertex) {
        vertices.append(vertex);
    }


    int newVertex(const StaticVector<ctype,3>& p) {

        if (freeVertexStack.size()){
            int newVertexIdx = freeVertexStack.back();
            freeVertexStack.pop_back();
            vertices(newVertexIdx) = p;
            return newVertexIdx;
        } 

        vertexArray.push_back(p);
        return vertexArray.size()-1;
    }

    int newEdge(int a, int b) {

        int newEdgeIdx;
        if (freeEdgeStack.size()) {
            newEdgeIdx = freeEdgeStack.back();
            freeEdgeStack.pop_back();
        } else {
            edgeArray.push_back(EdgeType());
            newEdgeIdx = edgeArray.size()-1;
        }

        EdgeType& newEdge = edges(newEdgeIdx);

        newEdge.from = a;
        newEdge.to   = b;

        newEdge.triangles.resize(0);
        
        return newEdgeIdx;
    }

    int createSpaceForTriangle(int a, int b, int c) {
        int newTri;
        if (freeTriangleStack.size()) {
            newTri = freeTriangleStack.back();
            freeTriangleStack.pop_back();
            triangleArray[newTri] = TriangleType(a,b,c);
        } else {
            triangleArray.push_back(TriangleType(a,b,c));
            newTri = triangleArray.size()-1;
        }

        return newTri;
    }

    void integrateTriangle(int triIdx) {

        TriangleType& tri = triangles(triIdx);
        
        // look for and possible create the new edges
        for (int i=0; i<3; i++){
            
            int pointA = tri.vertices[i];
            int pointB = tri.vertices[(i+1)%3];
            
            int thisEdge = findEdge(pointA, pointB);
            
            if (thisEdge == -1){
                
                // this edge doesn't exist yet
                int newEdgeIdx = newEdge(pointA, pointB);
                
                vertices(pointA).edges.push_back(newEdgeIdx);
                vertices(pointB).edges.push_back(newEdgeIdx);
                
                edges(newEdgeIdx).triangles.push_back(triIdx);
                tri.edges[i] = newEdgeIdx;
            }
            else{
                if (!edges(thisEdge).isConnectedToTriangle(triIdx))
                    edges(thisEdge).triangles.push_back(triIdx);

                tri.edges[i] = thisEdge;
            }
        }
    }
    //@}

    /**@name Topological Queries */
    //@{

    /** Given two vertices, this routine returns the edge that connects them,
        or -1, if no such edge exists. */
    int findEdge(unsigned int a, unsigned int b) const {
        assert(a>=0 && a<getNumVertices() && b>=0 && b<getNumVertices());
        for (int i=0; i<vertices(a).degree(); i++)
            if (edges(vertices(a).edges[i]).from == b ||
                edges(vertices(a).edges[i]).to   == b)
                return vertices(a).edges[i];
        
        return -1;
    }

    /** Given three vertices, this routine returns the triangle that connects them,
        or -1, if no such triangle exists. */
    int findTriangle(int a, int b, int c) const {
        assert(a>=0 && a<getNumVertices());
        assert(b>=0 && b<getNumVertices());
        assert(c>=0 && c<getNumVertices());

        int oneEdge = findEdge(a, b);
        
        if (oneEdge==-1)
            return -1;
        for (int i=0; i<edges(oneEdge).triangles.size(); i++)
            if (triangles(edges(oneEdge).triangles[i]).isConnectedTo(c)){
                assert(edges(oneEdge).triangles[i]<getNumTriangles());
                return edges(oneEdge).triangles[i];
            }
        return -1;
    }

    /// Tests whether the two edges are connected by a common triangle
    int findCommonTriangle(int a, int b) const {
        assert(a!=b);
        assert(a>=0 && a<getNumVertices());
        assert(b>=0 && b<getNumVertices());

        for (int i=0; i<edges(a).triangles.size(); i++)
            for (int j=0; j<edges(b).triangles.size(); j++)
                if (edges(a).triangles[i] == edges(b).triangles[j])
                    return edges(a).triangles[i];
        
        return -1;
    }

    /// 
    std::vector<int> getTrianglesPerVertex(int v) const {

        const VertexType& cV = vertices(v);

        // A temporary set for fast searching and insertion
        std::set<int> resultSet;

        for (size_t i=0; i<cV.edges.size(); i++) {

            const EdgeType& cE = edges(cV.edges[i]);
            resultSet.insert(cE.triangles.begin(),cE.triangles.end());

        }

        // copy set to std::vector;
        std::vector<int> result(resultSet.begin(), resultSet.end());

        return result;
    }

    ///
    std::vector<int> getNeighbors(int v) const {

        const VertexType& cV = vertices(v);
        std::vector<int> result;
        
        for (int i=0; i<cV.edges.size(); i++) {
            
            const EdgeType& cE = edges(cV.edges[i]);
            result.push_back(cE.theOtherVertex(v));

        }

        return result;
    }

    int getNeighboringTriangle(int tri, int side) const {
        assert(side>=0 && side<3);
        int neighboringTri = -1;
        int cE = triangles(tri).edges[side];
                
        if (edges(cE).triangles.size()==2) {
            neighboringTri = (edges(cE).triangles[0]==tri) 
                ? edges(cE).triangles[1]
                : edges(cE).triangles[0];
        }
        return neighboringTri;
    }
    
    //@}

    /**@name Geometrical Queries */
    //@{

    /// gives the smallest interior angle of a triangle
    ctype minInteriorAngle(int n) const {
        ctype minAngle = 2*M_PI;
        const std::tr1::array<int, 3>& p = triangles(n).vertices;

        for (int i=0; i<3; i++){
            StaticVector<ctype,3> a = vertices(p[(i+1)%3]) - vertices(p[i]);
            StaticVector<ctype,3> b = vertices(p[(i+2)%3]) - vertices(p[i]);

            ctype angle = acosf(a.dot(b) / (a.length() * b.length()));
            if (angle<minAngle)
                minAngle = angle;
        }

        return minAngle;
    }

    /// returns the aspect ratio
    ctype aspectRatio(int n) const {

        const std::tr1::array<int, 3>& p = triangles(n).vertices;

        const ctype a = (vertices(p[1]) - vertices(p[0])).length();
        const ctype b = (vertices(p[2]) - vertices(p[1])).length();
        const ctype c = (vertices(p[0]) - vertices(p[2])).length();

        const ctype aR = 2*a*b*c/((-a+b+c)*(a-b+c)*(a+b-c));

        // might be negative due to inexact arithmetic
        return fabs(aR);
    }

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
                                  ctype eps=0) const {
        bool parallel;
        StaticVector<ctype,3> where;
        return intersectionTriangleEdge(tri, edge, where, parallel, eps);
    }

    /** Tests whether this triangle intersects the given edge, and returns the intersection
        point if there is one. If not, the variable @c where is untouched. */
    bool intersectionTriangleEdge(int tri, 
                                  const McEdge *edge, 
                                  StaticVector<ctype,3>& where, 
                                  bool& parallel, ctype eps=0) const{

        const TriangleType& cT = triangles(tri);

        const StaticVector<ctype,3> &p = vertices(edge->from);
        const StaticVector<ctype,3> &q = vertices(edge->to);
        const StaticVector<ctype,3> &a = vertices(cT.vertices[0]);
        const StaticVector<ctype,3> &b = vertices(cT.vertices[1]);
        const StaticVector<ctype,3> &c = vertices(cT.vertices[2]);

        // Cramer's rule
        ctype det = StaticMatrix<ctype,3>(b-a, c-a, p-q).det();
        if (det<-eps || det>eps){
            
            // triangle and edge are not parallel
            parallel = false;

            ctype nu = StaticMatrix<ctype,3>(b-a, c-a, p-a).det() / det;
            if (nu<-eps || nu>1+eps) return false;

            ctype lambda = StaticMatrix<ctype,3>(p-a, c-a, p-q).det() / det;
            if (lambda<-eps) return false;

            ctype mu = StaticMatrix<ctype,3>(b-a, p-a, p-q).det() / det;
            if (mu<-eps) return false;

            if (lambda+mu > 1+eps) 
                return false;
            else {
                where = p + nu*(q-p);
                return true;
            }

        } else {

            // triangle and edge are parallel
            parallel = true;
            
            ctype alpha = StaticMatrix<ctype,3>(b-a, c-a, p-a).det();
            if (alpha<-eps || alpha>eps)
                return false;
            else {

                int i;

                // 2D intersection test

                // project onto the coordinate plane that is 'most parallel' to the triangle
                StaticVector<ctype,3> normal = (b-a).cross(c-a);

                StaticVector<ctype,2> a2D, b2D, c2D, p2D, q2D;

                for (i=0; i<3; i++) 
                    if (normal[i]<0)
                        normal[i] = -normal[i];

                for (i=0; i<3; i++)
                    if (normal[i]>=normal[(i+1)%3] && normal[i]>=normal[(i+2)%3]) {

                        a2D = StaticVector<ctype,2>(a[(i+1)%3], a[(i+2)%3]);
                        b2D = StaticVector<ctype,2>(b[(i+1)%3], b[(i+2)%3]);
                        c2D = StaticVector<ctype,2>(c[(i+1)%3], c[(i+2)%3]);
                        p2D = StaticVector<ctype,2>(p[(i+1)%3], p[(i+2)%3]);
                        q2D = StaticVector<ctype,2>(q[(i+1)%3], q[(i+2)%3]);

                    }

                if (pointInTriangle(p2D, a2D, b2D, c2D, eps) ||
                    pointInTriangle(q2D, a2D, b2D, c2D, eps))
                    return true;
                else
                    return (lineIntersection2D(p2D, q2D, a2D, b2D, eps) ||
                            lineIntersection2D(p2D, q2D, b2D, c2D, eps) ||
                            lineIntersection2D(p2D, q2D, a2D, c2D, eps));

            }
                
        }

    }
        
protected:

    /// Tests whether the point is inside the triangle given by the three argument points.
    static bool pointInTriangle(const StaticVector<ctype,2>& p,
                                const StaticVector<ctype,2>& a, 
                                const StaticVector<ctype,2>& b, 
                                const StaticVector<ctype,2>& c, ctype eps=0) {
         StaticVector<ctype,3> localBarycentricCoords;

        // McMat3f(this->x, b.x, c.x,  this->y, b[1], c[1],  1, 1, 1).det();
        ctype area0 = p[0] * (b[1]-c[1]) - b[0] * (p[1] - c[1]) + c[0] * (p[1] - b[1]);

        // McMat3f(a[0], this->x, c[0],  a[1], this->y, c[1],  1, 1, 1).det();
        ctype area1 = a[0] * (p[1]-c[1]) - p[0] * (a[1] - c[1]) + c[0] * (a[1] - p[1]);

        // McMat3f(a[0], b[0], c[0],  a[1], b[1], c[1],  1, 1, 1).det();
        ctype areaTotal = a[0] * (b[1]-c[1]) - b[0] * (a[1] - c[1]) + c[0] * (a[1] - b[1]);

        localBarycentricCoords[0] = area0/areaTotal;
        localBarycentricCoords[1] = area1/areaTotal;
        localBarycentricCoords[2] = 1 - localBarycentricCoords[0] - localBarycentricCoords[1];

        return (localBarycentricCoords[0]>=-eps && localBarycentricCoords[1]>=-eps && localBarycentricCoords[2]>=-eps);
    }

    static bool lineIntersection2D(const StaticVector<ctype,2> &p1, const StaticVector<ctype,2> &p2, const StaticVector<ctype,2> &p3, const StaticVector<ctype,2> &p4, ctype eps=0) {
        const StaticVector<ctype,2> A = p2 - p1;
        const StaticVector<ctype,2> B = p3 - p4;
        const StaticVector<ctype,2> C = p1 - p3;
        
        ctype det = A[1]*B[0] - A[0]*B[1];

        // 1D intersection
        if (det>=-eps && det<=eps)
            return (  ((p3-p1).length() + (p3-p2).length()) / (p1-p2).length() < 1+eps ||
                      ((p4-p1).length() + (p4-p2).length()) / (p1-p2).length() < 1+eps ||
                      ((p2-p3).length() + (p2-p4).length()) / (p3-p4).length() < 1+eps ||
                      ((p1-p3).length() + (p1-p4).length()) / (p3-p4).length() < 1+eps   );
        
        ctype mu     = (A[0]*C[1] - A[1]*C[0]) / det;
        ctype lambda = (B[1]*C[0] - B[0]*C[1]) / det;

        return (mu>-eps && mu<1+eps && lambda>-eps && lambda<1+eps);
    }

    //@}

public:
    /// \todo The copying could be sped up considerably...
    void garbageCollection()
    {
        int i, j;
#ifndef NDEBUG
        std::cout << "This is the SurfaceBase garbage collection..." << std::endl;
        std::cout << freeVertexStack.size() << " vertices, "
                  << freeEdgeStack.size()   << " edges, "
                  << freeTriangleStack.size() << " triangles removed" << std::endl;
#endif

        std::vector<bool> isInvalid;

        // clean up vertices
        if (freeVertexStack.size()) {
            
            int offset = 0;

            std::vector<int> vertexOffsets(vertexArray.size());
            isInvalid.resize(vertexArray.size());

            for (i=0; i<isInvalid.size(); i++)
                isInvalid[i] = false;

            for (i=0; i<freeVertexStack.size(); i++)
                isInvalid[freeVertexStack[i]] = true;

            for (i=0; i<vertexArray.size(); i++){
                vertexOffsets[i] = offset;
                
                if (isInvalid[i]) 
                    offset++;
            }

            ////////////////////
            for (i=0; i<vertexOffsets.size(); i++)
                vertexArray[i-vertexOffsets[i]] = vertexArray[i];
            
            vertexArray.resize(vertexArray.size()-offset);
            
            
            // Adjust edges
            for (i=0; i<edgeArray.size(); i++) {
                edgeArray[i].from -= vertexOffsets[edgeArray[i].from];
                edgeArray[i].to   -= vertexOffsets[edgeArray[i].to];
            }                

            // Adjust triangles
            for (i=0; i<triangleArray.size(); i++) 
                for (j=0; j<3; j++) 
                    triangleArray[i].vertices[j] -= vertexOffsets[triangleArray[i].vertices[j]];

            freeVertexStack.resize(0);                
        }

        // clean up edges
        if (freeEdgeStack.size()) {
            
            int offset = 0;

            std::vector<int> edgeOffsets(edgeArray.size());
            isInvalid.resize(edgeArray.size());

            for (i=0; i<isInvalid.size(); i++)
                isInvalid[i] = false;
            
            for (i=0; i<freeEdgeStack.size(); i++)
                isInvalid[freeEdgeStack[i]] = true;
            
            for (i=0; i<edgeArray.size(); i++){
                edgeOffsets[i] = offset;
                
                if (isInvalid[i]) 
                    offset++;
            }
            //printf("1.3\n");
            ////////////////////
            for (i=0; i<edgeOffsets.size(); i++)
                edgeArray[i-edgeOffsets[i]] = edgeArray[i];
            
            edgeArray.resize(edgeArray.size()-offset);
            
            // Adjust vertices
            for (i=0; i<vertexArray.size(); i++) 
                for (j=0; j<vertexArray[i].degree(); j++)
                    vertexArray[i].edges[j] -= edgeOffsets[vertexArray[i].edges[j]];
            
            // Adjust triangles
            for (i=0; i<triangleArray.size(); i++) 
                for (j=0; j<3; j++) 
                    if (triangleArray[i].edges[j]>=0)
                        triangleArray[i].edges[j] -= edgeOffsets[triangleArray[i].edges[j]];

            freeEdgeStack.resize(0);
        }

        // clean up triangles
        if (freeTriangleStack.size()) {
            
            int offset = 0;

            std::vector<int> triangleOffsets(triangleArray.size());
            isInvalid.resize(triangleArray.size());

            for (i=0; i<isInvalid.size(); i++)
                isInvalid[i] = false;

            for (i=0; i<freeTriangleStack.size(); i++)
                isInvalid[freeTriangleStack[i]] = true;

            for (i=0; i<triangleArray.size(); i++){
                triangleOffsets[i] = offset;
                
                if (isInvalid[i]) 
                    offset++;
            }

            ////////////////////
            for (i=0; i<triangleOffsets.size(); i++)
                triangleArray[i-triangleOffsets[i]] = triangleArray[i];
            
            triangleArray.resize(triangleArray.size()-offset);
            
            
            // Adjust edges
            for (i=0; i<edgeArray.size(); i++) 
                for (j=0; j<edgeArray[i].triangles.size(); j++) 
                    edgeArray[i].triangles[j] -= triangleOffsets[edgeArray[i].triangles[j]];
                
            freeTriangleStack.resize(0);
        }


#ifndef NDEBUG
        printf("   ...Garbage collection finished!\n");
#endif
    }
    
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

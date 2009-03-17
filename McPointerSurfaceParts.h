#ifndef MC_POINTER_SURFACE_PARTS
#define MC_POINTER_SURFACE_PARTS

#include <mclib/McVec2f.h>
#include <mclib/McVec3f.h>
#include <psurface/Box.h>
#include <mclib/McMat3f.h>
#include <vector>
#include <mclib/McSmallArray.h>
#include <mclib/McSArray.h>


class McTriangle;
class McVertex;
template<class VertexType>  class McEdge;

///////////////////////////////////////////////////////////////
/** This is the base class for vertices in a McSurfaceBase.  The template argument
    @c VertexType that you instantiate McSurfaceBase with has to be derived from it.  
    You cannot use it directly because this class has itself as a template parameter.  (You 
    would get 
    <tt>McVertex &lt McVertex &lt McVertex &lt ... &gt &gt &gt, YourEdgeType, YourTriangleType &gt </tt>)
    However, any derived class will do.
    @see McPointerSurface, McEdge, McTriangle
*/
///////////////////////////////////////////////////////////////
class MCLIB_API McVertex: public McVec3f
{
public:
    McVertex() {}
    McVertex(const McVec3f &a) : McVec3f(a) {}

    ~McVertex() {}

    std::vector<int> edges;

    /// the number of edges connected to this vertex
    int degree() const {return edges.size();}

    ///
    void removeReferenceTo(int edge){
        /** \todo Can be implemented more elegantly using the find() method */
        for (int i=0; i<edges.size(); i++)
            if (edges[i] == edge){
                //edges.remove(i);
                edges.erase(edges.begin() + i);
                break;
            }
    }

    ///
    void replaceReferenceTo(int a, int b){
        for (int i=0; i<edges.size(); i++)
            if (edges[i] == a){
                edges[i] = b;
                break;
            }
    }
};

///////////////////////////////////////////////////////////////
/** The base class for an edge in a McPointerSurface 
    @see McPointerSurface, McVertex, McTriangle
*/
///////////////////////////////////////////////////////////////
template<class VertexType>
class MCLIB_API McEdge
{
public:
    ///
    McEdge(){}

    ///
    ~McEdge(){}
    
    ///
    McEdge(int a, int b)
        {
            to = a;
            from = b;
        }

    ///
    int numTriangles() const {return triangles.size();}

    ///
    int theOtherVertex(int point) const {
        return (point == to) ? from : to;
    }

    ///
    bool isConnectedTo(int vertex) const {
        return to==vertex || from==vertex;
    }

    ///
    bool isConnectedToTriangle(int tri) const {
        for (int i=0; i<triangles.size(); i++)
            if (triangles[i]==tri)
                return true;

        return false;
    }

    ///
    bool removeReferenceTo(int tri){
        for(int i=0; i<triangles.size(); i++)
            if (triangles[i] == tri){
                triangles.remove(i);
                return true;
            }

        return false;
    }

    ///
    bool replaceReferenceTo(int a, int b){
        for(int i=0; i<triangles.size(); i++)
            if (triangles[i] == a){
                triangles[i] = b;
                return true;
            }

        return false;
    }

    /// The connection to a is replaced by the connection to b.
    bool replaceReferenceToVertex(int a, int b){
        if (from==a){
            from = b;
            return true;
        } else if (to==a){
            to = b;
            return true;
        } else
            return false;
    }


#if 0
    /** tests whether this edge intersects a given triangle.  No code here,
        the equivalent function in McTriangle is called. */
    bool intersects(const TriangleType* triangle, float eps=0) const {
        return triangle->intersects(this, eps);
    }

    /// tests for intersection with a box, for use with a McOctree
    bool intersect(const McBox3f &box, void* userData=0) const {
        VertexType* vertices = (VertexType*)userData;
        if (box.contains(vertices[this->from]) || box.contains(vertices[this->to]))
            return true;

        return (intersectsXYPatch( McBox2f(box[0], box[1], box[2], box[3]), box[4], vertices) ||
                intersectsXYPatch( McBox2f(box[0], box[1], box[2], box[3]), box[5], vertices) ||
                intersectsXZPatch( McBox2f(box[0], box[1], box[4], box[5]), box[2], vertices) ||
                intersectsXZPatch( McBox2f(box[0], box[1], box[4], box[5]), box[3], vertices) ||
                intersectsYZPatch( McBox2f(box[2], box[3], box[4], box[5]), box[0], vertices) ||
                intersectsYZPatch( McBox2f(box[2], box[3], box[4], box[5]), box[1], vertices));
     }

    bool intersectsXYPatch(const McBox2f& rect, float z, const VertexType* vertices) const {

        const McVec3f& f = vertices[from];
        const McVec3f& t = vertices[to];

        if ( (f.z < z && t.z < z) || (f.z > z && t.z > z))
            return false;

        float lambda = (z - f.z) / (t.z - f.z);
        
        McVec2f intersection = McVec2f(f.x + lambda*(t.x - f.x), f.y + lambda*(t.y - f.y));

        return rect.intersect(intersection);
    }

    bool intersectsXZPatch(const McBox2f& rect, float y, const VertexType* vertices) const {

        const McVec3f& f = vertices[from];
        const McVec3f& t = vertices[to];

        if ( (f.y < y && t.y <y) || (f.y > y && t.y > y))
            return false;

        float lambda = (y - f.y) / (t.y - f.y);
        
        McVec2f intersection = McVec2f(f.x + lambda*(t.x - f.x), f.z + lambda*(t.z - f.z));

        return rect.intersect(intersection);
    }

    bool intersectsYZPatch(const McBox2f& rect, float x, const VertexType* vertices) const {

        const McVec3f& f = vertices[from];
        const McVec3f& t = vertices[to];

        if ( (f.x < x && t.x <x) || (f.x > x && t.x > x))
            return false;

        float lambda = (x - f.x) / (t.x - f.x);
        
        McVec2f intersection = McVec2f(f.y + lambda*(t.y - f.y), f.z + lambda*(t.z - f.z));

        return rect.intersect(intersection);
    }
#endif    

    ///////////////////////////////////////////
    ///
    int from, to;

    ///
    McSmallArray<int, 2> triangles;

};


///////////////////////////////////////////////////////////////
/** The base class for triangles in a McPointerSurface
    @see McPointerSurface, McVertex, McEdge
 */
///////////////////////////////////////////////////////////////

class MCLIB_API McTriangle 
{
public:
    /// default constructor
    McTriangle() {
        vertices[0] = vertices[1] = vertices[2] = -1;
        edges[0]  = edges[1]  = edges[2]  = -1;
    }

    /// creates a triangle using three given vertices
    McTriangle(int vertexIdx[3]) {

        vertices[0] = vertexIdx[0];
        vertices[1] = vertexIdx[1];
        vertices[2] = vertexIdx[2];

        edges[0] = edges[1] = edges[2] = -1;
    }

    ///
    McTriangle(int a, int b, int c) {
        vertices[0] = a;
        vertices[1] = b;
        vertices[2] = c;

        edges[0] = edges[1] = edges[2] = -1;
    }

    ~McTriangle() {}

    ///
    bool isConnectedTo(int vertex) const {
        return vertices[0]==vertex || vertices[1]==vertex || vertices[2]==vertex;
    }

    ///
    bool isCorrectlyOriented(const McTriangle& other) const {
        for (int i=0; i<3; i++)
            for (int j=0; j<3; j++) 
                if (vertices[i]==other.vertices[j] && vertices[(i+1)%3]==other.vertices[(j+1)%3])
                    return false;
                else if (vertices[i]==other.vertices[(j+1)%3] && vertices[(i+1)%3]==other.vertices[j])
                    return true;

        assert(false);
        return false;
    }

    ///
    int getThirdVertex(int a, int b) const {
        for (int i=0; i<3; i++)
            if (vertices[i]!=a && vertices[i]!=b)
                return vertices[i];

        return -1;
    }

    ///
    int getThirdVertex(int e) const {
        for (int i=0; i<3; i++)
            if (edges[i]==e)
                return vertices[(i+2)%3];

        return -1;
    }

    /** \brief Determines which of the three corners of the triangle
     * is a given vertex v.
     *
     * \return 0, 1, or 2 if everything is okay. -1 if the triangle 
     * doesn't touch v.
     */
    int getCorner(int v) const {
        for (int i=0; i<3; i++)
            if (vertices[i]==v)
                return i;
        return -1;
    }

    /** \brief Determines which of the three edges of the triangle
     * is a given edge v.
     *
     * \return 0, 1, or 2 if everything is okay. -1 if the triangle 
     * doesn't touch v.
     */
    int getEdge(int v) const {
        for (int i=0; i<3; i++)
            if (edges[i]==v)
                return i;
        return -1;
    }

    ///
    int getOppositeEdge(int a) const {
        assert(vertices[0]==a || vertices[1]==a || vertices[2]==a);
        for (int i=0; i<3; i++)
            if (vertices[i]==a)
                return edges[(i+1)%3];

        return -1;
    }

    int getCommonEdgeWith(const McTriangle& other) const {
        int i, j;
        for (i=0; i<3; i++)
            for (j=0; j<3; j++)
                if (edges[i]==other.edges[j])
                    return edges[i];

        return -1;
    }

    ///
    bool replaceReferenceToEdge(int a, int b) {
        for (int i=0; i<3; i++)
            if (edges[i]==a){
                edges[i] = b;
                return true;
            }

        return false;
    }

    ///
    bool replaceReferenceToVertex(int a, int b) {
        for (int i=0; i<3; i++)
            if (vertices[i]==a){
                vertices[i] = b;
                return true;
            }

        return false;
    }


    //////////////////////////////////////////////////////////////
public:

    McSArray<int, 3> vertices;

    McSArray<int, 3> edges;

};



#endif

#ifndef MC_POINTER_SURFACE_PARTS
#define MC_POINTER_SURFACE_PARTS

#include <vector>

#include <psurface/StaticVector.h>
#include <psurface/StaticMatrix.h>
#include <psurface/Box.h>


/** This is the base class for vertices in a McSurfaceBase.  

\tparam ctype The type used for coordinates

    @see McPointerSurface, McEdge, McTriangle
*/
template <class ctype>
class McVertex: public StaticVector<ctype,3>
{
public:

    /** \brief Export coordinate type */
    typedef ctype coordtype;

    /** \brief Default constructor */
    McVertex() {}

    /** \brief Constructor for a given position */
    McVertex(const StaticVector<ctype,3> &a) 
        : StaticVector<ctype,3>(a) 
    {}

    std::vector<int> edges;

    /// the number of edges connected to this vertex
    int degree() const {return edges.size();}

    ///
    void removeReferenceTo(int edge){
        typename std::vector<int>::iterator it = std::find(edges.begin(), edges.end(), edge);
        if (it != edges.end())
            edges.erase(it);
    }

    ///
    void replaceReferenceTo(int a, int b){
        typename std::vector<int>::iterator it = std::find(edges.begin(), edges.end(), a);
        if (it != edges.end())
            *it = b;
    }
};

///////////////////////////////////////////////////////////////
/** The base class for an edge in a McPointerSurface 
    @see McPointerSurface, McVertex, McTriangle
*/
///////////////////////////////////////////////////////////////

class McEdge
{
public:
    ///
    McEdge(){}

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
        for (size_t i=0; i<triangles.size(); i++)
            if (triangles[i]==tri)
                return true;

        return false;
    }

    ///
    bool removeReferenceTo(int tri){
        std::vector<int>::iterator it = std::find(triangles.begin(), triangles.end(), tri);
        if (it != triangles.end()) {
            triangles.erase(it);
            return true;
        }

        return false;
    }

    ///
    bool replaceReferenceTo(int a, int b){
        for(size_t i=0; i<triangles.size(); i++)
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

    ///////////////////////////////////////////
    ///
    int from, to;

    ///
    std::vector<int> triangles;

};


///////////////////////////////////////////////////////////////
/** The base class for triangles in a McPointerSurface
    @see McPointerSurface, McVertex, McEdge
 */
///////////////////////////////////////////////////////////////

class McTriangle 
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

    std::tr1::array<int, 3> vertices;

    std::tr1::array<int, 3> edges;

};



#endif

#ifndef VTKIO_H
#define VTKIO_H

#include <vector>
#include <string>

namespace psurface{

// Forward declaration
namespace VTK {
  class VTUWriter;
}

template<class ctype,int dim>
class VTKIO{
    private:
    /// Psurface object to be read to vtu file.
    PSurface<dim,ctype>* par;

    /// Global coordinate of nodes
    std::vector<StaticVector<ctype,3> > nodePositions;
    /// Parameter Edge endpoints in global index
    std::vector<StaticVector<int,2> > parameterEdgeArray;

    ///Node type
    std::vector<typename Node<ctype>::NodeType>  nodeType;

    /// Number of triangle vertices
    int numVertices;
    /// Number of base grid triangles
    int numTriangles;
    /// Number of nodes
    int numNodes;
    /// Number of ParamEdge
    int numParamEdges;
    /// Number of total Edge Points
    int numEdgePoints;

    /// Total number of cells
    int ncells;
    /// Total number of points
    int nvertices;

    public:
    VTKIO(PSurface<dim,ctype>* psurface);

    ///write the parametrization into vtu file
    void createVTU(std::string filename, bool basegrid);

    ///write data file to stream
    void writeDataFile(std::ostream& s, bool basegrid);

    ///write the data section in vtu
    void writeAllData(VTK::VTUWriter& writer, bool basegrid);

    /// write point data
    void writePointData(VTK::VTUWriter& writer, bool basegrid);

    /// write the positions of vertices
    void writeGridPoints(VTK::VTUWriter& writer, bool basegrid);

    /// write the connectivity array
    void writeGridCells(VTK::VTUWriter& writer, bool basegrid);
};
}// namespace psurface
#endif

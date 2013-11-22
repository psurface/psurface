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

    public:
    VTKIO(PSurface<dim,ctype>* psurface);

    ///write the parametrization into vtu file
    void createVTU(const std::string& element_filename, const std::string& graph_filename);

private:
    ///write data file to stream
    void writeElementDataFile(std::ostream& s);

    /// write the positions of vertices
    void writeElementGridPoints(VTK::VTUWriter& writer);

    /// write the connectivity array
    void writeElementGridCells(VTK::VTUWriter& writer);

    /// write cell data
    void writeElementGridCellData(VTK::VTUWriter& writer);


    ///write data file to stream
    void writeGraphDataFile(std::ostream& s);

    /// write point data
    void writeGraphNodeTypes(VTK::VTUWriter& writer);

    /// write the positions of vertices
    void writeGraphGridPoints(VTK::VTUWriter& writer);

    /// write the connectivity array
    void writeGraphGridCells(VTK::VTUWriter& writer);
};
}// namespace psurface
#endif

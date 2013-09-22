#ifndef VTKIO_H
#define VTKIO_H
#include <vector>
//#include "psurface_convert_new.h"
#include "psurfaceAPI.h"
#include "vtuwriter.hh"
namespace psurface{
template<class ctype,int dim>
class VTKIO{
    public:
    PSurface<dim,ctype>* par;
    Surface* surf;
    std::vector<StaticVector<ctype,3> > nodePositions;
    std::vector<StaticVector<int,2> > parameterEdgeArray;
    std::vector<int>  numNodesAndEdgesArray;
    std::vector<int>  nodeType;

    int numVertices;
    int numTriangles;
    int numNodes;
    int numParamEdges;
    int numEdgePoints;
    
    int ncells;
    int nvertices;
    const VTK::OutputType outputtype = VTK::ascii;

    public:
    VTKIO(PSurface<dim,ctype>* psurface);

    ///write the psurface into vtu file
    void creatVTU(const char *filename, bool basegrid);

   ///write data file to stream
    void writeDataFile(std::ostream& s, bool basegrid);
 
    ///write the data section in vtu
    void writeAllData(VTK::VTUWriter& writer, bool basegrid); 

    /// write point data
    void writePointData(VTK::VTUWriter& writer, bool basegrid);

    //! write the positions of vertices
    void writeGridPoints(VTK::VTUWriter& writer, bool basegrid);

    //! write the connectivity array
    void writeGridCells(VTK::VTUWriter& writer, bool basegrid);
};
}
#endif

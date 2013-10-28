#ifndef HDF5IO_H
#define HDF5IO_H
#include <stdio.h>
#include <stdlib.h>

namespace psurface{
template<class ctype,int dim>
class Hdf5IO{
    private:
    /// Psurface object hdf5IO deals.
    PSurface<dim,ctype>* par;
    /// Global coordinate of nodes
    std::vector<StaticVector<ctype,3> > nodePositions;
    /// Parameter Edge endpoints in global index
    std::vector<StaticVector<int,2> > parameterEdgeArray;

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

    /// Writes the parametrization in hdf5 format with all data we need to read it in paraview.
    void writeHdf5Data(const std::string& filename);
    ///  Writes the parametrization in hdf5 format
    void writeBaseHdf5Data(const std::string& filename);

    public:
    Hdf5IO(PSurface<dim,ctype>* psurface);

    /** \brief Writes the parametrization in hdf5 format and create related xdmf file
     *
     *  \param base if true, we create hdf5 file that only contain necessary data
     *  of the psurface. If false, we create hdf5 file that contains all data
     *  needed to display the psurface object in paraview.
     */
    void createHdfAndXdmf(const std::string&  xdf_filename, const std::string&  hdf_filename, bool base);

    /// Writes xdmf file so we can read hdf5 throught paraview
    void writeXdmf(const std::string&  xdf_filename, const std::string&  hdf_filename);

    /// Reads the parametrization from Hdf5 object
    void initCompletePSurface(Surface* surf, const std::string&  filename);
};
}// namespace psurface
#endif

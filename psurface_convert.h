#ifndef PSURFACE_CONVERT_H
#define PSURFACE_CONVERT_H
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <memory>
#include <tr1/memory>
#include "StaticVector.h"
#include "Domains.h"
#include "PSurface.h"
#include "PSurfaceFactory.h"
#include "vtuwriter.hh"

//#ifdef PSURFACE_STANDALONE
#include "TargetSurface.h"
//#else
//#include <hxsurface/Surface.h>
//#endif
#if defined HAVE_AMIRAMESH || !defined PSURFACE_STANDALONE
#include <amiramesh/AmiraMesh.h>
#endif
namespace psurface{

template<class ctype,int dim>
class PsurfaceConvert{
  typedef typename PSurface<2,ctype>::Patch PATCH;
  typedef typename std::vector<StaticVector<ctype,3> >::iterator v_iterator;
  typedef typename std::vector<StaticVector<int,3> >::iterator t_iterator;
  typedef typename std::vector<StaticVector<int,2> >::iterator e_iterator;

  enum NodeTypes{INTERIOR_NODE,
  INTERSECTION_NODE,
  CORNER_NODE,
  TOUCHING_NODE,
  GHOST_NODE};
  public:
  ///The psurface object we need to read
  PSurface<dim,ctype>* par;
  /**@name data structure to store psurface nodes, triangles and edges on origrinal surface*/
  //@{
  ///the position of corner nodes in global coordinate
  std::vector<StaticVector<ctype,3> > baseGridVertexCoordsArray;
  ///the position of other nodes in global coordinate
  //namely, interior_node, intersection_node, touching_node. we do not consider ghost node
  std::vector<StaticVector<ctype,3> > nodePositions;
  ///the position of nodes in local coordinate on each triangle
  std::vector<StaticVector<ctype,2> > domainPositions;

  ///triangles
  std::vector<StaticVector<int, 3> > baseGridTriArray;

  ///edges
  std::vector<StaticVector<int,2> > parameterEdgeArray;
  std::vector<StaticVector<int,2> > parameterEdgeArrayLocal;

  //edgepoints
  std::vector<int>  edgePointsArray;

  ///nodetype
  std::vector<int> nodeType;

  ///nodeNumber
  std::vector<int> nodeNumber;

  ///patches
  std::vector<PATCH> patches;


  ///triangle index;
  std::vector<int> triId;
  ///coordinates of the image position
  std::vector<StaticVector<ctype,3> > imagePos;
  std::vector<StaticVector<ctype,3> > iPos;
  ///local information on planer graphs
  std::vector<int> numNodesAndEdgesArray;
  //@}

  int numVertices;
  int numTriangles;
  int numNodes;
  int numParamEdges;
  int numEdgePoints;
  ///total number of all type of nodes
  int ncells;
  ///total number of triangles and edges
  int nvertices;
  const VTK::OutputType outputtype = VTK::ascii;

  public:
  /**create hdf5 file and xdmf file. hdf5 file store the basic data in psurface. xdmf is used to read the hdf5 file by 
  *paraview.
  */
  bool creatHdfAndXdmf(const std::string&  xdf_filename, const std::string&  hdf_filename, bool base);

  bool writeBaseHdf5Data(const std::string& filename);

  ///write the data array into hdf5 data structure
  bool writeHdf5Data(const std::string& filename);

  ///writhe the xdmf file which store the structure information of hdf5 file.
  bool writeXdmf(const std::string&  xdf_filename, const std::string&  hdf_filename);
  
    ///write the psurface into vtu file
  bool creatVTU(const char *filename, bool basegrid);

   ///write data file to stream
  void writeDataFile(std::ostream& s, bool basegrid);
 
    ///write the data section in vtu
  void writeAllData(VTK::VTUWriter& writer, bool basegrid); 

  //! write point data
  virtual void writePointData(VTK::VTUWriter& writer, bool basegrid);

  //! write the positions of vertices
  void writeGridPoints(VTK::VTUWriter& writer, bool basegrid);

  //! write the connectivity array
  virtual void writeGridCells(VTK::VTUWriter& writer, bool basegrid);
  
  ///initialize PsurfaceConvert from the psurface object
  PsurfaceConvert(PSurface<dim,ctype>* psurface);

  ///initialize PsurfaceConvert from file
  ///read data array from hdf5 or gmsh file into the data arrary in class PsurfaceConvert
  PsurfaceConvert(const std::string&  filename, bool a);

  bool readHdf5Data(const std::string&  filename);
 
  ///read psurface_convert from Gmsh file
  bool readGmsh(const std::string&  filename);

  ///creat surface from psurface_convert
  ///baseTriangle == 1  we create the psurface subject which only have base grid triangles
  bool initPsurface(PSurface<2,ctype>* psurf, Surface* surf, bool baseTriangle);
 
  ///Create PSurface which only have BaseGridVertex and BaseTriangle
  bool initBaseCasePSurface(PSurface<2,ctype>* psurf, Surface* surf);

  bool initCompletePSurface(PSurface<2,ctype>* psurf, Surface* surf);

};
}
#endif

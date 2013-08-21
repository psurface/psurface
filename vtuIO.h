//the code to write psurface object into vtk file
//vtuwrite
#ifndef VTUIO_H
#define VTUIO_H
#include <vector>
#include "StaticVector.h"
#include "Domains.h"
#include "PSurface.h"
#include "PSurfaceFactory.h"
#include "dune/vtuwriter.hh"
#include "dune/shared_ptr.hh"
#include <fstream>
using namespace psurface;
enum NodeTypes{INTERIOR_NODE,
INTERSECTION_NODE,
CORNER_NODE,
TOUCHING_NODE,
GHOST_NODE};

Dune::VTK::OutputType outputtype = Dune::VTK::ascii;
template<class ctype,int dim>
class VtkPWriter{
/**
 *@brief class to write psurface subject into vtu format.
 *@tParams ctype the type used for the coordinates of psurface.
 *@tParams dim  the dimension of the psurface.
 */
  public:
  ///The psurface object we need to read
  PSurface<dim,ctype>* par; //initialize from the input
  /**@name data structure to store psurface nodes, triangles and edges*/
  //@{
  ///the position of corner nodes in global coordinate
  std::vector<StaticVector<ctype,3> > baseGridVertexCoordsArray;
  ///the position of other nodes in global coordinate
  std::vector<StaticVector<ctype,3> > domainPositions;
  ///triangles
  std::vector<StaticVector<int, 3> > baseGridTriArray;
  ///edges
  std::vector<StaticVector<int,2> >parameterEdgeArray;
  ///patches
  std::vector<int> patches;
  ///nodetype
  std::vector<int> nodeType;
  ///coordinates of the image position
  std::vector<StaticVector<ctype,3> > imagePos;
  ///number of local nodes on each triangle(do not include corner nodes)
  std::vector<int> nodeperTriangle;
  ///edge number on each triangle(do not include triangle edges)
  std::vector<int> edgeperTriangle;
  //@}

  int numVertices;
  int numTriangles;
  int numNodes;
  int numParamEdges;
  ///total number of all type of nodes
  int ncells;
  ///total number of triangles and edges
  int nvertices;
  //output type(set as ascii)
  const Dune::VTK::OutputType outputtype = Dune::VTK::ascii;
  /**@brief compute the nodes index of inner nodes/intersection nodes/ghost node/touching node
   * @param a trianlge index
   * @param b nodes index on the triangle
   * @return global index of the nodes
   */
  int indexMapBack(int a,int b)
  {
    if( b > 3)
    {
      std::vector<int>::iterator pE;
      int index = numVertices;
      int i;
      for(i = 0; i < a; i++)
      {
        for(pE = nodeperTriangle.begin(); pE != nodeperTriangle.end(); pE++)
        index += (*pE);
      }
      index += b - 3;
      return index;
    }
    else
    return par->triangles(a).vertices[b];
  };
  int indexMapBack(int a)
  {
     return a;
  }
  void computeCellsAndVertices()
  {
    nvertices = numVertices + numNodes;
    ncells = numTriangles + numParamEdges;
  }
  ///write the psurface into unstructured grid
  void writeUnstructureMesh(const char *filename){
  std::ofstream file;
  file.open(filename);
  if (! file.is_open()) printf("error\n");
  writeDataFile( file );
  file.close();
  return;
  }
  ///write data file to stream
  void writeDataFile(std::ostream& s)
  {
    Dune::VTK::FileType fileType = Dune::VTK::unstructuredGrid;

    Dune::VTK::VTUWriter writer(s, outputtype,fileType);//Most inportant structure used here

    writer.beginMain(numTriangles + numParamEdges, numVertices + numNodes);
    writeAllData(writer);
    writer.endMain();
  }
  ///write the data section in vtu
  void writeAllData(Dune::VTK::VTUWriter& writer) {
  //PointData
  writePointData(writer);
  // CellData
  writeCellData(writer);
  // Points
  writeGridPoints(writer);
  // Cells
  writeGridCells(writer);
  }
  typedef typename std::vector<StaticVector<ctype,3> >::iterator v_iterator;
  typedef typename std::vector<StaticVector<int,3> >::iterator t_iterator;
  typedef typename std::vector<StaticVector<int,2> >::iterator e_iterator;
  //! write point data
  virtual void writePointData(Dune::VTK::VTUWriter& writer)
  {
    std::string scalars = "nodetype";
    std::string vectors = "imageposition";
    std::vector<int>::iterator pt;
    writer.beginPointData(scalars, vectors);
    {
      Dune::shared_ptr<Dune::VTK::DataArrayWriter<ctype> > p
      (writer.makeArrayWriter<ctype>(scalars, 1, nvertices));
      for (pt = nodeType.begin(); pt!= nodeType.end(); ++pt)
      p->write(*pt);
    }
    v_iterator pi;
    {
      Dune::shared_ptr<Dune::VTK::DataArrayWriter<ctype> > p
      (writer.makeArrayWriter<ctype>(vectors, 3, nvertices));
      for (pi = imagePos.begin(); pi!= imagePos.end(); ++pi)
      {
        for(int l = 0; l < 3; l++) p->write((*pi)[l]);
      }
    }
    writer.endPointData();
  }
  //! write cell data
  virtual void writeCellData(Dune::VTK::VTUWriter& writer)
  {
    //std::vector<StaticVector<int,3> >::iterator pp;
    std::vector<int>::iterator pp;
    std::string scalars = "patches";
    std::string vectors = "";
    writer.beginCellData(scalars, vectors);
    {
        Dune::shared_ptr<Dune::VTK::DataArrayWriter<ctype> > p
      (writer.makeArrayWriter<ctype>("patches", 1, ncells));
      for (pp = patches.begin(); pp!= patches.end(); ++pp)
      p->write(*pp);
    }
    writer.endCellData();
  }
  //! write the positions of vertices
  void writeGridPoints(Dune::VTK::VTUWriter& writer)
  {
    writer.beginPoints();
    v_iterator vp;
    v_iterator dp;
    {
      Dune::shared_ptr<Dune::VTK::DataArrayWriter<ctype> > p
      (writer.makeArrayWriter<ctype>("Coordinates", 3, nvertices));
      if(!p->writeIsNoop()) {
      for(vp =  baseGridVertexCoordsArray.begin(); vp != baseGridVertexCoordsArray.end(); vp++)
      for(int l = 0; l < 3; l++) p->write((*vp)[l]);
      for(dp = domainPositions.begin(); dp < domainPositions.end();dp++)
      for(int l = 0; l < 3; l++) p->write((*dp)[l]);
      }
    }
    //	p.reset();
    writer.endPoints();
  }
  //! write the connectivity array
  virtual void writeGridCells(Dune::VTK::VTUWriter& writer)
  {
    t_iterator tp;
    e_iterator ep;
    writer.beginCells();
    // connectivity
    {
      Dune::shared_ptr<Dune::VTK::DataArrayWriter<int> > p1
      (writer.makeArrayWriter<int>("connectivity", 1, 3*numTriangles +2*numParamEdges));
      if(!p1->writeIsNoop())
      {
        for( tp = baseGridTriArray.begin(); tp!=baseGridTriArray.end(); ++tp)
        for( int l = 0; l < 3; l++)
        p1->write((*tp)[l]);
        for( ep = parameterEdgeArray.begin(); ep != parameterEdgeArray.end(); ep++)
        {
          for(int l = 0; l < 2; l++)
          p1->write((*ep)[l]);
        }
      }
    }
    // offsets
    {
      Dune::shared_ptr<Dune::VTK::DataArrayWriter<int> > p2
      (writer.makeArrayWriter<int>("offsets", 1, ncells));
      if(!p2->writeIsNoop()) {
      int offset = 0;
      for( tp = baseGridTriArray.begin(); tp!=baseGridTriArray.end(); ++tp)
      {
        offset += 3;
        p2->write(offset);
      }
      for( ep = parameterEdgeArray.begin(); ep != parameterEdgeArray.end(); ep++)
      {
        offset += 2;
        p2->write(offset);
      }
      }
    }
    // types
    {
      int offset = 0;
      {
      Dune::shared_ptr<Dune::VTK::DataArrayWriter<unsigned char> > p3
      (writer.makeArrayWriter<unsigned char>("types", 1, ncells));
      if(!p3->writeIsNoop())
      {
        for( tp = baseGridTriArray.begin(); tp!=baseGridTriArray.end(); ++tp)
        p3->write(5); //vtktype of triangle
      }
      for( ep = parameterEdgeArray.begin(); ep != parameterEdgeArray.end(); ep++)
      {
        offset += 2;
        p3->write(3);//vtktype of edges
      }
      }
    }
    writer.endCells();
  }
  ///Read information from the psurface object and stored it in local vector containers in vtkPWriter
  void update(PSurface<2,ctype>* par0)
  {
    par = par0;
    //interate over elementu
    int i,j,k;
    numVertices  = par->getNumVertices();
    numTriangles = par->getNumTriangles();
    //base grid vertices
    for (size_t i=0; i< par->getNumVertices(); i++)
    {
      baseGridVertexCoordsArray.push_back(par->vertices(i));
      nodeType.push_back(CORNER_NODE);
    }
    std::vector<StaticVector<ctype,3> > imagearray;
    imagearray.resize(numVertices);
    StaticVector<int,3> vertices;
    StaticVector<ctype,3> imagepos;
    //base grid triangles
    for (size_t i=0; i<par->getNumTriangles(); i++)
    {
      for (int j=0; j<3; j++)
      {
        vertices[j] = par->triangles(i).vertices[j];
        imagepos = par->imagePos(i,j);
        imagearray[vertices[j]] = imagepos;
      }
      baseGridTriArray.push_back(vertices);
    }
    //imageine position of each vertices
    for (size_t i=0; i< par->getNumVertices(); i++)
    imagePos.push_back(imagearray[i]);
    //patchs on each Triangles
    //	StaticVector<int,3> patch;
    for (size_t i=0; i<par->getNumTriangles(); i++)
    {
      const DomainTriangle<ctype>& cT0 = par->triangles(i);
      patches.push_back(cT0.patch);
//		 patch[0] = cT0.patch.innerRegion;
//		 patch[1] = cT0.patch.outerRegion;
//		 patch[2] = ct0.patch.boundaryId;
//		 Patches.push_back(patch);
    }
    //plane graph on each base grid triangle, saved as a list of nodes and a list of edges.
    int numNode = 0; //total number of intersection + touching + interior nodes + corner node
    int numEdge = 0; // num of all edges of the graphs on all triagules
    for(i = 0; i < numTriangles;i++)
    {
      const DomainTriangle<ctype>& cT = par->triangles(i);
      int numIntersectionNodes;
      int numTouchingNodes;
      int numInteriorNodes;
      cT.countNodes(numIntersectionNodes, numTouchingNodes, numInteriorNodes);
      int nEdges = cT.getNumRegularEdges();
      numNode += numIntersectionNodes;
      numNode += numTouchingNodes;
      numNode += numInteriorNodes;
      numEdge += nEdges - 3;
      nodeperTriangle.push_back(numIntersectionNodes + numTouchingNodes + numInteriorNodes);
      edgeperTriangle.push_back(nEdges - 3);
    }
    numNodes = numNode;
    numParamEdges = numEdge;
    computeCellsAndVertices();
    //patchs on each edges are defined to be negetive 1
    for (i =0; i< numEdge; i++)
    {
      patches.push_back(-1);
    }
    StaticVector<ctype,3> nodePosition;
    StaticVector<int,2> edgearray;
    int arrayIdx           = 0;
    int edgeArrayIdx       = 0;
    int edgePointsArrayIdx = 0;
    for (i=0; i<numTriangles; i++) {
    int localArrayIdx = 0;
    const DomainTriangle<ctype>& cT = par->triangles(i);
    std::vector<int> newIdx(cT.nodes.size());
    // the three remaining types are saved separately, in order to avoid
    // explicitly saving the type for each node.
    for (size_t cN=0; cN<cT.nodes.size(); cN++) {
    if (cT.nodes[cN].isINTERSECTION_NODE()||cT.nodes[cN].isTOUCHING_NODE()||cT.nodes[cN].isINTERIOR_NODE()){
    if(cT.nodes[cN].isINTERSECTION_NODE()) nodeType.push_back(INTERSECTION_NODE);
    else if(cT.nodes[cN].isTOUCHING_NODE())
    nodeType.push_back(TOUCHING_NODE);
    else nodeType.push_back(INTERIOR_NODE);
    //get world coordinate from barycentric coordinate
    for(j = 0; j < 3; j++)
    {
      int indexx = par->triangles(i).vertices[j];
      int k;
      ctype corner_coords[3];
      for(k = 0; k <3 ; k++) corner_coords[i] = par->vertices(indexx)[k];
      ctype domain_position[2];
      for(k = 0; k < 2; k++) domain_position[k] = cT.nodes[cN].domainPos()[k];
      nodePosition[j] = corner_coords[0]*domain_position[0];

      domainPositions.push_back(nodePosition);
      newIdx[cN] = localArrayIdx;
      arrayIdx++;
      localArrayIdx ++;
    }
    //get the image positon of each nodes
    imagepos = par->imagePos(i,cN);
    imagePos.push_back(imagepos);
    }
    }
  // the parameterEdges for this triangle
    bool triangle_edge = 0;
    typename PlaneParam<ctype>::UndirectedEdgeIterator cE;
    for (cE = cT.firstUndirectedEdge(); cE.isValid(); ++cE){
    if (cE.isRegularEdge()) {
    edgearray[0] = indexMapBack(i,cE.from());
    edgearray[1] = indexMapBack(i,cE.to());
    if((cE.from() == 0 || cE.from() == 1 || cE.from() == 2) && (cE.to() == 0 || cE.to() == 1 || cE.to() == 2))
    triangle_edge = 1;
    else
    {
      parameterEdgeArray.push_back(edgearray);
      edgeArrayIdx++;
    }
    }
    }
    }
  return;
  }
  };
/**
 *@brief class to write 1d psurface subject into vtu format(a special case of VtkPWriter.
 *@tParams ctype the type used for the coordinates of psurface.
 */
template <class ctype>
class VtkPWriter<ctype,1>
{
public:
  ///The psurface object we need to read
  PSurface<1,ctype>* par1d;
  /**@name nodes and edges on original surface*/
  //@{
  ///Corner Nodes
  std::vector<StaticVector<ctype,2> > cornerNodes;
  ///Inner Nodes0
  std::vector<StaticVector<ctype,2> > innerNodes;
  ///Edges
  std::vector<StaticVector<int,2> > edges;
  ///Normal direction
  std::vector<StaticVector<ctype,2> > domainNorm;
  //@}
  /**@name nodes and edges on image surface*/
  //@{
  ///Target Nodes
  std::vector<StaticVector<ctype, 2> > targetNodes;
  ///Target Edges
  std::vector<std::tr1::array<int, 2> > targetSegments;
  //@}
  int numNodes,numCornerNodes,numInnerNodes;
  int numEdges;
  int numTargetSegments, numTargetNodes;
  //output type(set as ascii)
    ///write the origin surface into unstructured grid
  void writeDomain(const char *filename){
  std::ofstream file;
  file.open(filename);
  if (! file.is_open()) printf("error\n");
  writeDomainData(file);
  file.close();
  return;
  }
  ///write the image surface into unstructured grid
  void writeTarget(const char *filename){
  std::ofstream file;
  file.open(filename);
  if (! file.is_open()) printf("error\n");
  writeTargetData(file);
  file.close();
  return;
  }
  ///write domain data into the file
  void writeDomainData(std::ostream& s)
  {
      Dune::VTK::FileType fileType = Dune::VTK::unstructuredGrid;

      Dune::VTK::VTUWriter writer(s, outputtype,fileType);//Most inportant structure used here

    writer.beginMain(numEdges, numCornerNodes+ numInnerNodes);
    writeAllData1d(0,writer);
    writer.endMain();
  }
  ///write image surface data into the file
  void writeTargetData(std::ostream& s)
  {
      Dune::VTK::FileType fileType = Dune::VTK::unstructuredGrid;

      Dune::VTK::VTUWriter writer(s, outputtype,fileType);//Most inportant structure used here

      writer.beginMain(numTargetSegments,numTargetNodes);
      writeAllData1d(1,writer);
      writer.endMain();
  }
  ///write the data section in vtu
  void writeAllData1d(int i,Dune::VTK::VTUWriter& writer) {
  //PointData
  writePointData1d(i,writer);
  //Points
  writeGridPoints1d(i,writer);
  //cells
  writeGridCells1d(i,writer);
  }
  //! write point data
  virtual void writePointData1d(int otype, Dune::VTK::VTUWriter& writer)
  {
    int i;
    if(otype == 0)
    {
      std::string scalars = "nodetype";
//      std::string vectors = "imageposition";
      std::string vectors = "";
      int i;
      writer.beginPointData(scalars, vectors);
      {
        Dune::shared_ptr<Dune::VTK::DataArrayWriter<ctype> > p
        (writer.makeArrayWriter<ctype>(scalars, 1, numNodes));
        for (i = 0; i < cornerNodes.size();i++)
        p->write(0);
        for (i = 0; i < innerNodes.size(); i++)
        p->write(1);
      }
      writer.endPointData();
    }
  }
  //! write the positions of vertices
  void writeGridPoints1d(int otype,Dune::VTK::VTUWriter& writer)
  {
    int i;
    if(otype == 0)
    {
    writer.beginPoints();
    {
      Dune::shared_ptr<Dune::VTK::DataArrayWriter<ctype> > p
      (writer.makeArrayWriter<ctype>("Coordinates", 3, numNodes));
      if(!p->writeIsNoop()) {
      for(i = 0; i < cornerNodes.size(); i++)
      {
        for(int l = 0; l < 2; l++) p->write(cornerNodes[l]);
        p->write(0);
      }
      for(i = 0; i < innerNodes.size(); i++)
      {
        for(int l = 0; l < 2; l++) p->write(innerNodes[l]);
        p->write(0);
      }
      }
    }
    //p.reset();
    writer.endPoints();
    }
    else
    {
    writer.beginPoints();
    {
      Dune::shared_ptr<Dune::VTK::DataArrayWriter<ctype> > p
      (writer.makeArrayWriter<ctype>("Coordinates", 3, numTargetNodes));
      if(!p->writeIsNoop()) {
      for(i = 0; i < targetNodes.size(); i++)
      {
        for(int l = 0; l < 2; l++) p->write(targetNodes[l]);
        p->write(0);
      }
      }
    }
    //	p.reset();
    writer.endPoints();
    }
  }
  //! write the connectivity array
  virtual void writeGridCells1d(int otype,Dune::VTK::VTUWriter& writer)
  {
    int i;
    writer.beginCells();
    if(otype == 0)
    {
    // connectivity
    {
      Dune::shared_ptr<Dune::VTK::DataArrayWriter<int> > p1
      (writer.makeArrayWriter<int>("connectivity", 1, 2*numEdges));
      if(!p1->writeIsNoop())
      {
        for(i = 0; i < edges.size();i++)
        {
          for(int l = 0; l < 2; l++) p1->write((edges[i])[l]);
        }
      }
    }
    // offsets
    {
      Dune::shared_ptr<Dune::VTK::DataArrayWriter<int> > p2
      (writer.makeArrayWriter<int>("offsets", 1, numEdges));
      if(!p2->writeIsNoop()) {
      int offset = 0;
      for(i = 0; i < edges.size();i++)
      {
        offset += 2;
        p2->write(offset);
      }
      }
    }
    // types
    {
      int offset = 0;
      {
      Dune::shared_ptr<Dune::VTK::DataArrayWriter<unsigned char> > p3
      (writer.makeArrayWriter<unsigned char>("types", 1, numEdges));
      for(i = 0; i < edges.size();i++)
      {
        p3->write(3);//vtktype of edges
      }
      }
    }
    }
    else
    {
    // connectivity
    {
      Dune::shared_ptr<Dune::VTK::DataArrayWriter<int> > p1
      (writer.makeArrayWriter<int>("connectivity", 1, 2*numTargetSegments));
      if(!p1->writeIsNoop())
      {
        std::tr1::array<int, 2>  segment;
        for(i = 0; i < targetSegments.size();i++)
        {
          segment = targetSegments[i];
          for(int l = 0; l < 2; l++) p1->write(segment[l]);
        }
      }
    }
    // offsets
    {
      Dune::shared_ptr<Dune::VTK::DataArrayWriter<int> > p2
      (writer.makeArrayWriter<int>("offsets", 1, numTargetSegments));
      if(!p2->writeIsNoop()) {
      int offset = 0;
      for(i = 0; i < targetSegments.size();i++)
      {
        offset += 2;
        p2->write(offset);
      }
      }
    }
    // types
    {
      int offset = 0;
      {
      Dune::shared_ptr<Dune::VTK::DataArrayWriter<unsigned char> > p3
      (writer.makeArrayWriter<unsigned char>("types", 1, numTargetSegments));
      for(i = 0; i < targetSegments.size();i++)
      {
        p3->write(3);//vtktype of edges
      }
      }
    }
    }
    writer.endCells();
  }
  ///Read information from the psurface object and stored it in local vector containers in vtkPWriter
 void update1d(PSurface<1,ctype>* par0)
 {
    par1d = par0;
    int i, j, k;
    numCornerNodes = par1d->domainVertices.size();
    numInnerNodes = 0;
    for(i = 0; i < par1d->domainSegments.size();i++)
    {
       numInnerNodes += par1d->domainSegments[i].nodes.size();
       StaticVector<ctype,2> nodeend;
       for(j = 0; j < par1d->domainSegments[i].nodes.size();j++)
       {
         int loc =  par1d->domainSegments[i].nodes.domainLocalPosition;
         StaticVector<ctype,2> lver;
         StaticVector<ctype,2> rver;
         StaticVector<ctype,2> innernode;
         lver =  par1d->domainVertices[par1d->domainSegments[i].nodes.rangeSegments[0]];
         rver =  par1d->domainVertices[par1d->domainSegments[i].nodes.rangeSegments[1]];
         nodeend[0] = loc*lver[0] + (1 - loc)*rver[0];
         nodeend[1] = loc*lver[1] + (1 - loc)*rver[1];
         innerNodes.push_back(nodeend);
       }
    }
    numNodes = numCornerNodes + numInnerNodes;
    for(i = 0; i < par1d->domainVertices.size();i++)
      cornerNodes.push_back(par1d->domainVertices[i]);
    StaticVector<int,2> edgeEnd;
    for(i = 0; i < par1d->domainSegments.size();i++)
    {
      edgeEnd[0] = par1d->DomainSegment[i].points[0];
      edgeEnd[1] = par1d->DomainSegment[i].points[1];
      edges.push_back(edgeEnd);
    }
    for(i = 0; i < par1d->targetVertices.size();i++)
      targetNodes.push_back(par1d->targetVertices[i]);
    StaticVector<int,2> vertexEdgeEnd;
    for(i = 0; i < par1d->targetSegments.size();i++)
    {
      vertexEdgeEnd[0] = par1d->targetSegments[i][0];
      vertexEdgeEnd[1] = par1d->targetSegments[i][1];
      targetSegments.push_back(vertexEdgeEnd);
    }
    numTargetSegments = par1d->targetSegments.size();
    targetNodes = par1d->targetVertices.size();
  }
};
#endif

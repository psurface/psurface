#ifndef PSURFACE_CONVERT_H
#define PSURFACE_CONVERT_H
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <hdf5.h>
#include <fstream>
#include <memory>
#include <tr1/memory>
#include <amiramesh/HxParamBundle.h>
#include <amiramesh/HxParamBase.h>
#include "StaticVector.h"
#include "Domains.h"
#include "PSurface.h"
#include "PSurfaceFactory.h"
#include "vtuwriter.hh"

#ifdef PSURFACE_STANDALONE
#include "TargetSurface.h"
#else
#include <hxsurface/Surface.h>
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


  public:
  bool creatHdfAndXdmf(const char* xdf_filename, const char* hdf_filename)
  {
     writeHdf5Data(hdf_filename);
     writeXdmf(xdf_filename, hdf_filename);
     return 0;
  }

  ///write the data array into hdf5 data structure
  bool writeHdf5Data(const char*filename)
  {
    printf("in writeHdf5Date function: hdf_file = %s\n", filename); 
    hid_t     file_id;
    file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    int i;
    hid_t     dataset_id, dataspace_id,datatype;
    hsize_t   dims[] ={numVertices};
    herr_t    status;

    // Create the coordinate data of all types of nodes.
    float *x = (float *) malloc(nvertices * sizeof(float));
    float *y = (float *) malloc(nvertices * sizeof(float));
    float *z = (float *) malloc(nvertices * sizeof(float));
    //vertex
    for(i = 0; i < numVertices; i++)
    {
      x[i] = (baseGridVertexCoordsArray[i])[0];
      y[i] = (baseGridVertexCoordsArray[i])[1];
      z[i] = (baseGridVertexCoordsArray[i])[2];
    }
    //nodes
    for(i = 0; i < numNodes; i++)
    {
      x[i + numVertices] = (nodePositions[i])[0];
      y[i + numVertices] = (nodePositions[i])[1];
      z[i + numVertices] = (nodePositions[i])[2];
    }
    // Write separate coordinate arrays for the x y and z coordinates.
    //coordinate data array of vertices
    dims[0] = nvertices;
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dataset_id = H5Dcreate2(file_id, "/X", H5T_NATIVE_FLOAT, dataspace_id,
         H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, x);
    status = H5Sclose(dataspace_id);
    status = H5Dclose(dataset_id);
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dataset_id = H5Dcreate2(file_id, "/Y", H5T_NATIVE_FLOAT, dataspace_id,
    H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, y);
    status = H5Sclose(dataspace_id);
    status = H5Dclose(dataset_id);
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dataset_id = H5Dcreate2(file_id, "/Z", H5T_NATIVE_FLOAT, dataspace_id,
    H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, z);
    status = H5Sclose(dataspace_id);
    status = H5Dclose(dataset_id);
    free(x);
    free(y);
    free(z);

    //Creat connection array of triangle and edges
    //triangles
    int  *topo_array = (int *) malloc(4*numTriangles*sizeof(int) + 4 * numParamEdges*sizeof(int));
    for(i = 0; i < numTriangles;i++)
    {
      topo_array[4*i] = 4; //topology number of triangle in xdmf
      topo_array[4*i+1] = (baseGridTriArray[i])[0];
      topo_array[4*i+2] = (baseGridTriArray[i])[1];
      topo_array[4*i+3] = (baseGridTriArray[i])[2];
    }
    //edges
    for(i = 0; i < numParamEdges; i++)
    {
       topo_array[4*numTriangles + 4*i] = 2;
       topo_array[4*numTriangles + 4*i + 1] = 2;
       topo_array[4*numTriangles + 4*i + 2] = (parameterEdgeArray[i])[0];
       topo_array[4*numTriangles + 4*i + 3] = (parameterEdgeArray[i])[1];
    }
    dims[0] = numTriangles*4 + numParamEdges*4 ;
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dataset_id = H5Dcreate(file_id, "/Topo", H5T_NATIVE_INT, dataspace_id,
    H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
    H5P_DEFAULT, topo_array);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);
    free(topo_array);

    ///local position on trianlge
    int *tri_n = (int *) malloc(numNodes * sizeof(int));
    float *dp = (float *) malloc(2*numNodes * sizeof(float));
    for(i = 0; i < numNodes; i++)
    {
        tri_n[i] = triId[i];
        dp[2*i] = (domainPositions[i])[0];
        dp[2*i+1] = (domainPositions[i])[1];
    }
    //triangle index
    dims[0] = numNodes;
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dataset_id = H5Dcreate(file_id, "/TriangleIndx", H5T_NATIVE_INT,
    dataspace_id, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
    H5P_DEFAULT, tri_n);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);
    //local position of nodes on triangle
    dims[0] = numNodes*2;
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dataset_id = H5Dcreate(file_id, "/LocalNodePosition", H5T_NATIVE_FLOAT, dataspace_id,
    H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
    H5P_DEFAULT, dp);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);
    free(tri_n);
    free(dp);

    //creat connection array of parameter edges
    int *parameter_edge_local_array = (int*) malloc(2*numParamEdges*sizeof(int));
    for(i = 0; i < numParamEdges; i++)
    {
       parameter_edge_local_array[2*i] =  (parameterEdgeArrayLocal[i])[0];
       parameter_edge_local_array[2*i+1] = (parameterEdgeArrayLocal[i])[1];
    }
    //parameter edge on local array
    dims[0] = numParamEdges*2;
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dataset_id = H5Dcreate(file_id, "/LocalParamEdge", H5T_NATIVE_INT, dataspace_id,
    H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
    H5P_DEFAULT, parameter_edge_local_array);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);
    free( parameter_edge_local_array);

    //PATCH
    int *patches_vec = (int *) malloc(patches.size()*3*sizeof(int));
    for(i = 0; i < patches.size();i++)
    {
        patches_vec[3*i] = (patches[i]).innerRegion;
        patches_vec[3*i + 1] = (patches[i]).outerRegion;
        patches_vec[3*i + 2] = (patches[i]).boundaryId;
    }
    //patches
    dims[0] = patches.size()*3;
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dataset_id = H5Dcreate(file_id, "/patches_vec", H5T_NATIVE_INT, dataspace_id,
    H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
    H5P_DEFAULT, patches_vec);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);
    free(patches_vec);

    // create the scalar data on triangle
    int *patch = (int *) malloc(ncells * sizeof(int));
    for (i = 0; i < numTriangles; i++)
      patch[i] = numNodesAndEdgesArray[11*i+4];
    for(i = 0; i < numParamEdges;i++)
        patch[i + numTriangles] = 0;

    //patches
    dims[0] = numTriangles + numParamEdges;
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dataset_id = H5Dcreate(file_id, "/patches", H5T_NATIVE_INT, dataspace_id,
    H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
    H5P_DEFAULT, patch);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);
    free(patch);

    int *edgepointsarray = (int*)malloc(edgePointsArray.size()*sizeof(int));
    for(i = 0; i < edgePointsArray.size();i++)
      edgepointsarray[i] = edgePointsArray[i];
    //edgepointsarray
    dims[0] = edgePointsArray.size();
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dataset_id = H5Dcreate(file_id, "/EdgePoints", H5T_NATIVE_INT, dataspace_id,
    H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
    H5P_DEFAULT,edgepointsarray);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);
    free(edgepointsarray);

    //create the scalar data on nodes
    //nodetype
    int *nodetype = (int *) malloc(nvertices * sizeof(int));
    for (i = 0; i < nvertices; i++)
    {
      nodetype[i] = nodeType[i];
    }

    // nodestype
    dims[0] = nvertices;
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dataset_id = H5Dcreate(file_id, "/NodeType", H5T_NATIVE_INT,
    dataspace_id, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
    H5P_DEFAULT, nodetype);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);
    free(nodetype);

    //node number
    int *nodenumber = (int*)malloc(numNodes*sizeof(int));
    for(i = 0; i < numNodes; i++)
        nodenumber[i] = nodeNumber[i];

    //nodenumber
    dims[0] = numNodes;
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dataset_id = H5Dcreate(file_id, "/NodeNumber", H5T_NATIVE_INT,
    dataspace_id, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
    H5P_DEFAULT, nodenumber);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);
    free(nodenumber);

    //image position
    float* imageX = (float *) malloc(nvertices * sizeof(float));
    float* imageY = (float *) malloc(nvertices * sizeof(float));
    float* imageZ = (float *) malloc(nvertices * sizeof(float));
    for (i = 0; i< nvertices; i++)
    {
      imageX[i] = (imagePos[i])[0];
      imageY[i] = (imagePos[i])[1];
      imageY[i] = (imagePos[i])[2];
    }
    //images
    //x
    dims[0] = nvertices;
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dataset_id = H5Dcreate(file_id, "/imageX", H5T_NATIVE_FLOAT,
    dataspace_id, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
    H5P_DEFAULT, imageX);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);
    //y
    dims[0] = nvertices;
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dataset_id = H5Dcreate(file_id, "/imageY", H5T_NATIVE_FLOAT,
    dataspace_id, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
    H5P_DEFAULT, imageY);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);
    //z
    dims[0] = nvertices;
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dataset_id = H5Dcreate(file_id, "/imageZ", H5T_NATIVE_FLOAT,
    dataspace_id, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
    H5P_DEFAULT, imageZ);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);
    free(imageX);
    free(imageY);
    free(imageZ);

    //ipos
    float *ipos = (float *) malloc(iPos.size()*3*sizeof(float));
    for(i = 0; i < iPos.size(); i++)
    {
        ipos[3*i] = (iPos[i])[0];
        ipos[3*i+1] = (iPos[i])[1];
        ipos[3*i+2] = (iPos[i])[2];
    }
    //iPos
    dims[0] = iPos.size()*3;
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dataset_id = H5Dcreate(file_id, "/iPos", H5T_NATIVE_FLOAT, dataspace_id,
    H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
    H5P_DEFAULT, ipos);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);
    free(ipos);

    //numNodesandEdges array
    int *num_Nodes_and_Edges_Array = (int *) malloc(11*numTriangles*sizeof(int));
    for(i = 0; i < 11*numTriangles;i++)
        num_Nodes_and_Edges_Array[i] = numNodesAndEdgesArray[i];
    //nodeNodesAndEdgesArray
    dims[0] = numTriangles*11;
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dataset_id = H5Dcreate(file_id, "/numNodesAndEdgesArray", H5T_NATIVE_INT,
    dataspace_id, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
    H5P_DEFAULT, num_Nodes_and_Edges_Array);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);
    free(num_Nodes_and_Edges_Array);

    //params
    int *psurfaceparams = (int *)malloc(4*sizeof(int));
    psurfaceparams[0] = numVertices;
    psurfaceparams[1] = numNodes;
    psurfaceparams[2] = numTriangles;
    psurfaceparams[3] = numParamEdges;
    dims[0] = 4;
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dataset_id = H5Dcreate(file_id, "/Params", H5T_NATIVE_INT,
    dataspace_id, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
    H5P_DEFAULT, psurfaceparams);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);
    free(psurfaceparams);

    //close the file
    status = H5Fclose(file_id);
    return 0;
  };

  ///writhe the xdmf file which store the structure information of hdf5 file.
  bool writeXdmf(const char* xdf_filename, const char* hdf_filename)
  {
    printf("in writeXdmd function: xdf_file = %s hdf_file = %s\n", xdf_filename,hdf_filename);
    FILE *xmf = 0;
    xmf = fopen(xdf_filename, "w");
    fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
    fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
    fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
    fprintf(xmf, "<Domain>\n");
    fprintf(xmf, "<Grid Name=\"psurface\" GridType=\"Uniform\">\n");
    fprintf(xmf, "<Topology TopologyType=\"Mixed\" NumberOfElements=\"%d\">\n", numTriangles + numParamEdges);
    fprintf(xmf, "<DataItem Dimensions = \"%d\" NumberType=\"int\" Format=\"HDF\">\n",  4*numTriangles + 4*numParamEdges);
    fprintf(xmf, "%s:/Topo\n", hdf_filename);
    fprintf(xmf, "</DataItem>\n");
    fprintf(xmf, "</Topology>\n");
    fprintf(xmf, "<Geometry GeometryType=\"X_Y_Z\">\n");
    fprintf(xmf, "<DataItem Dimensions=\"%d\" NumberType=\"Float\" Format=\"HDF\">\n", nvertices);
    fprintf(xmf, "%s:/X\n", hdf_filename);
    fprintf(xmf, "</DataItem>\n");
    fprintf(xmf, "<DataItem Dimensions=\"%d\" NumberType=\"Float\" Format=\"HDF\">\n", nvertices);
    fprintf(xmf, "%s:/Y\n", hdf_filename);
    fprintf(xmf, "</DataItem>\n");
    fprintf(xmf, "<DataItem Dimensions=\"%d\" NumberType=\"Float\" Format=\"HDF\">\n", nvertices);
    fprintf(xmf, "%s:/Z\n", hdf_filename);
    fprintf(xmf, "</DataItem>\n");
    fprintf(xmf, "</Geometry>\n");
    fprintf(xmf, "<Attribute Name=\"Patches\" AttributeType=\"Scalar\" Center=\"Cell\">\n");
    fprintf(xmf, "<DataItem Dimensions=\"%d\" NumberType=\"int\" Format=\"HDF\">\n", ncells);
    fprintf(xmf, "%s:/patches\n", hdf_filename);
    fprintf(xmf, "</DataItem>\n");
    fprintf(xmf, "</Attribute>\n");
    fprintf(xmf, "<Attribute Name=\"Nodetype\" AttributeType=\"Scalar\" Center=\"Node\">\n");
    fprintf(xmf, "<DataItem Dimensions=\"%d\" NumberType=\"int\" Format=\"HDF\">\n", nvertices);
    fprintf(xmf, "%s:/NodeType\n", hdf_filename);
    fprintf(xmf, "</DataItem>\n");
    fprintf(xmf, "</Attribute>\n");
    fprintf(xmf, "<Attribute Name=\"imageX\" AttributeType=\"Scalar\" Center=\"Node\">\n");
    fprintf(xmf, "<DataItem Dimensions=\"%d\" NumberType=\"float\" Format=\"HDF\">\n", nvertices);
    fprintf(xmf, "%s:/imageX\n", hdf_filename);
    fprintf(xmf, "</DataItem>\n");
    fprintf(xmf, "</Attribute>\n");
    fprintf(xmf, "<Attribute Name=\"imageY\" AttributeType=\"Scalar\" Center=\"Node\">\n");
    fprintf(xmf, "<DataItem Dimensions=\"%d\" NumberType=\"float\" Format=\"HDF\">\n", nvertices);
    fprintf(xmf, "%s:/imageY\n", hdf_filename);
    fprintf(xmf, "</DataItem>\n");
    fprintf(xmf, "</Attribute>\n");
    fprintf(xmf, "<Attribute Name=\"imageZ\" AttributeType=\"Scalar\" Center=\"Node\">\n");
    fprintf(xmf, "<DataItem Dimensions=\"%d\" NumberType=\"float\" Format=\"HDF\">\n", nvertices);
    fprintf(xmf, "%s:/imageZ\n", hdf_filename);
    fprintf(xmf, "</DataItem>\n");
    fprintf(xmf, "</Attribute>\n");
    fprintf(xmf, "</Grid>\n");
    fprintf(xmf, "</Domain>\n");
    fprintf(xmf, "</Xdmf>\n");
    fclose(xmf);
    return true;
  };
  
  ///write the psurface into vtu file
  bool creatVTU(const char *filename)
  {
    std::ofstream file;
    file.open(filename);
    if (! file.is_open()) printf("%c does not exits!\n", filename);
    writeDataFile( file );
    file.close();
    return 0;
  }
  ///write data file to stream
  void writeDataFile(std::ostream& s)
  {
    VTK::OutputType outputtype;
    outputtype = VTK::ascii;
    VTK::FileType fileType = VTK::unstructuredGrid;

    VTK::VTUWriter writer(s, outputtype,fileType);//Most inportant structure used here

    writer.beginMain(numTriangles + numParamEdges, numVertices + numNodes);
    writeAllData(writer);
    writer.endMain();
  }
  ///write the data section in vtu
  void writeAllData(VTK::VTUWriter& writer) {
    //PointData
    writePointData(writer);
    // Points
    writeGridPoints(writer);
    // Cells
    writeGridCells(writer);
  }

  
  //! write point data
  virtual void writePointData(VTK::VTUWriter& writer)
  {
    std::string scalars = "nodetype";
    std::string vectors = "imageposition";
    std::vector<int>::iterator pt;
    writer.beginPointData(scalars, vectors);
    {
      std::tr1::shared_ptr<VTK::DataArrayWriter<ctype> > p
      (writer.makeArrayWriter<ctype>(scalars, 1, nvertices));
      for( int i = 0; i < nodeType.size();i++)
      p->write(nodeType[i]);
    }
    v_iterator pi;
    {
      std::tr1::shared_ptr<VTK::DataArrayWriter<ctype> > p
      (writer.makeArrayWriter<ctype>(vectors, 3, nvertices));
      for(int i = 0; i < imagePos.size(); i++)
      {
        for(int l = 0; l < 3; l++) p->write((imagePos[i])[l]);
      }
    }
    writer.endPointData();
  }


  //! write the positions of vertices
  void writeGridPoints(VTK::VTUWriter& writer)
  {
    writer.beginPoints();
    v_iterator vp;
    v_iterator dp;
    {
      std::tr1::shared_ptr<VTK::DataArrayWriter<ctype> > p
      (writer.makeArrayWriter<ctype>("Coordinates", 3, nvertices));
      if(!p->writeIsNoop()) {
      for(int i = 0; i < baseGridVertexCoordsArray.size(); i++)
        for(int l = 0; l < 3; l++) p->write((baseGridVertexCoordsArray[i])[l]);
      for(int i = 0; i < domainPositions.size();i++)
      {
        for(int l = 0; l < 3; l++) p->write((nodePositions[i])[l]);
      }
      }
    }
    //	p.reset();
    writer.endPoints();
  }

  //! write the connectivity array
  virtual void writeGridCells(VTK::VTUWriter& writer)
  {
    t_iterator tp;
    e_iterator ep;
    writer.beginCells();
    // connectivity
    {
      std::tr1::shared_ptr<VTK::DataArrayWriter<int> > p1
      (writer.makeArrayWriter<int>("connectivity", 1, 3*numTriangles +2*numParamEdges));
      if(!p1->writeIsNoop())
      {
        for(int i = 0; i < baseGridTriArray.size(); i++)
        for( int l = 0; l < 3; l++)
        p1->write((baseGridTriArray[i])[l]);
        for( int i = 0; i <  parameterEdgeArray.size(); i++)
        {
          for(int l = 0; l < 2; l++)
          p1->write((parameterEdgeArray[i])[l]);
        }
      }
    }
    // offsets
    {
      std::tr1::shared_ptr<VTK::DataArrayWriter<int> > p2
      (writer.makeArrayWriter<int>("offsets", 1, ncells));
      if(!p2->writeIsNoop()) {
      int offset = 0;
      for(int i = 0; i < baseGridTriArray.size(); i++)
      {
        offset += 3;
        p2->write(offset);
      }
      for(int i = 0; i < parameterEdgeArray.size(); i++)
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
      std::tr1::shared_ptr<VTK::DataArrayWriter<unsigned char> > p3
      (writer.makeArrayWriter<unsigned char>("types", 1, ncells));
      if(!p3->writeIsNoop())
      {
        for(int i = 0; i < baseGridTriArray.size();i++)
        p3->write(5); //vtktype of triangle
      }
      for(int i = 0; i < parameterEdgeArray.size(); i++) 
      {
        offset += 2;
        p3->write(3);//vtktype of edges
      }
      }
    }
    writer.endCells();
  }

  ///initializers
  ///initialize PsurfaceConvert from the psurface object
  PsurfaceConvert(PSurface<dim,ctype>* psurface)
  {
    par = psurface;
    ///interate over element
    int i,j,k;
    numVertices  = par->getNumVertices();
    numTriangles = par->getNumTriangles();
    ///coordinate of vertices
    for (size_t i=0; i< par->getNumVertices(); i++)
    {
      baseGridVertexCoordsArray.push_back(par->vertices(i));
      nodeType.push_back(CORNER_NODE);
    }

    ///base grid triangles
    StaticVector<int,3> vertices;
    for (size_t i=0; i<par->getNumTriangles(); i++)
    {
      for(j=0; j<3; j++)
        vertices[j] = par->triangles(i).vertices[j];
      baseGridTriArray.push_back(vertices);
    }

    ///patch
    patches.resize(par->patches.size());
    for (size_t i=0; i<par->patches.size(); i++)  patches[i] = par->patches[i];
    
    ///imagePosition on vertices
    std::vector<StaticVector<float,3> > imagearray;
    imagearray.resize(numVertices);
    for (size_t i=0; i<par->getNumTriangles(); i++)
    {
      for (int j=0; j<3; j++)
      {
        vertices[j] = par->triangles(i).vertices[j];
        imagearray[vertices[j]] = par->imagePos(i,j);
      }
    }

    
    for (size_t i=0; i< par->getNumVertices(); i++)
    {
      imagePos.push_back(imagearray[i]);
    }

    ///nodes image positions(only for intersection points)
    for(i = 0; i < par->iPos.size();i++) iPos.push_back(par->iPos[i]);

    numNodes      = 0;
    numParamEdges = 0;
    numEdgePoints = 0;

    numNodesAndEdgesArray.resize(11*numTriangles);
    for (i=0; i<numTriangles; i++) {
        const DomainTriangle<ctype>& cT = par->triangles(i);

        int numIntersectionNodes;
        int numTouchingNodes;
        int numInteriorNodes;

        cT.countNodes(numIntersectionNodes, numTouchingNodes, numInteriorNodes);
        //        int numEdges = cT.getNumRegularEdges() - 3;
        int numEdges = cT.getNumRegularEdges();
        numNodesAndEdgesArray[11*i+0] = numIntersectionNodes;
        numNodesAndEdgesArray[11*i+1] = numTouchingNodes;
        numNodesAndEdgesArray[11*i+2] = numInteriorNodes;
        numNodesAndEdgesArray[11*i+3] = numEdges;
        numNodesAndEdgesArray[11*i+4] = cT.patch;

        numNodesAndEdgesArray[11*i+5] = cT.edgePoints[0].size()-2;
        numNodesAndEdgesArray[11*i+6] = cT.edgePoints[1].size()-2;
        numNodesAndEdgesArray[11*i+7] = cT.edgePoints[2].size()-2;

        numNodesAndEdgesArray[11*i+8] = cT.nodes[cT.cornerNode(0)].getNodeNumber();
        numNodesAndEdgesArray[11*i+9] = cT.nodes[cT.cornerNode(1)].getNodeNumber();
        numNodesAndEdgesArray[11*i+10] = cT.nodes[cT.cornerNode(2)].getNodeNumber();

        numNodes += numIntersectionNodes;
        numNodes += numTouchingNodes;
        numNodes += numInteriorNodes;

        numEdgePoints += cT.edgePoints[0].size() + cT.edgePoints[1].size() + cT.edgePoints[2].size() - 6;

        numParamEdges += numEdges;
    }

    ncells = numTriangles + numParamEdges;
    nvertices = numVertices + numNodes;
    /////////////////////////////////////////////////////////////////////
    StaticVector<float,3> imagepos;

    //plane graph on each base grid triangle, saved as a list of nodes and a list of edges.
    int arrayIdx           = 0;
    int edgeArrayIdx       = 0;
    int edgePointsArrayIdx = 0;
    ctype cCoords[3][3];

    triId.resize(numNodes);
    imagePos.resize(numVertices + numNodes);
    nodeType.resize(numVertices + numNodes);
    nodePositions.resize(numNodes);
    domainPositions.resize(numNodes);
    nodeNumber.resize(numNodes);
    parameterEdgeArray.resize(numParamEdges);
    parameterEdgeArrayLocal.resize(numParamEdges);
    edgePointsArray.resize(numEdgePoints);
    for(i = 0; i < numVertices; i++) nodeType[i] = CORNER_NODE;
    for (i=0; i<numTriangles; i++) {
        const DomainTriangle<ctype>& cT = par->triangles(i);

        std::vector<int> newIdx(cT.nodes.size());
        std::vector<int> newIdxlocal(cT.nodes.size());
        for(j = 0; j < 3; j++)
          for(k = 0; k < 3; k++)
            cCoords[j][k] = par->vertices(par->triangles(i).vertices[j])[k];

        // the cornerNode are not saved, because everything about them
        // can be deduced from the base grid
        // the three remaining types are saved separately, in order to avoid
        // explicitly saving the type for each node.
        for(size_t cN = 0; cN < 3; cN++)
        {
            if(!cT.nodes[cN].isCORNER_NODE())printf("error in the corner node indx!\n");
            newIdx[cN] = par->triangles(i).vertices[cN] - numVertices;
        }
        newIdxlocal[cT.cornerNode(0)] = 0;
        newIdxlocal[cT.cornerNode(1)] = 1;
        newIdxlocal[cT.cornerNode(2)] = 2;
        int localArrayIdx = 3;
        for (size_t cN=0; cN<cT.nodes.size(); cN++) {
            if (cT.nodes[cN].isINTERSECTION_NODE()){

                domainPositions[arrayIdx] = cT.nodes[cN].domainPos();
                for(k = 0; k < 3; k++)
                (nodePositions[arrayIdx])[k] = cCoords[0][k]*cT.nodes[cN].domainPos()[0]
                                              +cCoords[1][k]*cT.nodes[cN].domainPos()[1]
                                              +cCoords[2][k]*(1 - cT.nodes[cN].domainPos()[0] - cT.nodes[cN].domainPos()[1]);
                nodeNumber[arrayIdx]     = cT.nodes[cN].getNodeNumber();
                nodeType[arrayIdx+ numVertices] = INTERSECTION_NODE;
                triId[arrayIdx] = i;
                imagePos[arrayIdx+ numVertices] = par->imagePos(i,cN);
                newIdx[cN] = arrayIdx;
                newIdxlocal[cN] = localArrayIdx;
                arrayIdx++;
                localArrayIdx++;
            }
        }

        for (size_t cN=0; cN<cT.nodes.size(); cN++) {

            if (cT.nodes[cN].isTOUCHING_NODE()){

                domainPositions[arrayIdx] = cT.nodes[cN].domainPos();
                for(k = 0; k < 3; k++)
                (nodePositions[arrayIdx])[k] = cCoords[0][k]*cT.nodes[cN].domainPos()[0]
                                              +cCoords[1][k]*cT.nodes[cN].domainPos()[1]
                                              +cCoords[2][k]*(1 - cT.nodes[cN].domainPos()[0] - cT.nodes[cN].domainPos()[1]);
                nodeNumber[arrayIdx]     = cT.nodes[cN].getNodeNumber();
                nodeType[arrayIdx+ numVertices] = TOUCHING_NODE;
                imagePos[arrayIdx+ numVertices] = par->imagePos(i,cN);
                triId[arrayIdx] = i;
                newIdx[cN] = arrayIdx;
                newIdxlocal[cN] = localArrayIdx;
                arrayIdx++;
                localArrayIdx++;
            }
        }

        for (size_t cN=0; cN<cT.nodes.size(); cN++) {

            if (cT.nodes[cN].isINTERIOR_NODE()){

                domainPositions[arrayIdx] = cT.nodes[cN].domainPos();
                for(k = 0; k < 3; k++)
                (nodePositions[arrayIdx])[k] = cCoords[0][k]*cT.nodes[cN].domainPos()[0]
                                              +cCoords[1][k]*cT.nodes[cN].domainPos()[1]
                                              +cCoords[2][k]*(1 - cT.nodes[cN].domainPos()[0] - cT.nodes[cN].domainPos()[1]);
                nodeNumber[arrayIdx]     = cT.nodes[cN].getNodeNumber();
                nodeType[arrayIdx+ numVertices] = INTERIOR_NODE;
                triId[arrayIdx] = i;
                imagePos[arrayIdx+ numVertices] = par->imagePos(i,cN);
                newIdxlocal[cN] = INTERIOR_NODE;
                newIdx[cN] = arrayIdx;
                newIdxlocal[cN] = localArrayIdx;

                arrayIdx++;
                localArrayIdx++;
            }
        }
        // the parameterEdges for this triangle
        typename PlaneParam<ctype>::UndirectedEdgeIterator cE;
        for (cE = cT.firstUndirectedEdge(); cE.isValid(); ++cE){
            if(cE.isRegularEdge())
            {
                //////////////////////////////////////////////////////
                parameterEdgeArray[edgeArrayIdx][0] = newIdx[cE.from()] + numVertices;
                parameterEdgeArray[edgeArrayIdx][1] = newIdx[cE.to()] + numVertices;
/*                float a[3], b[3];
                if( parameterEdgeArray[edgeArrayIdx][0] < numVertices)
                {
                    for(int i = 0; i < 3; i++) a[i] = (baseGridTriArray[parameterEdgeArray[edgeArrayIdx][0]])[i];
                }
                else
                {
                    for(int i = 0; i < 3; i++) a[i] = (nodePositions[parameterEdgeArray[edgeArrayIdx][0] - numVertices])[i];
                }

                if( parameterEdgeArray[edgeArrayIdx][1] < numVertices)
                {
                    for(int i = 0; i < 3; i++) b[i] = (baseGridTriArray[parameterEdgeArray[edgeArrayIdx][1]])[i];
                }
                else
                {
                    for(int i = 0; i < 3; i++) b[i] = (nodePositions[parameterEdgeArray[edgeArrayIdx][1] - numVertices])[i];
                }
                if(parameterEdgeArray[edgeArrayIdx][0] < numVertices && parameterEdgeArray[edgeArrayIdx][1] < numVertices)
                printf("triangle edges\n");
                printf("edges[%d] = [%d %d] first point = [%f  %f  %f] second point = [%f  %f  %f]\n",edgeArrayIdx, parameterEdgeArray[edgeArrayIdx][0], parameterEdgeArray[edgeArrayIdx][1],a[0], a[1], a[2], b[0], b[1], b[2]);
*/

                //////////////////////////////////////////////////////
                parameterEdgeArrayLocal[edgeArrayIdx][0] = newIdxlocal[cE.from()];
                parameterEdgeArrayLocal[edgeArrayIdx][1] = newIdxlocal[cE.to()];
                edgeArrayIdx++;
            }
        }

        // the edgePoints for this triangle
        for (j=0; j<3; j++){
            for (size_t k=1; k<cT.edgePoints[j].size()-1; k++)
            {
                edgePointsArray[edgePointsArrayIdx++] = newIdxlocal[cT.edgePoints[j][k]];
            }
        }
    }
  };

  ///initialize PsurfaceConvert from file
  ///read data array from hdf5 or gmsh file into the data arrary in class PsurfaceConvert
  PsurfaceConvert(const char* filename, bool a)
  {
      if(a) //initialize from hdf5 data
          readHdf5Data(filename);
      else
          readGmsh(filename);
  }   
  
  bool readHdf5Data(const char* filename)
  {
      hid_t file;
      hid_t dataset;
      hid_t filespace;
      hid_t       memspace;
      hsize_t dims[2];
      hsize_t dimz[1];
      herr_t status,status_n;
      int rank;
      int i,j,k;
      file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

      //read params
      dataset = H5Dopen(file,"/Params", H5P_DEFAULT);
      filespace = H5Dget_space(dataset);
      rank = H5Sget_simple_extent_ndims(filespace);
      status_n  = H5Sget_simple_extent_dims(filespace, dimz, NULL);
      memspace = H5Screate_simple(1,dimz,NULL);
      int psurfaceparams[dimz[0]];
      status = H5Dread(dataset, H5T_NATIVE_INT, memspace, filespace,
      H5P_DEFAULT, psurfaceparams);
      numVertices = psurfaceparams[0];
      numNodes = psurfaceparams[1];
      numTriangles = psurfaceparams[2];
      numParamEdges = psurfaceparams[3];
      nvertices = numVertices + numNodes;
      ncells = numTriangles + numParamEdges;
      H5Dclose(dataset);
      H5Sclose(filespace);
      H5Sclose(memspace);
      //read xyz coordinate of vertices
      baseGridVertexCoordsArray.resize(numVertices);
      nodePositions.resize(numNodes);
      for(i = 0; i < 3; i++)
      {
          if(i == 0)  dataset = H5Dopen(file,"/X", H5P_DEFAULT);
          else if(i == 1) dataset = H5Dopen(file,"/Y", H5P_DEFAULT);
          else dataset = H5Dopen(file,"/Z", H5P_DEFAULT);
          filespace = H5Dget_space(dataset);
          rank = H5Sget_simple_extent_ndims(filespace);
          status_n  = H5Sget_simple_extent_dims(filespace, dimz, NULL);

          memspace = H5Screate_simple(1,dimz,NULL);
          ctype coords[dimz[0]];
          status = H5Dread(dataset, H5T_NATIVE_FLOAT, memspace, filespace,
              H5P_DEFAULT, coords);
          H5Dclose(dataset);
          H5Sclose(filespace);
          H5Sclose(memspace);
          for(j = 0; j < numVertices ;j++) (baseGridVertexCoordsArray[j])[i] = coords[j];
          for(j = 0; j < numNodes; j++) (nodePositions[j])[i] = coords[numVertices + j];
        }

        //triangle
        dataset = H5Dopen(file,"/Topo", H5P_DEFAULT);
        filespace = H5Dget_space(dataset);
        rank = H5Sget_simple_extent_ndims(filespace);
        status_n  = H5Sget_simple_extent_dims(filespace, dimz, NULL);

        memspace = H5Screate_simple(1,dimz,NULL);
        int tri[dimz[0]];
        status = H5Dread(dataset, H5T_NATIVE_INT, memspace, filespace,
        H5P_DEFAULT, tri);
        H5Dclose(dataset);
        H5Sclose(filespace);
        H5Sclose(memspace);
        baseGridTriArray.resize(numTriangles);
        for(j = 0; j < numTriangles;j++)
          for(i = 0; i < 3; i++)
            (baseGridTriArray[j])[i] = tri[4*j+ 1 +i];

        parameterEdgeArray.resize(numParamEdges);
        for(j = 0; j < numParamEdges; j++)
          for(i = 0; i < 2; i++)
            (parameterEdgeArray[j])[i] = tri[4*j + 2 + i + 4*numTriangles];

          //Parameter Edge
          dataset = H5Dopen(file,"/LocalParamEdge", H5P_DEFAULT);
          filespace = H5Dget_space(dataset);
          rank = H5Sget_simple_extent_ndims(filespace);
          status_n  = H5Sget_simple_extent_dims(filespace, dimz, NULL);
          memspace = H5Screate_simple(1,dimz,NULL);
          int lparam[dimz[0]];
          status = H5Dread(dataset, H5T_NATIVE_INT, memspace, filespace,
          H5P_DEFAULT,lparam);
          H5Dclose(dataset);
          H5Sclose(filespace);
          H5Sclose(memspace);
          parameterEdgeArrayLocal.resize(numParamEdges);
          for(j = 0; j < dimz[0]/2;j++)
              for(i = 0; i < 2; i++)
                  (parameterEdgeArrayLocal[j])[i] = lparam[2*j + i];

          //nodetype
          dataset = H5Dopen(file,"/NodeType", H5P_DEFAULT);
          filespace = H5Dget_space(dataset);
          rank = H5Sget_simple_extent_ndims(filespace);
          status_n  = H5Sget_simple_extent_dims(filespace, dimz, NULL);

          memspace = H5Screate_simple(1,dimz,NULL);
          int nodetype[dimz[0]];
          status = H5Dread(dataset, H5T_NATIVE_INT, memspace, filespace,
             H5P_DEFAULT,nodetype);
          H5Dclose(dataset);
          H5Sclose(filespace);
          H5Sclose(memspace);
          nodeType.resize(nvertices);
          for(i = 0; i < nvertices;i++) nodeType[i] = nodetype[i];

          //nodeNumber
          dataset = H5Dopen(file,"/NodeNumber",H5P_DEFAULT);
          filespace = H5Dget_space(dataset);
          rank = H5Sget_simple_extent_ndims(filespace);
          status_n  = H5Sget_simple_extent_dims(filespace, dimz, NULL);
          memspace = H5Screate_simple(1,dimz,NULL);
          int nodenumber[dimz[0]];
          status = H5Dread(dataset, H5T_NATIVE_INT, memspace, filespace,
             H5P_DEFAULT,nodenumber);
          H5Dclose(dataset);
          H5Sclose(filespace);
          H5Sclose(memspace);
          nodeNumber.resize(numNodes);
          for(i = 0; i < numNodes;i++) nodeNumber[i] = nodenumber[i];

          //iPos
          dataset = H5Dopen(file,"/iPos", H5P_DEFAULT);
          filespace = H5Dget_space(dataset);
          rank = H5Sget_simple_extent_ndims(filespace);
          status_n  = H5Sget_simple_extent_dims(filespace, dimz, NULL);

          memspace = H5Screate_simple(1,dimz,NULL);
          float coords[dimz[0]];
          status = H5Dread(dataset, H5T_NATIVE_FLOAT, memspace, filespace,
           H5P_DEFAULT, coords);
          H5Dclose(dataset);
          H5Sclose(filespace);
          H5Sclose(memspace);
          iPos.resize(dimz[0]/3);
          for(j = 0; j < dimz[0]/3; j++)
          {
              (iPos[j])[0] = coords[3*j];
              (iPos[j])[1] = coords[3*j+1];
              (iPos[j])[2] = coords[3*j+2];
          }

          //nodeNodesAndEdgesArray
          dataset = H5Dopen(file,"/numNodesAndEdgesArray", H5P_DEFAULT);
          filespace = H5Dget_space(dataset);
          rank = H5Sget_simple_extent_ndims(filespace);
          status_n  = H5Sget_simple_extent_dims(filespace, dimz, NULL);

          memspace = H5Screate_simple(1,dimz,NULL);
          int nodearray[dimz[0]];
          status = H5Dread(dataset, H5T_NATIVE_INT, memspace, filespace,
             H5P_DEFAULT,nodearray);
          H5Dclose(dataset);
          H5Sclose(filespace);
          H5Sclose(memspace);
          numNodesAndEdgesArray.resize(11*numTriangles);
          for(i = 0; i < 11*numTriangles;i++)
              numNodesAndEdgesArray[i] = nodearray[i];

          //triangle index
          dataset = H5Dopen(file,"/TriangleIndx", H5P_DEFAULT);
          filespace = H5Dget_space(dataset);
          rank = H5Sget_simple_extent_ndims(filespace);
          status_n  = H5Sget_simple_extent_dims(filespace, dimz, NULL);

          memspace = H5Screate_simple(1,dimz,NULL);
          int tridx[dimz[0]];
          status = H5Dread(dataset, H5T_NATIVE_INT, memspace, filespace,
               H5P_DEFAULT,tridx);
          H5Dclose(dataset);
          H5Sclose(filespace);
          H5Sclose(memspace);
          triId.resize(numNodes);
          for(i = 0; i < numNodes;i++)
              triId[i] = tridx[i];

          //local position of nodes on triangle
          dataset = H5Dopen(file,"/LocalNodePosition", H5P_DEFAULT);
          filespace = H5Dget_space(dataset);
          rank = H5Sget_simple_extent_ndims(filespace);
          status_n  = H5Sget_simple_extent_dims(filespace, dimz, NULL);
          memspace = H5Screate_simple(1,dimz,NULL);
          float localpos[dimz[0]];
          status = H5Dread(dataset, H5T_NATIVE_FLOAT, memspace, filespace,
                H5P_DEFAULT,localpos);
          H5Dclose(dataset);
          H5Sclose(filespace);
          H5Sclose(memspace);
          domainPositions.resize(numNodes);
          for(i = 0; i < numNodes;i++)
            for(j = 0; j < 2;j++)
              (domainPositions[i])[j] = localpos[2*i + j];

          //edgepointsarray
          dataset = H5Dopen(file,"/EdgePoints", H5P_DEFAULT);
          filespace = H5Dget_space(dataset);
          rank = H5Sget_simple_extent_ndims(filespace);
          status_n  = H5Sget_simple_extent_dims(filespace, dimz, NULL);

          memspace = H5Screate_simple(1,dimz,NULL);
          int edgep[dimz[0]];
          status = H5Dread(dataset, H5T_NATIVE_INT, memspace, filespace,
             H5P_DEFAULT,edgep);
          H5Dclose(dataset);
          H5Sclose(filespace);
          H5Sclose(memspace);
          edgePointsArray.resize(dimz[0]);
          for(i = 0; i < dimz[0];i++)
              edgePointsArray[i] = edgep[i];

          //patches
          dataset = H5Dopen(file,"/patches_vec", H5P_DEFAULT);
          filespace = H5Dget_space(dataset);
          rank = H5Sget_simple_extent_ndims(filespace);
          status_n  = H5Sget_simple_extent_dims(filespace, dimz, NULL);

          memspace = H5Screate_simple(1,dimz,NULL);
          int patch[dimz[0]];
          status = H5Dread(dataset, H5T_NATIVE_INT, memspace, filespace,
             H5P_DEFAULT,patch);
          H5Dclose(dataset);
          H5Sclose(filespace);
          H5Sclose(memspace);
          patches.resize(dimz[0]/3);
          for(i = 0; i < dimz[0]/3;i++)
          {
            (patches[i]).innerRegion = patch[3*i];
            (patches[i]).outerRegion = patch[3*i + 1];
            (patches[i]).boundaryId = patch[3*i + 2];
          }
          H5Fclose(file);
          return 0;
  }

  bool readGmsh(const char* filename)
  {
      FILE* file = fopen(filename,"r");
      int number_of_real_vertices = 0;
      int element_count = 0;
      // process header
      double version_number;
      int file_type, data_size;
      char buf[512];

      readfile(file,1,"%s\n",buf);
      if (strcmp(buf,"$MeshFormat")!=0)
          throw(std::runtime_error("expected $MeshFormat in first line\n"));

      readfile(file,3,"%lg %d %d\n",&version_number,&file_type,&data_size);
      if( (version_number < 2.0) || (version_number > 2.2) )
          throw(std::runtime_error("can only read Gmsh version 2 files\n"));
      readfile(file,1,"%s\n",buf);
      if (strcmp(buf,"$EndMeshFormat")!=0)
          throw(std::runtime_error("expected $EndMeshFormat\n"));

      // node section
      int number_of_nodes;
      readfile(file,1,"%s\n",buf);
      if (strcmp(buf,"$Nodes")!=0)
        throw(std::runtime_error("expected $Nodes\n"));
      readfile(file,1,"%d\n",&number_of_nodes);

      std::vector<StaticVector<int, 3> > triArray;
      std::vector<StaticVector<ctype,3> > coordsArray;
       
      // read nodes
      int id;
      double x[3];
      for( int i = 1; i <= number_of_nodes; ++i )
      {
          readfile(file,4, "%d %lg %lg %lg\n", &id, &x[0], &x[1], &x[2] );
          if( id != i ) 
              throw(std::runtime_error("id does not match in reading gmsh"));
            
          // just store node position
          StaticVector<ctype,3> vertex;
          for(int j = 0 ; j < 3; j++)
              vertex[j] = x[j];
          coordsArray.push_back(vertex);
      }
      

      readfile(file,1,"%s\n",buf);
      if (strcmp(buf,"$EndNodes")!=0)
        printf("expected $EndNodes\n");

      // element section
      readfile(file,1,"%s\n",buf);
      if (strcmp(buf,"$Elements")!=0)
          throw(std::runtime_error("expected $Elements\n"));
      int number_of_elements;
      readfile(file,1,"%d\n",&number_of_elements);

      for (int i=1; i<=number_of_elements; i++)
      {
        int id, elm_type, number_of_tags;
        readfile(file,3,"%d %d %d ",&id,&elm_type,&number_of_tags);
        for (int k=1; k<=number_of_tags; k++)
        {
          int blub;
          readfile(file,1,"%d ",&blub);
        }
      
        if(elm_type != 2) 
        {
            skipline(file);
            continue;
        }
        std::string formatString = "%d";
        for (int i= 1; i< 3; i++)
          formatString += " %d";
        formatString += "\n";

        // '10' is the largest number of dofs we may encounter in a .msh file
        StaticVector<int, 3> elementDofs(3);

        readfile(file, 3, formatString.c_str(), &(elementDofs[0]),&(elementDofs[1]),&(elementDofs[2]));
        triArray.push_back(elementDofs);
      }

      //remove vetices in the code that is not the corner of triangle
      
      bool nodeInTri[number_of_nodes];
      int  newNodeIndex[number_of_nodes];
      for(int i = 0; i < number_of_nodes; i++) nodeInTri[i] = 0;
      for(int i = 0; i < triArray.size(); i++) 
      {
          nodeInTri[(triArray[i])[0] - 1] = 1;
          nodeInTri[(triArray[i])[1] - 1] = 1;
          nodeInTri[(triArray[i])[2] - 1] = 1;
      }
      int newIndx = 0;
      for(int i = 0; i < number_of_nodes; i++)
      {
          if(!nodeInTri[i]) 
              newNodeIndex[i] = -1;
          else
          {
              newNodeIndex[i] = newIndx;
              newIndx++;
          }
      }
      
      //push the nodes and element into baseGridVertexCoordsArray and baseGridTriArray
      for(int i = 0; i < number_of_nodes; i++)
          if(nodeInTri[i])
            baseGridVertexCoordsArray.push_back(coordsArray[i]);
      
      baseGridTriArray.resize(triArray.size());
      for(int i = 0; i < triArray.size(); i++) 
      {
          (baseGridTriArray[i])[0] = newNodeIndex[(triArray[i])[0] - 1];
          (baseGridTriArray[i])[1] = newNodeIndex[(triArray[i])[1] - 1];
          (baseGridTriArray[i])[2] = newNodeIndex[(triArray[i])[2] - 1];
      }
      numVertices = baseGridVertexCoordsArray.size();
      numTriangles = baseGridTriArray.size();
      return 0;
  }
 
  bool initPsurface(PSurface<2,ctype>* psurf, Surface* surf, bool baseTriangle)
  {
      if(baseTriangle)
          return initBaseCasePSurface(psurf, surf);
      else
          return initCompletePSurface(psurf, surf);
  }

  ///Create PSurface which only have BaseGridVertex and BaseTriangle
  bool initBaseCasePSurface(PSurface<2,ctype>* psurf, Surface* surf)
  {
    PSurfaceFactory<2,ctype> factory(psurf);
    ///(Assume) Target surface already exists
    factory.setTargetSurface(surf);
      
    //set param for psurface
    AmiraMesh am;
    am.parameters.set("ContentType", "Parametrization");
    am.parameters.remove(am.parameters[0]);
    psurf->getPaths(am.parameters);

    //patches
    PSurface<2, float>::Patch Patch; //= new PSurface<2, float>::Patch;
    Patch.innerRegion = 0;
    Patch.outerRegion = 1;
    Patch.boundaryId =  2;
    psurf->patches.push_back(Patch);
    Patch.innerRegion = 0;
    Patch.outerRegion = 1;
    Patch.boundaryId =  1;
    psurf->patches.push_back(Patch);

    ///insert vertex
    StaticVector<ctype,3> newVertex;
    for(int i = 0; i < numVertices; i++)
    {
      for(int j = 0; j < 3; j++) newVertex[j] = (baseGridVertexCoordsArray[i])[j];
      factory.insertVertex(newVertex);
    }

    ///insert image node position
    psurf->iPos.resize(numVertices);
    for (size_t i=0; i< numVertices; i++)
      for (int j=0; j<3; j++)
        psurf->iPos[i][j] =(baseGridVertexCoordsArray[i])[j];

    ///insert trianlges and the plain graph onto it.
    int edgeCounter=0, edgePointCounter=0;
    int nodeArrayIdx = 0;
    for (int i=0; i<numTriangles; i++){
      std::tr1::array<unsigned int, 3> triangleVertices = {(baseGridTriArray[i])[0],(baseGridTriArray[i])[1],(baseGridTriArray[i])[2]};
    int newTriIdx = factory.insertSimplex(triangleVertices);
    psurf->triangles(newTriIdx).patch = 0;
    /// get the parametrization on this triangle

    ///nodes
    psurf->triangles(newTriIdx).nodes.resize(3);
    int nodenumber;
    // three corner nodes
    StaticVector<ctype,2> domainPos(1, 0);
    nodenumber =  0;
    psurf->triangles(newTriIdx).nodes[0].setValue(domainPos, nodenumber, Node<ctype>::CORNER_NODE);
    psurf->triangles(newTriIdx).nodes[0].makeCornerNode(0, nodenumber);

    domainPos = StaticVector<ctype,2>(0, 1);
    nodenumber = 1;
    psurf->triangles(newTriIdx).nodes[1].setValue(domainPos, nodenumber, Node<ctype>::CORNER_NODE);
    psurf->triangles(newTriIdx).nodes[1].makeCornerNode(1, nodenumber);

    domainPos = StaticVector<ctype,2>(0, 0);
    nodenumber = 2;
    psurf->triangles(newTriIdx).nodes[2].setValue(domainPos, nodenumber, Node<ctype>::CORNER_NODE);
    psurf->triangles(newTriIdx).nodes[2].makeCornerNode(2, nodenumber);

    /// the parameterEdges
    for(int j = 0; j < 3; j++)
      psurf->triangles(newTriIdx).addEdge(j, (j+1)%3);
    
    /// the edgePoints arrays on each triangular edge
    for (int j=0; j<3; j++){
      psurf->triangles(newTriIdx).edgePoints[j].resize(2);
      psurf->triangles(newTriIdx).edgePoints[j][0]     = j;
      psurf->triangles(newTriIdx).edgePoints[j].back() = (j+1)%3;
    }
    }
    psurf->hasUpToDatePointLocationStructure = false;
    psurf->setupOriginalSurface();
    return true;
  } 
 
  ///Creat Psurface
  bool initCompletePSurface(PSurface<2,ctype>* psurf, Surface* surf)
  {
    /// Create PSurface factory
    PSurfaceFactory<2,ctype> factory(psurf);
    ///(Assume) Target surface already exists
    factory.setTargetSurface(surf);

    //set param for psurface
    AmiraMesh am;
    am.parameters.set("ContentType", "Parametrization");
    am.parameters.remove(am.parameters[0]);
    psurf->getPaths(am.parameters);

    //patches
    PSurface<2, float>::Patch Patch; //= new PSurface<2, float>::Patch;
    for (size_t i=0; i<patches.size(); i++)
    {
      Patch.innerRegion = (patches[i]).innerRegion;
      Patch.outerRegion = (patches[i]).outerRegion ;
      Patch.boundaryId = (patches[i]).boundaryId;
      psurf->patches.push_back(Patch);
    }

    ///insert vertex
    StaticVector<ctype,3> newVertex;
    for(int i = 0; i < numVertices; i++)
    {
       for(int j = 0; j < 3; j++) newVertex[j] = (baseGridVertexCoordsArray[i])[j];
       factory.insertVertex(newVertex);
    }

    ///insert image node position
    psurf->iPos.resize(iPos.size());
    for (size_t i=0; i<psurf->iPos.size(); i++)
      for (int j=0; j<3; j++)
        psurf->iPos[i][j] =((iPos[i])[j]);



    ///insert trianlges and the plain graph onto it.
    int edgeCounter=0, edgePointCounter=0;
    int nodeArrayIdx = 0;
    for (int i=0; i<numTriangles; i++){
      std::tr1::array<unsigned int, 3> triangleVertices = {(baseGridTriArray[i])[0],(baseGridTriArray[i])[1],(baseGridTriArray[i])[2]};
      int newTriIdx = factory.insertSimplex(triangleVertices);
    psurf->triangles(newTriIdx).patch = numNodesAndEdgesArray[11*i+4];
    /// get the parametrization on this triangle
    int numIntersectionNodes = numNodesAndEdgesArray[11*i+0];
    int numTouchingNodes     = numNodesAndEdgesArray[11*i+1];
    int numInteriorNodes     = numNodesAndEdgesArray[11*i+2];
    int numParamEdges        = numNodesAndEdgesArray[11*i+3];

    ///nodes
    psurf->triangles(newTriIdx).nodes.resize(numIntersectionNodes + numTouchingNodes + numInteriorNodes + 3);
    int nodenumber;
    // three corner nodes
    StaticVector<ctype,2> domainPos(1, 0);
    nodenumber =  numNodesAndEdgesArray[11*i+8];
    psurf->triangles(newTriIdx).nodes[0].setValue(domainPos, nodenumber, Node<ctype>::CORNER_NODE);
    psurf->triangles(newTriIdx).nodes[0].makeCornerNode(0, nodenumber);

    domainPos = StaticVector<ctype,2>(0, 1);
    nodenumber = numNodesAndEdgesArray[11*i+9];
    psurf->triangles(newTriIdx).nodes[1].setValue(domainPos, nodenumber, Node<ctype>::CORNER_NODE);
    psurf->triangles(newTriIdx).nodes[1].makeCornerNode(1, nodenumber);

    domainPos = StaticVector<ctype,2>(0, 0);
    nodenumber = numNodesAndEdgesArray[11*i+10];
    psurf->triangles(newTriIdx).nodes[2].setValue(domainPos, nodenumber, Node<ctype>::CORNER_NODE);
    psurf->triangles(newTriIdx).nodes[2].makeCornerNode(2, nodenumber);

     int nodeCounter = 3;
    ///the intersection nodes
    for (int j=0; j<numIntersectionNodes; j++, nodeCounter++, nodeArrayIdx++){
      nodenumber = nodeNumber[nodeArrayIdx];
      psurf->triangles(newTriIdx).nodes[nodeCounter].setValue(domainPositions[nodeArrayIdx], nodenumber, Node<ctype>::INTERSECTION_NODE);

    }

    // the touching nodes
    for (int j=0; j<numTouchingNodes; j++, nodeCounter++, nodeArrayIdx++){
      int nodenumber    = nodeNumber[nodeArrayIdx];
      psurf->triangles(newTriIdx).nodes[nodeCounter].setValue(domainPositions[nodeArrayIdx], nodenumber, Node<ctype>::TOUCHING_NODE);
    }

    // the interior nodes
    for (int j=0; j<numInteriorNodes; j++, nodeCounter++, nodeArrayIdx++){
      int nodenumber    = nodeNumber[nodeArrayIdx];
      psurf->triangles(newTriIdx).nodes[nodeCounter].setValue(domainPositions[nodeArrayIdx], nodenumber, Node<ctype>::INTERIOR_NODE);
    }

    /// the parameterEdges
    for(int j = 0; j < numParamEdges;j++, edgeCounter++)
    {
      psurf->triangles(newTriIdx).addEdge(parameterEdgeArrayLocal[edgeCounter][0],parameterEdgeArrayLocal[edgeCounter][1]);
    }
    /// the edgePoints arrays on each triangular edge
    for (int j=0; j<3; j++){
      psurf->triangles(newTriIdx).edgePoints[j].resize(numNodesAndEdgesArray[11*i+5+j] + 2);
      psurf->triangles(newTriIdx).edgePoints[j][0]     = j;
      psurf->triangles(newTriIdx).edgePoints[j].back() = (j+1)%3;


      for (int k=0; k<numNodesAndEdgesArray[11*i+5+j]; k++){
          psurf->triangles(newTriIdx).edgePoints[j][k+1] = edgePointsArray[edgePointCounter];
          edgePointCounter++;
      }
    }
    }
    psurf->hasUpToDatePointLocationStructure = false;
    psurf->setupOriginalSurface();
    return true;
  };

  void skipline(FILE * file)
  {
    int c;
    do {
      c = fgetc(file);
    } while(c != '\n' && c != EOF);
  }

  void readfile(FILE * file, int cnt, const char * format,
                  void* t1, void* t2 = 0, void* t3 = 0, void* t4 = 0,
                  void* t5 = 0, void* t6 = 0, void* t7 = 0, void* t8 = 0,
                  void* t9 = 0, void* t10 = 0)
  {
      off_t pos = ftello(file);
      int c = fscanf(file, format, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10);
      if (c != cnt)
          throw(std::runtime_error("error in readfile\n"));
  }
};
}
#endif

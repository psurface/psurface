//the code to write psurface object into vtk file
//vtuwrite
//#ifndef HDF5IO_H
//#define HDF5IO_H
#include <vector>
#include "StaticVector.h"
#include "Domains.h"
#include <stdio.h>
#include <stdlib.h>
#include <hdf5.h>
#include <amiramesh/HxParamBundle.h>
#include <amiramesh/HxParamBase.h>
#define PSURFACE_STANDALONE

#ifdef PSURFACE_STANDALONE
#include "TargetSurface.h"
#else
#include <hxsurface/Surface.h>
#endif
enum NodeTypes {INTERIOR_NODE=0,
                INTERSECTION_NODE=1,
                CORNER_NODE=2,
                TOUCHING_NODE=3,
                GHOST_NODE=4};
using namespace psurface;
namespace psurface{
/**
 *@brief class to write psurface subject into HDF5 format data
 *@tParams ctype the type used for the coordinates of psurface.
 *This class produce two files.
 *One is the hdf5 data, which contains the information of the triangles,
 *edges and nodes, stored in related arrays.
 *Another one is the XDMF file, which contains information on how to
 *organize the data in hdf5 file so it could be read by paraview or
 *Visit.
 */
template<class ctype,int dim>
class Hdf5IO{
   typedef typename PSurface<2,ctype>::Patch PATCH;
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

  /**@name data structure to store psurface nodes, triangles and edges on target surface*/
  //@{
  //vertices
  std::vector<StaticVector<ctype,3> > VertexCoordsArray;
  ///triangles
  std::vector<StaticVector<int, 3> > TriArray;
  int nV;
  int nTri;
  //@}
  //
  ///read data array from hdf5 data structure into the data arrary in hdf5IO
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

  ///write the data array into hdf5 data structure
  bool writeHdf5Data(const char*filename)
  {
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
  void writeXdmf(const char* xdf_filename, const char* hdf_filename)
  {
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
  };

  ///Read information from the psurface object and stored it in local vector containers in hdf5Writer
  void update()
  {
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
      imagePos.push_back(imagearray[i]);

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



  bool initFromHDF5(PSurface<2,ctype>* psurf, Surface* surf, const char* filename)
  {
    readHdf5Data(filename);
    /// Create PSurface factory
    PSurfaceFactory<2,ctype> factory(psurf);
    ///(Assume) Target surface already exists
    factory.setTargetSurface(surf);

    //set param for psurface

/*  {
    HxParamBundle parameters("Parameters");
//    HxParamBundle* paramBaseList = new  HxParamBundle;
//    HxParamBase paramBaseList(0,"Materials");
//    HxParamBase paramBaseList;
    HxParamBundle *paramBaseList1 = new HxParamBundle;
    paramBaseList1->setName("Materials");
    parameters.insert(paramBaseList1);

    HxParamBundle* paramBaseList2 = new HxParamBundle;
    paramBaseList2->setName("BoundaryIds");
    parameters.insert(paramBaseList2);

//    HxParamBundle* paramBaseList3 = new HxParamBundle;
//    paramBaseList3->setName("ParameterInfo");
//    parameters.insert(paramBaseList3);

//    parameters.move(*paramBaseList3);

    psurf->getPaths(parameters);
 //   delete &parameters;
 //   delete paramBaseList1;
 //   delete paramBaseList2;
 //   delete paramBaseList3;
  }
  */
    AmiraMesh am;
    am.parameters.set("ContentType", "Parametrization");
    am.parameters.remove(am.parameters[0]);
    psurf->getPaths(am.parameters);
//    {
//        HxParameter* p = new HxParameter(buf, tmp.size(), &tmp[0]);
//        char buf[64];
//        std::vector<int> tmp;
//        parameters.insert(p);
//    }*/
    //    paramBase =
//    paramBaseList = new HxParamBundle;
//    paramBaseList->insert(paramBase);

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
};
template <class ctype>
class Hdf5IO<ctype,1>
{
  public:
  PSurface<1,ctype>* par1d;
  //stuff on original surface
  std::vector<StaticVector<ctype,2> > cornerNodes;
  std::vector<StaticVector<ctype,2> > innerNodes;
  std::vector<StaticVector<int,2> > edges;
  std::vector<StaticVector<ctype,2> > domainNorm;
  //stuff on target surface
  std::vector<StaticVector<ctype, 2> > targetNodes;
  std::vector<std::tr1::array<int, 2> > targetSegments;
  int numNodes,numCornerNodes,numInnerNodes;
  int numEdges;
  int numTargetSegments, numTargetNodes;

  ///write the data array into hdf5 data structure
  void writeHdf5Data1d()
  {
    hid_t     file_id;
    file_id = H5Fcreate("psurface1d.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    int i;
    // Create the domain coordinate data.
    float *x = (float *) malloc(numNodes* sizeof(float));
    float *y = (float *) malloc(numNodes * sizeof(float));
    std::vector<StaticVector<float,2> >::iterator pv;
    pv = cornerNodes.begin();
    for(i = 0; i < numCornerNodes; i++)
    {
      x[i] = (*pv)[0];
      y[i] = (*pv)[1];
      pv++;
    }
    pv = innerNodes.begin();
    for(i = 0; i < numInnerNodes;i++)
    {
      x[i] = (*pv)[0];
      y[i] = (*pv)[1];
      pv++;
    }
    //Creat connection array
    int  *edge_array = (int *) malloc(numEdges * 2*sizeof(int));
    std::vector<StaticVector<int, 2> > ::iterator pt;
    pt = edges.begin();
    for(i = 0; i < numEdges;i++)
    {
      edge_array[2*i] = (*pt)[0];
      edge_array[2*i+1] = (*pt)[1];
      pt++;
    }
    // Create the image coordinate data.
    float *ix = (float *) malloc(numTargetNodes * sizeof(float));
    float *iy = (float *) malloc(numTargetNodes * sizeof(float));
    pv = targetNodes.begin();
    for(i = 0; i < numTargetNodes; i++)
    {
      x[i] = (*pv)[0];
      y[i] = (*pv)[1];
      pv++;
    }
    //Creat connection array
    int  *image_edge_array = (int *) malloc(numTargetSegments * 2*sizeof(int));
    pt = targetSegments.begin();
    for(i = 0; i < numTargetSegments;i++)
    {
      image_edge_array[2*i] = (*pt)[0];
      image_edge_array[2*i+1] = (*pt)[1];
      pt++;
    }

    // Write the data file.
    hid_t     dataset_id, dataspace_id,datatype;
    hsize_t   dims[] ={numNodes};
    herr_t    status;
    // Write separate coordinate arrays for the x and ycoordinates.
    dims[0] = numNodes;
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
    //Write separate coordinate arrays for the x and ycoordinates of image surface.
    dims[0] = numTargetNodes;
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dataset_id = H5Dcreate2(file_id, "/ImageX", H5T_NATIVE_FLOAT, dataspace_id,
         H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, ix);
    status = H5Sclose(dataspace_id);
    status = H5Dclose(dataset_id);
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dataset_id = H5Dcreate2(file_id, "/ImageY", H5T_NATIVE_FLOAT, dataspace_id,
    H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, iy);
    status = H5Sclose(dataspace_id);
    status = H5Dclose(dataset_id);
    //topology information of original surface
    dims[0] = numEdges*2;
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dataset_id = H5Dcreate(file_id, "/Edges", H5T_NATIVE_INT, dataspace_id,
    H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
    H5P_DEFAULT, edge_array);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);
    //topology information of image surface
    dims[0] = numTargetSegments*2;
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dataset_id = H5Dcreate(file_id, "/TargetSegments", H5T_NATIVE_INT, dataspace_id,
    H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
    H5P_DEFAULT, image_edge_array);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);
    // Free the data.
    free(x);
    free(y);
    free(ix);
    free(iy);
    free(edge_array);
    free(image_edge_array);
    status = H5Fclose(file_id);
  };
  ///writhe the xdmf file which store the structure information of hdf5 file.
  void writedomainxdmf()
  {
    FILE *xmf = 0;
    xmf = fopen("domainxdmf2d.xmf", "w");
    fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
    fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
    fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
    fprintf(xmf, "<Domain>\n");
    fprintf(xmf, "<Grid Name=\"psurface\" GridType=\"Uniform\">\n");
    fprintf(xmf, "<Topology TopologyType=\"polyline\" NumberOfElements=\"%d\">\n", numEdges);
    fprintf(xmf, "<DataItem Dimensions = \"%d %d\" NumberType=\"Float\" Format=\"HDF\">\n", numEdges, 2);
    fprintf(xmf, "psurface.h5:/Edges\n");
    fprintf(xmf, "</DataItem>\n");
    fprintf(xmf, "</Topology>\n");
    fprintf(xmf, "<Geometry GeometryType=\"X_Y\">\n");
    fprintf(xmf, "<DataItem Dimensions=\"%d\" NumberType=\"Float\" Format=\"HDF\">\n", numNodes);
    fprintf(xmf, "psurface.h5:/X\n");
    fprintf(xmf, "</DataItem>\n");
    fprintf(xmf, "<DataItem Dimensions=\"%d\" NumberType=\"Float\" Format=\"HDF\">\n", numNodes);
    fprintf(xmf, "psurface.h5:/Y\n");
    fprintf(xmf, "</DataItem>\n");
    fprintf(xmf, "</Geometry>\n");
    fprintf(xmf, "</Grid>\n");
    fprintf(xmf, "</Domain>\n");
    fprintf(xmf, "</Xdmf>\n");
    fclose(xmf);
  };

  ///writhe the xdmf file which store the structure information of hdf5 file.
  void writeimagexdmf()
  {
    FILE *xmf = 0;
    xmf = fopen("imagexdmf2d.xmf", "w");
    fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
    fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
    fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
    fprintf(xmf, "<Domain>\n");
    fprintf(xmf, "<Grid Name=\"psurface\" GridType=\"Uniform\">\n");
    fprintf(xmf, "<Topology TopologyType=\"polyline\" NumberOfElements=\"%d\">\n", numTargetSegments);
    fprintf(xmf, "<DataItem Dimensions = \"%d %d\" NumberType=\"Float\" Format=\"HDF\">\n", numTargetSegments, 2);
    fprintf(xmf, "psurface.h5:/Edges\n");
    fprintf(xmf, "</DataItem>\n");
    fprintf(xmf, "</Topology>\n");
    fprintf(xmf, "<Geometry GeometryType=\"X_Y\">\n");
    fprintf(xmf, "<DataItem Dimensions=\"%d\" NumberType=\"Float\" Format=\"HDF\">\n", numTargetNodes);
    fprintf(xmf, "psurface.h5:/X\n");
    fprintf(xmf, "</DataItem>\n");
    fprintf(xmf, "<DataItem Dimensions=\"%d\" NumberType=\"Float\" Format=\"HDF\">\n", numTargetNodes);
    fprintf(xmf, "psurface.h5:/Y\n");
    fprintf(xmf, "</DataItem>\n");
    fprintf(xmf, "</Geometry>\n");
    fprintf(xmf, "</Grid>\n");
    fprintf(xmf, "</Domain>\n");
    fprintf(xmf, "</Xdmf>\n");
    fclose(xmf);
  };

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
};
//#endif

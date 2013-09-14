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
#include "ppsurface_convert.h"
//#ifdef PSURFACE_STANDALONE
#include "TargetSurface.h"
//#else
//#include <hxsurface/Surface.h>
//#endif

using namespace psurface;

void writeIntDataToFile(hid_t* file_id, hid_t* dataset_id, hid_t* dataspace_id, hid_t* datatype, hsize_t* dims, herr_t* status, int dimen, const char*  name, int* address)
{
    *dataspace_id = H5Screate_simple(1, dims, NULL);
    *dataset_id = H5Dcreate(*file_id, name, H5T_NATIVE_INT, *dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    *status = H5Dwrite(*dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, address);
    *status = H5Dclose(*dataset_id);
    *status = H5Sclose(*dataspace_id);
}

void writeIntDataToFile(hid_t* file_id, hid_t* dataset_id, hid_t* dataspace_id, hid_t* datatype, hsize_t* dims, herr_t* status, int dimen, const char*  name, int** address)
{
    *dataspace_id = H5Screate_simple(2, dims, NULL);
    *dataset_id = H5Dcreate(*file_id, name, H5T_NATIVE_INT, *dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    *status = H5Dwrite(*dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, address);
    *status = H5Dclose(*dataset_id);
    *status = H5Sclose(*dataspace_id);
}

void writeFloatDataToFile(hid_t* file_id, hid_t* dataset_id, hid_t* dataspace_id, hid_t* datatype, hsize_t* dims, herr_t* status, int dimen, const char* name, float* address)
  {
    *dataspace_id = H5Screate_simple(1, dims, NULL);
    *dataset_id = H5Dcreate(*file_id, name, H5T_NATIVE_FLOAT, *dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    *status = H5Dwrite(*dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, address);
    *status = H5Dclose(*dataset_id);
    *status = H5Sclose(*dataspace_id);
  }

void writeFloatDataToFile(hid_t* file_id, hid_t* dataset_id, hid_t* dataspace_id, hid_t* datatype, hsize_t* dims, herr_t* status, int dimen, const char* name, float** address)
  {
    *dataspace_id = H5Screate_simple(2, dims, NULL);
    *dataset_id = H5Dcreate(*file_id, name, H5T_NATIVE_FLOAT, *dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    *status = H5Dwrite(*dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, address);
    *status = H5Dclose(*dataset_id);
    *status = H5Sclose(*dataspace_id);
  }

  void writeDoubletDataToFile(hid_t* file_id, hid_t* dataset_id, hid_t* dataspace_id, hid_t* datatype, hsize_t* dims, herr_t* status, int dimen, const char* name, double* address)
  {
    *dataspace_id = H5Screate_simple(dimen, dims, NULL);
    *dataset_id = H5Dcreate(*file_id, name, H5T_NATIVE_DOUBLE, *dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    *status = H5Dwrite(*dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, address);
    *status = H5Dclose(*dataset_id);
    *status = H5Sclose(*dataspace_id);
  }

  void writeDoubleDataToFile(hid_t* file_id, hid_t* dataset_id, hid_t* dataspace_id, hid_t* datatype, hsize_t* dims, herr_t* status, int dimen, const char* name, double** address)
  {
    *dataspace_id = H5Screate_simple(dimen, dims, NULL);
    *dataset_id = H5Dcreate(*file_id, name, H5T_NATIVE_DOUBLE, *dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    *status = H5Dwrite(*dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, address);
    *status = H5Dclose(*dataset_id);
    *status = H5Sclose(*dataspace_id);
  }

  void readIntDataFromFile(hid_t* file, hid_t* dataset, hid_t* filespace, hid_t* memspace,  hsize_t* dimz,const char* dataname, int*& data)
  {
      int rank;
      herr_t status,status_n;

      *dataset = H5Dopen(*file, dataname, H5P_DEFAULT);
      *filespace = H5Dget_space(*dataset);
      rank = H5Sget_simple_extent_ndims(*filespace);
      status_n  = H5Sget_simple_extent_dims(*filespace, dimz, NULL);
      *memspace = H5Screate_simple(rank,dimz,NULL);
      if(rank == 1)
          data = (int *) malloc(dimz[0]*sizeof(int));      
      else
          data = (int *) malloc(dimz[0]*dimz[1]*sizeof(int)); 

      status = H5Dread(*dataset, H5T_NATIVE_INT, *memspace, *filespace, H5P_DEFAULT, data);
  }

  void readFloatDataFromFile(hid_t* file, hid_t* dataset, hid_t* filespace, hid_t* memspace, hsize_t* dimz,const char* dataname, float*& data)
  {
      int rank;
      herr_t status,status_n;
      *dataset = H5Dopen(*file, dataname, H5P_DEFAULT);
      *filespace = H5Dget_space(*dataset);
      rank = H5Sget_simple_extent_ndims(*filespace);
      status_n  = H5Sget_simple_extent_dims(*filespace, dimz, NULL);
      *memspace = H5Screate_simple(rank,dimz,NULL);
      if(rank == 1)
          data = (float *) malloc(dimz[0]*sizeof(float));        
      else
          data = (float *) malloc(dimz[0]*dimz[1]*sizeof(float));      
      
      status = H5Dread(*dataset, H5T_NATIVE_FLOAT, *memspace, *filespace, H5P_DEFAULT, data);
  }

  void readDoubleDataFromFile(hid_t* file, hid_t* dataset, hid_t* filespace, hid_t* memspace,  hsize_t* dimz,const char* dataname, double*& data)
  {
      int rank;
      herr_t status,status_n;
      *dataset = H5Dopen(*file, dataname, H5P_DEFAULT);
      *filespace = H5Dget_space(*dataset);
      rank = H5Sget_simple_extent_ndims(*filespace);
      status_n  = H5Sget_simple_extent_dims(*filespace, dimz, NULL);
      *memspace = H5Screate_simple(rank,dimz,NULL);
      if(rank == 1)
          data = (double *) malloc(dimz[0]*sizeof(double));      
      else
          data = (double *) malloc(dimz[0]*dimz[1]*sizeof(double));
      status = H5Dread(*dataset, H5T_NATIVE_DOUBLE, *memspace, *filespace, H5P_DEFAULT, data);
  }

  inline void hdf_close(hid_t dataset, hid_t filespace, hid_t memspace)
  {
      H5Dclose(dataset);
      H5Sclose(filespace);
      H5Sclose(memspace);
  }

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
  /**create hdf5 file and xdmf file. hdf5 file store the basic data in psurface. xdmf is used to read the hdf5 file by 
  *paraview.
  */
  template<class ctype,int dim>
  bool psurface::PsurfaceConvert<ctype,dim>::creatHdfAndXdmf(const std::string&  xdf_filename, const std::string&  hdf_filename, bool base)
  {
    if(base)
    {
      writeHdf5Data(hdf_filename);
      writeXdmf(xdf_filename, hdf_filename);
    }
    else 
      writeBaseHdf5Data(hdf_filename);

     return 0;
  }

  template<class ctype,int dim>
  ///write the data array into hdf5 data structure
  bool psurface::PsurfaceConvert<ctype,dim>::writeHdf5Data(const std::string& filename)
  {
    hid_t     file_id;
    file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    int i,j;
    hid_t     dataset_id, dataspace_id,datatype;
    hsize_t   dims[] ={numVertices};
    hsize_t   dimz[2];
    herr_t    status;

    //1) 'Patches'
    int *patches_vec = (int *) malloc(patches.size()*3*sizeof(int));
    for(i = 0; i < patches.size();i++)
    {
        patches_vec[3*i] = (patches[i]).innerRegion;
        patches_vec[3*i + 1] = (patches[i]).outerRegion;
        patches_vec[3*i + 2] = (patches[i]).boundaryId;
    }
    dims[0] = patches.size()*3;
    writeIntDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dims, &status, 2, "/patches_vec", patches_vec);
    free(patches_vec);    

    //2) 'BaseGridVertexCoords'
    ctype **basecoords;
    basecoords = new float *[numVertices];
    for(i = 0; i < numVertices; i++)
      basecoords[i] = new float[3];
    for(i = 0; i < numVertices; i++)
    {
      basecoords[i][0] = (baseGridVertexCoordsArray[i])[0];
      basecoords[i][1] = (baseGridVertexCoordsArray[i])[1];
      basecoords[i][2] = (baseGridVertexCoordsArray[i])[2];
    }
    dimz[0] = numVertices;
    dimz[1] = 3;
    writeFloatDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dimz, &status, 2, "/BaseCoords", basecoords);
  
    //3) 'BaseGridTriangles'
    int **base_tri;
    base_tri = new int *[numTriangles];
    for(i = 0; i < numTriangles;i++) 
        base_tri[i] = new int[4];
    for(i = 0; i < numTriangles;i++)
    {
      base_tri[i][0] = 4; //topology number of triangle in xdmf
      base_tri[i][1] = (baseGridTriArray[i])[0];
      base_tri[i][2] = (baseGridTriArray[i])[1];
      base_tri[i][3] = (baseGridTriArray[i])[2];
    }
    dimz[0] = numTriangles;
    dimz[1] = 4;
    writeIntDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dimz, &status, 2, "/BaseTri", base_tri);
  
    //4) NodePositions(x, y, and z-coordinates of the image position).
    //ipos
    ctype **ipos;
    ipos = new float *[iPos.size()];
    for(i = 0; i < iPos.size();i++)
        ipos[i] = new float[3];
    for(i = 0; i < iPos.size(); i++)
    {
        ipos[i][0] = (iPos[i])[0];
        ipos[i][1] = (iPos[i])[1];
        ipos[i][2] = (iPos[i])[2];
    }
    dimz[0] = iPos.size();
    dimz[1] = 3;
    writeFloatDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dimz, &status, 2, "/iPos", ipos);
    free(ipos);

    //5) 'NumNodesAndParameterEdgesPerTriangle'
    int **num_nodes_and_edges_array;
    num_nodes_and_edges_array = new int *[numTriangles];
    for(i = 0; i < numTriangles; i++)
        num_nodes_and_edges_array[i] = new int[11];
    for(i = 0; i < numTriangles;i++)
        for(j = 0; j < 11; j++)
            num_nodes_and_edges_array[i][j] = numNodesAndEdgesArray[11*i + j];
    dimz[0] = numTriangles;
    dimz[1] = 11;
    writeIntDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dimz, &status, 2, "/numNodesAndEdgesArray", num_nodes_and_edges_array);
    free(num_nodes_and_edges_array);

    //6) 'Nodes'
    // barycentric coordinates on the respective triangle and 
    // x, y, and z-coordinates of the image position.
    //nodes on plane surface of triangle
    ctype **nodecoords;
    nodecoords = new float *[numNodes];
    for(i = 0; i< numNodes; i++)
            nodecoords[i]= new float[3];
    for(i = 0; i < numNodes; i++)
    {
      nodecoords[i][0] = (nodePositions[i])[0];
      nodecoords[i][1] = (nodePositions[i])[1];
      nodecoords[i][2] = (nodePositions[i])[2];
    }
    dimz[0] = numNodes;
    dimz[1] = 3;
    writeFloatDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dimz, &status, 2, "/NodeCoords", nodecoords);
    free(nodecoords);

    //NodeData 
    ctype **dp;
    dp = new float *[numNodes];
    for(i = 0; i < numNodes; i++)
        dp[i] = new float[2];
    for(i = 0; i < numNodes; i++)
    {
        dp[i][0] = (domainPositions[i])[0];
        dp[i][1] = (domainPositions[i])[1];
    }
    dimz[0] = numNodes;
    dimz[1] = 2;
    writeFloatDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dimz, &status, 2, "/LocalNodePosition", dp);
    free(dp);

    //image position
    ctype **imageposition;
    imageposition = new float *[nvertices];
    for(i = 0; i < nvertices; i++)
        imageposition[i] = new float[3];
    for (i = 0; i< nvertices; i++)
    {
      imageposition[i][0] = (imagePos[i])[0];
      imageposition[i][1] = (imagePos[i])[1];
      imageposition[i][2] = (imagePos[i])[2];
    }
    dimz[0] = nvertices;
    dimz[1] = 3;
    writeFloatDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dimz, &status, 2 , "/ImagePosition", imageposition);
    free(imageposition);

    //7)NodeNumbers
    int *nodenumber = (int*)malloc(numNodes*sizeof(int));
    for(i = 0; i < numNodes; i++)
        nodenumber[i] = nodeNumber[i];
    writeIntDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dims, &status, numNodes, "/NodeNumber", nodenumber);
    free(nodenumber);

    //8) 'ParameterEdges'
    //connection array(in global index) of parameter edges
    int **parameter_edge_local_array;
    parameter_edge_local_array = new int *[numParamEdges];
    for(i = 0; i < numParamEdges; i++)
        parameter_edge_local_array[i] = new int[2];
    for(i = 0; i < numParamEdges; i++)
    {
       parameter_edge_local_array[i][0] =  (parameterEdgeArrayLocal[i])[0];
       parameter_edge_local_array[i][1] = (parameterEdgeArrayLocal[i])[1];
    }
    dimz[0] = numParamEdges;
    dimz[1] = 2;
    writeIntDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dimz, &status, 2, "/LocalParamEdge", parameter_edge_local_array);    
    free( parameter_edge_local_array);

    //param edge
    int  **paramedge;
    paramedge = new int *[numParamEdges];
    for(i = 0; i < numParamEdges; i++)
        paramedge[i] = new int[4];
    for(i = 0; i < numParamEdges; i++)
    {
       paramedge[i][0] = 2;
       paramedge[i][1] = 2;
       paramedge[i][2] = (parameterEdgeArray[i])[0];
       paramedge[i][3] = (parameterEdgeArray[i])[1];
    }
    dimz[0] = numParamEdges;
    dimz[1] = 4;
    writeIntDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dimz, &status, 2, "/ParamEdge", paramedge);
    free(paramedge);

    //9) 'EdgePoints'
    int edgepointsarray[edgePointsArray.size()];
    for(i = 0; i < edgePointsArray.size();i++)
      edgepointsarray[i] = edgePointsArray[i];
    dims[0] = edgePointsArray.size();
    writeIntDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dims, &status, 1, "/EdgePoints", edgepointsarray);

    //Supportive data
    //params
    int psurfaceparams[4];
    psurfaceparams[0] = numVertices;
    psurfaceparams[1] = numNodes;
    psurfaceparams[2] = numTriangles;
    psurfaceparams[3] = numParamEdges;

    dims[0] = 4;
    writeIntDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dims, &status, 1, "/Params", psurfaceparams);

    //nodetype
    int nodetype[nvertices];
    for (i = 0; i < nvertices; i++)
      nodetype[i] = nodeType[i];

    dims[0] = nvertices;
    writeIntDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dims, &status, 1, "/NodeType", nodetype);
    
    //close the file
    status = H5Fclose(file_id);
    return 0;
  };

  template<class ctype,int dim>
  ///write the data array into hdf5 data structure(This function store exactly the same data as amiramesh does) 
  bool psurface::PsurfaceConvert<ctype,dim>::writeBaseHdf5Data(const std::string& filename)
  {
    hid_t     file_id;
    file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    int i,j;
    hid_t     dataset_id, dataspace_id,datatype;
    hsize_t   dims[] ={numVertices};
    hsize_t   dimz[2];
    herr_t    status;
    
    //1) 'Patches'
    int *patches_vec = (int *) malloc(patches.size()*3*sizeof(int));
    for(i = 0; i < patches.size();i++)
    {
        patches_vec[3*i] = (patches[i]).innerRegion;
        patches_vec[3*i + 1] = (patches[i]).outerRegion;
        patches_vec[3*i + 2] = (patches[i]).boundaryId;
    }
    dims[0] = patches.size()*3;
    writeIntDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dims, &status, 2, "/patches_vec", patches_vec);
    free(patches_vec);    

    //2) 'BaseGridVertexCoords'
    ctype **basecoords;
    basecoords = new float *[numVertices];
    for(i = 0; i < numVertices; i++)
      basecoords[i] = new float[3];
    for(i = 0; i < numVertices; i++)
    {
      basecoords[i][0] = (baseGridVertexCoordsArray[i])[0];
      basecoords[i][1] = (baseGridVertexCoordsArray[i])[1];
      basecoords[i][2] = (baseGridVertexCoordsArray[i])[2];
    }
    dimz[0] = numVertices;
    dimz[1] = 3;
    writeFloatDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dimz, &status, 2, "/BaseCoords", basecoords);
  
    //3) 'BaseGridTriangles'
    int **base_tri;
    base_tri = new int *[numTriangles];
    for(i = 0; i < numTriangles;i++) 
        base_tri[i] = new int[4];
    for(i = 0; i < numTriangles;i++)
    {
      base_tri[i][0] = 4; //topology number of triangle in xdmf
      base_tri[i][1] = (baseGridTriArray[i])[0];
      base_tri[i][2] = (baseGridTriArray[i])[1];
      base_tri[i][3] = (baseGridTriArray[i])[2];
    }
    dimz[0] = numTriangles;
    dimz[1] = 4;
    writeIntDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dimz, &status, 2, "/BaseTri", base_tri);
  
    //4) NodePositions(x, y, and z-coordinates of the image position).
    //ipos
    ctype **ipos;
    ipos = new float *[iPos.size()];
    for(i = 0; i < iPos.size();i++)
        ipos[i] = new float[3];
    for(i = 0; i < iPos.size(); i++)
    {
        ipos[i][0] = (iPos[i])[0];
        ipos[i][1] = (iPos[i])[1];
        ipos[i][2] = (iPos[i])[2];
    }
    dimz[0] = iPos.size();
    dimz[1] = 3;
    writeFloatDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dimz, &status, 2, "/iPos", ipos);
    free(ipos);

    //5) 'NumNodesAndParameterEdgesPerTriangle'
    int **num_nodes_and_edges_array;
    num_nodes_and_edges_array = new int *[numTriangles];
    for(i = 0; i < numTriangles; i++)
        num_nodes_and_edges_array[i] = new int[11];
    for(i = 0; i < numTriangles;i++)
        for(j = 0; j < 11; j++)
            num_nodes_and_edges_array[i][j] = numNodesAndEdgesArray[11*i + j];
    dimz[0] = numTriangles;
    dimz[1] = 11;
    writeIntDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dimz, &status, 2, "/numNodesAndEdgesArray", num_nodes_and_edges_array);
    free(num_nodes_and_edges_array); 

    //6) 'Nodes'
    //NodeData 
    ctype **dp;
    dp = new float *[numNodes];
    for(i = 0; i < numNodes; i++)
        dp[i] = new float[2];
    for(i = 0; i < numNodes; i++)
    {
        dp[i][0] = (domainPositions[i])[0];
        dp[i][1] = (domainPositions[i])[1];
    }
    dimz[0] = numNodes;
    dimz[1] = 2;
    writeFloatDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dimz, &status, 2, "/LocalNodePosition", dp);
    free(dp);

    //7)NodeNumbers
    int *nodenumber = (int*)malloc(numNodes*sizeof(int));
    for(i = 0; i < numNodes; i++)
        nodenumber[i] = nodeNumber[i];
    writeIntDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dims, &status, numNodes, "/NodeNumber", nodenumber);
    free(nodenumber);
    
    //8) 'ParameterEdges'
    //connection array(in global index) of parameter edges
    int **parameter_edge_local_array;
    parameter_edge_local_array = new int *[numParamEdges];
    for(i = 0; i < numParamEdges; i++)
        parameter_edge_local_array[i] = new int[2];
    for(i = 0; i < numParamEdges; i++)
    {
       parameter_edge_local_array[i][0] =  (parameterEdgeArrayLocal[i])[0];
       parameter_edge_local_array[i][1] = (parameterEdgeArrayLocal[i])[1];
    }
    dimz[0] = numParamEdges;
    dimz[1] = 2;
    writeIntDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dimz, &status, 2, "/LocalParamEdge", parameter_edge_local_array);    
    free( parameter_edge_local_array);

    //9) 'EdgePoints'
    int *edgepointsarray = (int*)malloc(edgePointsArray.size()*sizeof(int));
    for(i = 0; i < edgePointsArray.size();i++)
      edgepointsarray[i] = edgePointsArray[i];
    writeIntDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dims, &status, edgePointsArray.size(), "/EdgePoints", edgepointsarray);
    free(edgepointsarray);

    //Supportive data
    //params
    int psurfaceparams[4];
    psurfaceparams[0] = numVertices;
    psurfaceparams[1] = numNodes;
    psurfaceparams[2] = numTriangles;
    psurfaceparams[3] = numParamEdges;

    dims[0] = 4;
    writeIntDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dims, &status, 1, "/Params", psurfaceparams);

    //close the file
    status = H5Fclose(file_id);
    return 0;
  };

  template<class ctype,int dim>
  ///writhe the xdmf file which store the structure information of hdf5 file.
  bool psurface::PsurfaceConvert<ctype,dim>::writeXdmf(const std::string&  xdf_filename, const std::string&  hdf_filename)
  {
    printf("in writeXdmd function: xdf_file = %s hdf_file = %s\n", xdf_filename.c_str(),hdf_filename.c_str());
    FILE *xmf = 0;
    xmf = fopen(xdf_filename.c_str(), "w");
    fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
    fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
    fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
    fprintf(xmf, "<Domain>\n");

    fprintf(xmf, "<Grid Name=\"basegrid\" GridType=\"Uniform\">\n");
    fprintf(xmf, "<Topology TopologyType=\"Mixed\" NumberOfElements=\"%d\">\n", numTriangles);
    fprintf(xmf, "<DataItem Dimensions = \"%d\" NumberType=\"int\" Format=\"HDF\">\n",  4*numTriangles);
    fprintf(xmf, "%s:/BaseTri\n", hdf_filename.c_str());
    fprintf(xmf, "</DataItem>\n");
    fprintf(xmf, "</Topology>\n");

    fprintf(xmf, "<Geometry GeometryType=\"XYZ\">\n");
    fprintf(xmf, "<DataItem Name = \"triangles\" Dimensions=\"%d 3\" NumberType=\"Float\" Format=\"HDF\">\n", numVertices);
    fprintf(xmf, "%s:/BaseCoords\n", hdf_filename.c_str());
    fprintf(xmf, "</DataItem>\n");
    fprintf(xmf, "</Geometry>\n");

    fprintf(xmf, "<Attribute Name=\"Patches\" AttributeType=\"Scalar\" Center=\"Cell\">\n");
      fprintf(xmf, "<DataItem ItemType=\"HyperSlab\" Dimensions=\"%d\" Type=\"HyperSlab\">\n",numTriangles);

        fprintf(xmf, "<DataItem  Dimensions=\"3\" Format=\"XML\">\n");
        fprintf(xmf, "%d\n", 4);
        fprintf(xmf, "%d\n", 11);
        fprintf(xmf, "%d\n", numTriangles);
        fprintf(xmf, "</DataItem>\n"); 

        fprintf(xmf, " <DataItem Dimensions=\"%d\" NumberType=\"int\" Format=\"HDF\">\n", 11*numTriangles);
        fprintf(xmf, "%s:/numNodesAndEdgesArray\n", hdf_filename.c_str());
        fprintf(xmf, "</DataItem>\n");

      fprintf(xmf, "</DataItem>\n");
    
    fprintf(xmf, "</Attribute>\n");

    fprintf(xmf, "</Grid>\n");

    fprintf(xmf, "<Grid Name=\"paramedge\" GridType=\"Uniform\">\n");
    fprintf(xmf, "<Topology TopologyType=\"Mixed\" NumberOfElements=\"%d\">\n", numParamEdges);
          fprintf(xmf, "<DataItem Dimensions = \"%d\" NumberType=\"int\" Format=\"HDF\">\n",  4*numParamEdges);
            fprintf(xmf, "%s:/ParamEdge\n", hdf_filename.c_str());
          fprintf(xmf, "</DataItem>\n");
    fprintf(xmf, "</Topology>\n");

    fprintf(xmf, "<Geometry GeometryType=\"XYZ\">\n");
    fprintf(xmf, "<DataItem ItemType= \"Function\"  Dimensions=\"%d 3\" Function=\" JOIN($0;$1) \">\n", nvertices);
    fprintf(xmf, "<DataItem Name = \"triangles\" Dimensions=\"%d 3\" NumberType=\"Float\" Format=\"HDF\">\n", numVertices);
    fprintf(xmf, "%s:/BaseCoords\n", hdf_filename.c_str());
    fprintf(xmf, "</DataItem>\n");
    fprintf(xmf, "<DataItem Name = \"nodecoords\" Dimensions=\"%d 3\" NumberType=\"Float\" Format=\"HDF\">\n", numNodes);
    fprintf(xmf, "%s:/NodeCoords\n", hdf_filename.c_str());
    fprintf(xmf, "</DataItem>\n");
    fprintf(xmf, "</DataItem>\n");
    fprintf(xmf, "</Geometry>\n");

      fprintf(xmf, "<Attribute Name=\"Nodetype\" AttributeType=\"Scalar\" Center=\"Node\">\n");
        fprintf(xmf, "<DataItem Dimensions=\"%d\" NumberType=\"int\" Format=\"HDF\">\n", nvertices);
          fprintf(xmf, "%s:/NodeType\n", hdf_filename.c_str());
        fprintf(xmf, "</DataItem>\n");
      fprintf(xmf, "</Attribute>\n");


      fprintf(xmf, "<Attribute Name=\"imageposition\" AttributeType=\"Vector\" Center=\"Node\">\n");
        fprintf(xmf, "<DataItem Dimensions=\"%d %d\" NumberType=\"ctype\" Format=\"HDF\">\n", nvertices, 3);
          fprintf(xmf, "%s:/ImagePosition\n", hdf_filename.c_str());
        fprintf(xmf, "</DataItem>\n");
      fprintf(xmf, "</Attribute>\n");
    fprintf(xmf, "</Grid>\n");

    fprintf(xmf, "</Domain>\n");
    fprintf(xmf, "</Xdmf>\n");
    fclose(xmf);
    return true;
  };
 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///write the psurface into vtu file
  template<class ctype,int dim>
  bool psurface::PsurfaceConvert<ctype,dim>::creatVTU(const char *filename, bool basegrid)
  {
    std::ofstream file;
    file.open(filename);
    if (! file.is_open()) printf("%c does not exits!\n", filename);
    writeDataFile(file, basegrid);
    file.close();
    return 0;
  }
  ///write data file to stream
  template<class ctype,int dim>
  void psurface::PsurfaceConvert<ctype,dim>::writeDataFile(std::ostream& s, bool basegrid)
  {
    VTK::FileType fileType = VTK::unstructuredGrid;

    VTK::VTUWriter writer(s, outputtype,fileType);//Most inportant structure used here
    
    if(basegrid)
      writer.beginMain(numTriangles, numVertices);
    else
      writer.beginMain(numTriangles + numParamEdges, numVertices + numNodes);
    
    writeAllData(writer, basegrid);
    writer.endMain();
  }

  ///write the data section in vtu
  template<class ctype,int dim>
  void psurface::PsurfaceConvert<ctype,dim>::writeAllData(VTK::VTUWriter& writer, bool basegrid) {
    //PointData
    writePointData(writer,basegrid);
    // Points
    writeGridPoints(writer,basegrid);
    // Cells
    writeGridCells(writer,basegrid);
  }

  
  //! write point data
  template<class ctype,int dim>
  void psurface::PsurfaceConvert<ctype,dim>::writePointData(VTK::VTUWriter& writer, bool basegrid)
  {
    std::string scalars = "nodetype";
    std::string vectors = "imageposition";
    std::vector<int>::iterator pt;
    writer.beginPointData(scalars, vectors);
    {
      if(basegrid)
      {
        std::tr1::shared_ptr<VTK::DataArrayWriter<ctype> > p
        (writer.makeArrayWriter<ctype>(scalars, 1, numVertices));
        for( int i = 0; i < numVertices;i++)
          p->write(nodeType[i]);
      }
      else 
      {
        std::tr1::shared_ptr<VTK::DataArrayWriter<ctype> > p
        (writer.makeArrayWriter<ctype>(scalars, 1, nvertices));
        for( int i = 0; i < nodeType.size();i++)
          p->write(nodeType[i]);
      }
    }
    v_iterator pi;
    {
      if( !basegrid)
      {
        std::tr1::shared_ptr<VTK::DataArrayWriter<ctype> > p
        (writer.makeArrayWriter<ctype>(vectors, 3, nvertices));
        for(int i = 0; i < imagePos.size(); i++)
        {
          for(int l = 0; l < 3; l++) p->write((imagePos[i])[l]);
        }
      }
    }
    writer.endPointData();
  }


  //! write the positions of vertices
  template<class ctype,int dim>
  void psurface::PsurfaceConvert<ctype,dim>::writeGridPoints(VTK::VTUWriter& writer, bool basegrid)
  {
    writer.beginPoints();
    v_iterator vp;
    v_iterator dp;
    {
      if(basegrid)
      {
        std::tr1::shared_ptr<VTK::DataArrayWriter<ctype> > p
        (writer.makeArrayWriter<ctype>("Coordinates", 3, numVertices));
        if(!p->writeIsNoop()) {
        for(int i = 0; i < baseGridVertexCoordsArray.size(); i++)
          for(int l = 0; l < 3; l++) p->write((baseGridVertexCoordsArray[i])[l]);
        }
      }
      else
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
    }
    //	p.reset();
    writer.endPoints();
  }

  //! write the connectivity array
  template<class ctype,int dim>
  void psurface::PsurfaceConvert<ctype,dim>::writeGridCells(VTK::VTUWriter& writer, bool basegrid)
  {
    t_iterator tp;
    e_iterator ep;
    writer.beginCells();
    // connectivity
    {
      if(basegrid)
      {
        std::tr1::shared_ptr<VTK::DataArrayWriter<int> > p1
        (writer.makeArrayWriter<int>("connectivity", 1, 3*numTriangles));
        if(!p1->writeIsNoop())
        {
          for(int i = 0; i < baseGridTriArray.size(); i++)
            for( int l = 0; l < 3; l++)
              p1->write((baseGridTriArray[i])[l]);
        }
      }
      else
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
    }
    // offsets
    {
      if(basegrid)
      {
        std::tr1::shared_ptr<VTK::DataArrayWriter<int> > p2
        (writer.makeArrayWriter<int>("offsets", 1, numTriangles));
        if(!p2->writeIsNoop()) {
        int offset = 0;
        for(int i = 0; i < baseGridTriArray.size(); i++)
        {
          offset += 3;
          p2->write(offset);
        }
        }
      }
      else
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
    }
    // types
    {
      int offset = 0;
      if(basegrid)
      {
        std::tr1::shared_ptr<VTK::DataArrayWriter<unsigned char> > p3
        (writer.makeArrayWriter<unsigned char>("types", 1, numVertices));
        if(!p3->writeIsNoop())
        {
          for(int i = 0; i < baseGridTriArray.size();i++)
          p3->write(5); //vtktype of triangle
        }
      }
      else
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
  
  ///initialize PsurfaceConvert from the psurface object
  template<class ctype,int dim>
  psurface::PsurfaceConvert<ctype,dim>::PsurfaceConvert(PSurface<dim,ctype>* psurface)
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
    std::vector<StaticVector<ctype,3> > imagearray;
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

    //nodes image positions(only for intersection points)
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
    StaticVector<ctype,3> imagepos;

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

  ///initialize PsurfaceConvert from file
  ///read data array from hdf5 or gmsh file into the data arrary in class PsurfaceConvert
    template<class ctype,int dim>
  psurface::PsurfaceConvert<ctype,dim>::PsurfaceConvert(const std::string&  filename, bool a)
  {
      if(a) //initialize from hdf5 data
          readHdf5Data(filename);
      else
          readGmsh(filename);
  }   

  template<class ctype,int dim>
  bool psurface::PsurfaceConvert<ctype,dim>::readHdf5Data(const std::string&  filename)
  {
      hid_t file;
      hid_t datatype, dataset;
      hid_t filespace;
      hid_t       memspace;
      hsize_t dims[2];
      hsize_t dimz[1];
      herr_t status,status_n;
      int rank;
      int i,j,k;
      file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

      //read params
      int *psurfaceparams;
      readIntDataFromFile(&file, &dataset, &filespace, &memspace, dimz,"/Params", psurfaceparams);
      numVertices = psurfaceparams[0];
      numNodes = psurfaceparams[1];
      numTriangles = psurfaceparams[2];
      numParamEdges = psurfaceparams[3];
      nvertices = numVertices + numNodes;
      ncells = numTriangles + numParamEdges;
      hdf_close(dataset, filespace, memspace);
      //read xyz coordinate of vertices
      {     baseGridVertexCoordsArray.resize(3*numVertices);
      ctype *coords;  
      readFloatDataFromFile(&file, &dataset, &filespace, &memspace, dimz, "/BaseCoords", coords);
      hdf_close(dataset, filespace, memspace);
      for(j = 0; j < numVertices ;j++) 
          for(i = 0; i < 3; i++)
              (baseGridVertexCoordsArray[j])[i] = coords[3*j + i];
      }
/*      {
      hsize_t     dims_out[2];
      dataset = H5Dopen(file, "/BaseCoords", H5P_DEFAULT);
      datatype  = H5Dget_type(dataset);
      filespace = H5Dget_space(dataset);
      status_n  = H5Sget_simple_extent_dims(filespace, dims_out, NULL);      
      memspace = H5Screate_simple(2,dims_out,NULL);
      float coords[dims_out[0]][dims_out[1]];
      status = H5Dread(dataset, H5T_NATIVE_FLOAT, memspace, filespace, H5P_DEFAULT,coords);
        for(j = 0; j < numVertices ;j++) 
          for(i = 0; i < 3; i++)
              (baseGridVertexCoordsArray[j])[i] = coords[j][i];
      }
*/
      //triangle
      int *tri;
      readIntDataFromFile(&file, &dataset, &filespace, &memspace, dimz,"/BaseTri", tri);
      hdf_close(dataset, filespace, memspace);
        
      baseGridTriArray.resize(numTriangles);
      for(j = 0; j < numTriangles;j++)
        for(i = 0; i < 3; i++)
          (baseGridTriArray[j])[i] = tri[4*j+ 1 + i];
      //Parameter Edge
      int *lparam;
      readIntDataFromFile(&file, &dataset, &filespace, &memspace, dimz,"/LocalParamEdge", lparam);
      hdf_close(dataset, filespace, memspace);
      parameterEdgeArrayLocal.resize(numParamEdges);
      for(j = 0; j < dimz[0]/2;j++)
         for(i = 0; i < 2; i++)
             (parameterEdgeArrayLocal[j])[i] = lparam[2*j + i];

      //nodeNumber
      int *nodenumber;
      readIntDataFromFile(&file, &dataset, &filespace, &memspace, dimz,"/NodeNumber", nodenumber);
      hdf_close(dataset, filespace, memspace);
      nodeNumber.resize(numNodes);
      for(i = 0; i < numNodes;i++) nodeNumber[i] = nodenumber[i];

      //iPos
      ctype *coords;
      readFloatDataFromFile(&file, &dataset, &filespace, &memspace, dimz,"/iPos",coords);
      hdf_close(dataset, filespace, memspace);
      iPos.resize(dimz[0]/3);
      for(j = 0; j < dimz[0]/3; j++)
      {
        (iPos[j])[0] = coords[3*j];
        (iPos[j])[1] = coords[3*j+1];
        (iPos[j])[2] = coords[3*j+2];
      }

      //nodeNodesAndEdgesArray
      int *nodearray;
      readIntDataFromFile(&file, &dataset, &filespace, &memspace, dimz, "/numNodesAndEdgesArray",nodearray);
      hdf_close(dataset, filespace, memspace);
      numNodesAndEdgesArray.resize(11*numTriangles);
      for(i = 0; i < 11*numTriangles;i++)
      numNodesAndEdgesArray[i] = nodearray[i];

      //local position of nodes on triangle
      ctype *localpos;
      readFloatDataFromFile(&file, &dataset, &filespace, &memspace, dimz, "/LocalNodePosition",localpos);
      hdf_close(dataset, filespace, memspace);
      domainPositions.resize(numNodes);
      for(i = 0; i < numNodes;i++)
      for(j = 0; j < 2;j++)
          (domainPositions[i])[j] = localpos[2*i + j];

      //edgepointsarray
      int *edgep;
      readIntDataFromFile(&file, &dataset, &filespace, &memspace, dimz, "/EdgePoints", edgep);
      hdf_close(dataset, filespace, memspace);
      edgePointsArray.resize(dimz[0]);
      for(i = 0; i < dimz[0];i++)
         edgePointsArray[i] = edgep[i];

      //patches
      int *patch;
      readIntDataFromFile(&file, &dataset, &filespace, &memspace, dimz, "/patches_vec", patch);
      hdf_close(dataset, filespace, memspace);
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

  ///read psurface_convert from Gmsh file
  template<class ctype,int dim>
  bool psurface::PsurfaceConvert<ctype,dim>::readGmsh(const std::string&  filename)
  {
      FILE* file = fopen(filename.c_str(),"r");
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

        StaticVector<int, 3> elementDofs(3);

        readfile(file, 3, formatString.c_str(), &(elementDofs[0]),&(elementDofs[1]),&(elementDofs[2]));
        triArray.push_back(elementDofs);
      }

      //remove vetices which is not the corner of triangle
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

  ///creat surface from psurface_convert
  ///baseTriangle == 1  we create the psurface subject which only have base grid triangles
  template<class ctype,int dim>
  bool psurface::PsurfaceConvert<ctype,dim>::initPsurface(PSurface<2,ctype>* psurf, Surface* surf, bool baseTriangle)
  {
      if(baseTriangle)
          return initBaseCasePSurface(psurf, surf);
      else
          return initCompletePSurface(psurf, surf);
  }

  ///Create PSurface which only have BaseGridVertex and BaseTriangle
  template<class ctype,int dim>
  bool psurface::PsurfaceConvert<ctype,dim>::initBaseCasePSurface(PSurface<2,ctype>* psurf, Surface* surf)
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
    PSurface<2,float>::Patch Patch; //= new PSurface<2, ctype>::Patch;
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
  template<class ctype,int dim>
  bool psurface::PsurfaceConvert<ctype,dim>::initCompletePSurface(PSurface<2,ctype>* psurf, Surface* surf)
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
    PSurface<2,float>::Patch Patch; //= new PSurface<2, ctype>::Patch;
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

template class PsurfaceConvert<float,2>;
//template class PsurfaceConvert<double,2>;




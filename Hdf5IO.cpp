#include <vector>
#include <string.h>
#include <hdf5.h>
#include <memory>
#include <tr1/memory>

#include "StaticVector.h"
#include "Domains.h"
#include "PSurface.h"
#include "PSurfaceFactory.h"
#include "Hdf5IO.h"
#include "psurface_convert.h"
using namespace psurface;
void writeIntDataToFile(hid_t* file_id, hid_t* dataset_id, hid_t* dataspace_id, hid_t* datatype, hsize_t* dims, herr_t* status, const char*  name, int address[])
{
    *dataspace_id = H5Screate_simple(1, dims, NULL);
    *dataset_id = H5Dcreate(*file_id, name, H5T_NATIVE_INT, *dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    *status = H5Dwrite(*dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, address);
    *status = H5Dclose(*dataset_id);
    *status = H5Sclose(*dataspace_id);
}

void writeIntDataToFile(hid_t* file_id, hid_t* dataset_id, hid_t* dataspace_id, hid_t* datatype, hsize_t* dims, herr_t* status, const char*  name, int address[][2])
{
    *dataspace_id = H5Screate_simple(2, dims, NULL);
    *dataset_id = H5Dcreate(*file_id, name, H5T_NATIVE_INT, *dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    *status = H5Dwrite(*dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, address);
    *status = H5Dclose(*dataset_id);
    *status = H5Sclose(*dataspace_id);
}

void writeIntDataToFile(hid_t* file_id, hid_t* dataset_id, hid_t* dataspace_id, hid_t* datatype, hsize_t* dims, herr_t* status, const char*  name, int address[][4])
{
    *dataspace_id = H5Screate_simple(2, dims, NULL);
    *dataset_id = H5Dcreate(*file_id, name, H5T_NATIVE_INT, *dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    *status = H5Dwrite(*dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, address);
    *status = H5Dclose(*dataset_id);
    *status = H5Sclose(*dataspace_id);
}

void writeIntDataToFile(hid_t* file_id, hid_t* dataset_id, hid_t* dataspace_id, hid_t* datatype, hsize_t* dims, herr_t* status, const char*  name, int address[][11])
{
    *dataspace_id = H5Screate_simple(2, dims, NULL);
    *dataset_id = H5Dcreate(*file_id, name, H5T_NATIVE_INT, *dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    *status = H5Dwrite(*dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, address);
    *status = H5Dclose(*dataset_id);
    *status = H5Sclose(*dataspace_id);
}

void writeFloatDataToFile(hid_t* file_id, hid_t* dataset_id, hid_t* dataspace_id, hid_t* datatype, hsize_t* dims, herr_t* status, const char* name, float address[])
  {
    *dataspace_id = H5Screate_simple(1, dims, NULL);
    *dataset_id = H5Dcreate(*file_id, name, H5T_NATIVE_FLOAT, *dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    *status = H5Dwrite(*dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, address);
    *status = H5Dclose(*dataset_id);
    *status = H5Sclose(*dataspace_id);
  }

void writeFloatDataToFile(hid_t* file_id, hid_t* dataset_id, hid_t* dataspace_id, hid_t* datatype, hsize_t* dims, herr_t* status, const char* name, float address[][2])
  {
    *dataspace_id = H5Screate_simple(2, dims, NULL);
    *dataspace_id = H5Screate_simple(2, dims, NULL);
    *dataset_id = H5Dcreate(*file_id, name, H5T_NATIVE_FLOAT, *dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    *status = H5Dwrite(*dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, address);
    *status = H5Dclose(*dataset_id);
  }

void writeFloatDataToFile(hid_t* file_id, hid_t* dataset_id, hid_t* dataspace_id, hid_t* datatype, hsize_t* dims, herr_t* status, const char* name, float address[][3])
  {
    *dataspace_id = H5Screate_simple(2, dims, NULL);
    *dataspace_id = H5Screate_simple(2, dims, NULL);
    *dataset_id = H5Dcreate(*file_id, name, H5T_NATIVE_FLOAT, *dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    *status = H5Dwrite(*dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, address);
    *status = H5Dclose(*dataset_id);
  }


  void writeDoubletDataToFile(hid_t* file_id, hid_t* dataset_id, hid_t* dataspace_id, hid_t* datatype, hsize_t* dims, herr_t* status, const char* name, double address[])
  {
    *dataspace_id = H5Screate_simple(1, dims, NULL);
    *dataset_id = H5Dcreate(*file_id, name, H5T_NATIVE_DOUBLE, *dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    *status = H5Dwrite(*dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, address);
    *status = H5Dclose(*dataset_id);
    *status = H5Sclose(*dataspace_id);
  }

  void writeDoubleDataToFile(hid_t* file_id, hid_t* dataset_id, hid_t* dataspace_id, hid_t* datatype, hsize_t* dims, herr_t* status, const char* name, double address[][2])
  {
    *dataspace_id = H5Screate_simple(2, dims, NULL);
    *dataset_id = H5Dcreate(*file_id, name, H5T_NATIVE_DOUBLE, *dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    *status = H5Dwrite(*dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, address);
    *status = H5Dclose(*dataset_id);
    *status = H5Sclose(*dataspace_id);
  }

  void writeDoubleDataToFile(hid_t* file_id, hid_t* dataset_id, hid_t* dataspace_id, hid_t* datatype, hsize_t* dims, herr_t* status, const char* name, double address[][3])
  {
    *dataspace_id = H5Screate_simple(2, dims, NULL);
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

  template<class ctype,int dim>
  void psurface::Hdf5IO<ctype,dim>::creatHdfAndXdmf(const std::string&  xdf_filename, const std::string&  hdf_filename, bool readablehdf)
  {
    if(readablehdf)
    {
      writeHdf5Data(hdf_filename);
      writeXdmf(xdf_filename, hdf_filename);
    }
    else 
      writeBaseHdf5Data(hdf_filename);

     return;
  }

  template<class ctype,int dim>
  ///write the data array into hdf5 data structure
  void psurface::Hdf5IO<ctype,dim>::writeHdf5Data(const std::string& filename)
  {
    hid_t     file_id;
    file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    hid_t     dataset_id, dataspace_id,datatype;
    hsize_t   dims[] ={numVertices};
    hsize_t   dimz[2];
    herr_t    status;
    
    int i, j, k;
    numVertices  = par->getNumVertices();
    numTriangles = par->getNumTriangles();

    //1) 'Patches'
    int patches_vec[par->patches.size()];
    for(i = 0; i < par->patches.size();i++)
    {
        patches_vec[3*i] = (par->patches[i]).innerRegion;
        patches_vec[3*i + 1] = (par->patches[i]).outerRegion;
        patches_vec[3*i + 2] = (par->patches[i]).boundaryId;
    }
    dims[0] = par->patches.size()*3;
    writeIntDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dims, &status, "/patches_vec", patches_vec);

    //2) 'BaseGridVertexCoords'
    ctype basecoords[numVertices][3];
    for(i = 0; i < numVertices; i++)
    {
      basecoords[i][0] = (par->vertices(i))[0];
      basecoords[i][1] = (par->vertices(i))[1];
      basecoords[i][2] = (par->vertices(i))[2];
    }
    dimz[0] = numVertices;
    dimz[1] = 3;
    writeFloatDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dimz, &status, "/BaseCoords", basecoords);
  
    //3) 'BaseGridTriangles'
    int base_tri[numTriangles][4];
    for(i = 0; i < numTriangles;i++)
    {
      base_tri[i][0] = 4; //topology number of triangle in xdmf
      base_tri[i][1] = par->triangles(i).vertices[0];
      base_tri[i][2] = par->triangles(i).vertices[1];
      base_tri[i][3] = par->triangles(i).vertices[2];
    }
    dimz[0] = numTriangles;
    dimz[1] = 4;
    writeIntDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dimz, &status, "/BaseTri", base_tri);
  
    //4) NodePositions(x, y, and z-coordinates of the image position).
    //ipos
    ctype ipos[par->iPos.size()][3];
    for(i = 0; i < par->iPos.size(); i++)
    {
        ipos[i][0] = (par->iPos[i])[0];
        ipos[i][1] = (par->iPos[i])[1];
        ipos[i][2] = (par->iPos[i])[2];
    }
    dimz[0] = par->iPos.size();
    dimz[1] = 3;
    writeFloatDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dimz, &status, "/iPos", ipos);

    //5) 'NumNodesAndParameterEdgesPerTriangle'
    int num_nodes_and_edges_array[numTriangles][11];
    numNodes      = 0;
    numParamEdges = 0;
    numEdgePoints = 0;
    for(i = 0; i < numTriangles;i++)
    {
        const DomainTriangle<ctype>& cT = par->triangles(i);

        int numIntersectionNodes;
        int numTouchingNodes;
        int numInteriorNodes;

        cT.countNodes(numIntersectionNodes, numTouchingNodes, numInteriorNodes);
        //        int numEdges = cT.getNumRegularEdges() - 3;
        int numEdges = cT.getNumRegularEdges();
        num_nodes_and_edges_array[i][0] = numIntersectionNodes;
        num_nodes_and_edges_array[i][1] = numTouchingNodes;
        num_nodes_and_edges_array[i][2] = numInteriorNodes;
        num_nodes_and_edges_array[i][3] = numEdges;
        num_nodes_and_edges_array[i][4] = cT.patch;

        num_nodes_and_edges_array[i][5] = cT.edgePoints[0].size()-2;
        num_nodes_and_edges_array[i][6] = cT.edgePoints[1].size()-2;
        num_nodes_and_edges_array[i][7] = cT.edgePoints[2].size()-2;

        num_nodes_and_edges_array[i][8] = cT.nodes[cT.cornerNode(0)].getNodeNumber();
        num_nodes_and_edges_array[i][9] = cT.nodes[cT.cornerNode(1)].getNodeNumber();
        num_nodes_and_edges_array[i][10] = cT.nodes[cT.cornerNode(2)].getNodeNumber();

        numNodes += numIntersectionNodes;
        numNodes += numTouchingNodes;
        numNodes += numInteriorNodes;

        numEdgePoints += cT.edgePoints[0].size() + cT.edgePoints[1].size() + cT.edgePoints[2].size() - 6;

        numParamEdges += numEdges;
    }

    dimz[0] = numTriangles;
    dimz[1] = 11;
    writeIntDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dimz, &status, "/numNodesAndEdgesArray", num_nodes_and_edges_array);

    /////////////////////////////////////////////////////////////////////    
    ncells = numTriangles + numParamEdges;
    nvertices = numVertices + numNodes;
    
    
    ctype nodePositions[numNodes][3];    
    ctype domainPositions[numNodes][2];
    ctype imagePos[nvertices][3];
    int   nodeNumber[numNodes];
    int   nodeType[numNodes];
    int   parameterEdgeArrayLocal[numParamEdges][2];
    int   parameterEdgeArray[numParamEdges][4];
    int   edgePointsArray[nvertices];

    for(i = 0; i < numVertices; i++) nodeType[i] = CORNER_NODE;

    int arrayIdx           = 0;
    int edgeArrayIdx       = 0;
    int edgePointsArrayIdx = 0;
    ctype cCoords[3][3];
    
    for (i=0; i<numTriangles; i++) {
        const DomainTriangle<ctype>& cT = par->triangles(i);

        std::vector<int> newIdx(cT.nodes.size());
        std::vector<int> newIdxlocal(cT.nodes.size());
        for(j = 0; j < 3; j++)
          for(k = 0; k < 3; k++)
            cCoords[j][k] = par->vertices(par->triangles(i).vertices[j])[k];

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


                for(k = 0; k < 3; k++)
                    nodePositions[arrayIdx][k] = cCoords[0][k]*cT.nodes[cN].domainPos()[0]
                                              +cCoords[1][k]*cT.nodes[cN].domainPos()[1]
                                              +cCoords[2][k]*(1 - cT.nodes[cN].domainPos()[0] - cT.nodes[cN].domainPos()[1]);
                
                nodeNumber[arrayIdx]     = cT.nodes[cN].getNodeNumber();
                nodeType[arrayIdx + numVertices] = INTERSECTION_NODE;
//                triId[arrayIdx] = i;
                for(k = 0; k < 3; k++)
                    imagePos[arrayIdx+ numVertices][k] = (par->imagePos(i,cN))[k];
                newIdx[cN] = arrayIdx;
                newIdxlocal[cN] = localArrayIdx;
                arrayIdx++;
                localArrayIdx++;
            }
        }

        for (size_t cN=0; cN<cT.nodes.size(); cN++) {

            if (cT.nodes[cN].isTOUCHING_NODE()){
                for(k = 0; i < 2; k++)
                    domainPositions[arrayIdx][i] = (cT.nodes[cN].domainPos())[k];
                for(k = 0; k < 3; k++)
                    nodePositions[arrayIdx][k] = cCoords[0][k]*cT.nodes[cN].domainPos()[0]
                                              +cCoords[1][k]*cT.nodes[cN].domainPos()[1]
                                              +cCoords[2][k]*(1 - cT.nodes[cN].domainPos()[0] - cT.nodes[cN].domainPos()[1]);
                nodeNumber[arrayIdx]     = cT.nodes[cN].getNodeNumber();
                nodeType[arrayIdx+ numVertices] = TOUCHING_NODE;
                for(k = 0; k < 3; k++)
                    imagePos[arrayIdx+ numVertices][k] = (par->imagePos(i,cN))[k];
//                triId[arrayIdx] = i;
                newIdx[cN] = arrayIdx;
                newIdxlocal[cN] = localArrayIdx;
                arrayIdx++;
                localArrayIdx++;
            }
        }

        for (size_t cN=0; cN<cT.nodes.size(); cN++) {

            if (cT.nodes[cN].isINTERIOR_NODE()){
               for(k = 0; i < 2; k++)
                    domainPositions[arrayIdx][i] = (cT.nodes[cN].domainPos())[k];

                for(k = 0; k < 3; k++)
                (nodePositions[arrayIdx])[k] = cCoords[0][k]*cT.nodes[cN].domainPos()[0]
                                              +cCoords[1][k]*cT.nodes[cN].domainPos()[1]
                                              +cCoords[2][k]*(1 - cT.nodes[cN].domainPos()[0] - cT.nodes[cN].domainPos()[1]);
                nodeNumber[arrayIdx]     = cT.nodes[cN].getNodeNumber();
                nodeType[arrayIdx+ numVertices] = INTERIOR_NODE;
//                triId[arrayIdx] = i;
                for(k = 0; k < 3; k++)
                    imagePos[arrayIdx+ numVertices][k] = (par->imagePos(i,cN))[k];

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
                parameterEdgeArray[edgeArrayIdx][0] = 2;
                parameterEdgeArray[edgeArrayIdx][1] = 2; 
                parameterEdgeArray[edgeArrayIdx][2] = newIdx[cE.from()] + numVertices;
                parameterEdgeArray[edgeArrayIdx][3] = newIdx[cE.to()] + numVertices;

                //////////////////////////////////////////////////////
                parameterEdgeArrayLocal[edgeArrayIdx][0] = newIdxlocal[cE.from()];
                parameterEdgeArrayLocal[edgeArrayIdx][1] = newIdxlocal[cE.to()];
                edgeArrayIdx++;
            }
        }

        // the edgePoints for this triangle
        for (j=0; j<3; j++){
            for (size_t k=1; k<cT.edgePoints[j].size()-1; k++)
                edgePointsArray[edgePointsArrayIdx++] = newIdxlocal[cT.edgePoints[j][k]];
        }
    }              
    
    //6) 'Nodes'
    // barycentric coordinates on the respective triangle and 
    // x, y, and z-coordinates of the image position.
    //nodes on plane surface of triangle
    dimz[0] = numNodes;
    dimz[1] = 3;
    writeFloatDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dimz, &status, "/NodeCoords", nodePositions);
                
    //NodeData                 
    dimz[0] = numNodes;
    dimz[1] = 2;
    writeFloatDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dimz, &status, "/LocalNodePosition", domainPositions);

    //image position
    dimz[0] = nvertices;
    dimz[1] = 3;
    writeFloatDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dimz, &status, "/ImagePosition", imagePos);

    //7)NodeNumbers
    dims[0] = numNodes;
    writeIntDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dims, &status, "/NodeNumber", nodeNumber);

    //8) 'ParameterEdges'
    //connection array(in global index) of parameter edges
    dimz[0] = numParamEdges;
    dimz[1] = 2;
    writeIntDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dimz, &status, "/LocalParamEdge", parameterEdgeArrayLocal);    

    //param edge
    dimz[0] = numParamEdges;
    dimz[1] = 4;
    writeIntDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dimz, &status, "/ParamEdge", parameterEdgeArray);

    //9) 'EdgePoints'
    dims[0] = edgePointsArrayIdx;
    printf("edgePointsArrayIdx = %d\n", edgePointsArrayIdx);
    writeIntDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dims, &status, "/EdgePoints", edgePointsArray);

    //Supportive data
    //params
    int psurfaceparams[4];
    psurfaceparams[0] = numVertices;
    psurfaceparams[1] = numNodes;
    psurfaceparams[2] = numTriangles;
    psurfaceparams[3] = numParamEdges;

    dims[0] = 4;
    writeIntDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dims, &status, "/Params", psurfaceparams);

    //nodetype
    dims[0] = nvertices;
    writeIntDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dims, &status, "/NodeType", nodeType);
    
    //close the file
    status = H5Fclose(file_id);
    return; 
  };

  template<class ctype,int dim>
  ///write the data array into hdf5 data structure(This function store exactly the same data as amiramesh does) 
  void psurface::Hdf5IO<ctype,dim>::writeBaseHdf5Data(const std::string& filename)
  {
    hid_t     file_id;
    file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    hid_t     dataset_id, dataspace_id,datatype;
    hsize_t   dims[] ={numVertices};
    hsize_t   dimz[2];
    herr_t    status;
    
    int i, j, k;
    numVertices  = par->getNumVertices();
    numTriangles = par->getNumTriangles();

    //1) 'Patches'
    int patches_vec[par->patches.size()];
    for(i = 0; i < par->patches.size();i++)
    {
        patches_vec[3*i] = (par->patches[i]).innerRegion;
        patches_vec[3*i + 1] = (par->patches[i]).outerRegion;
        patches_vec[3*i + 2] = (par->patches[i]).boundaryId;
    }
    dims[0] = par->patches.size()*3;
    writeIntDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dims, &status, "/patches_vec", patches_vec);

    //2) 'BaseGridVertexCoords'
    ctype basecoords[numVertices][3];
    for(i = 0; i < numVertices; i++)
    {
      basecoords[i][0] = (par->vertices(i))[0];
      basecoords[i][1] = (par->vertices(i))[1];
      basecoords[i][2] = (par->vertices(i))[2];
    }
    dimz[0] = numVertices;
    dimz[1] = 3;
    writeFloatDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dimz, &status, "/BaseCoords", basecoords);
  
    //3) 'BaseGridTriangles'
    int base_tri[numTriangles][4];
    for(i = 0; i < numTriangles;i++)
    {
      base_tri[i][0] = 4; //topology number of triangle in xdmf
      base_tri[i][1] = par->triangles(i).vertices[0];
      base_tri[i][2] = par->triangles(i).vertices[1];
      base_tri[i][3] = par->triangles(i).vertices[2];
    }
    dimz[0] = numTriangles;
    dimz[1] = 4;
    writeIntDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dimz, &status, "/BaseTri", base_tri);
  
    //4) NodePositions(x, y, and z-coordinates of the image position).
    //ipos
    ctype ipos[par->iPos.size()][3];
    for(i = 0; i < par->iPos.size(); i++)
    {
        ipos[i][0] = (par->iPos[i])[0];
        ipos[i][1] = (par->iPos[i])[1];
        ipos[i][2] = (par->iPos[i])[2];
    }
    dimz[0] = par->iPos.size();
    dimz[1] = 3;
    writeFloatDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dimz, &status, "/iPos", ipos);

    //5) 'NumNodesAndParameterEdgesPerTriangle'
    int num_nodes_and_edges_array[numTriangles][11];
    numNodes      = 0;
    numParamEdges = 0;
    numEdgePoints = 0;
    for(i = 0; i < numTriangles;i++)
    {
        const DomainTriangle<ctype>& cT = par->triangles(i);

        int numIntersectionNodes;
        int numTouchingNodes;
        int numInteriorNodes;

        cT.countNodes(numIntersectionNodes, numTouchingNodes, numInteriorNodes);
        //        int numEdges = cT.getNumRegularEdges() - 3;
        int numEdges = cT.getNumRegularEdges();
        num_nodes_and_edges_array[i][0] = numIntersectionNodes;
        num_nodes_and_edges_array[i][1] = numTouchingNodes;
        num_nodes_and_edges_array[i][2] = numInteriorNodes;
        num_nodes_and_edges_array[i][3] = numEdges;
        num_nodes_and_edges_array[i][4] = cT.patch;

        num_nodes_and_edges_array[i][5] = cT.edgePoints[0].size()-2;
        num_nodes_and_edges_array[i][6] = cT.edgePoints[1].size()-2;
        num_nodes_and_edges_array[i][7] = cT.edgePoints[2].size()-2;

        num_nodes_and_edges_array[i][8] = cT.nodes[cT.cornerNode(0)].getNodeNumber();
        num_nodes_and_edges_array[i][9] = cT.nodes[cT.cornerNode(1)].getNodeNumber();
        num_nodes_and_edges_array[i][10] = cT.nodes[cT.cornerNode(2)].getNodeNumber();

        numNodes += numIntersectionNodes;
        numNodes += numTouchingNodes;
        numNodes += numInteriorNodes;

        numEdgePoints += cT.edgePoints[0].size() + cT.edgePoints[1].size() + cT.edgePoints[2].size() - 6;

        numParamEdges += numEdges;
    }

    dimz[0] = numTriangles;
    dimz[1] = 11;
    writeIntDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dimz, &status, "/numNodesAndEdgesArray", num_nodes_and_edges_array);

    /////////////////////////////////////////////////////////////////////    
    ncells = numTriangles + numParamEdges;
    nvertices = numVertices + numNodes;
    
    
    ctype domainPositions[numNodes][2];
    int   nodeNumber[numNodes];
    int   parameterEdgeArrayLocal[numParamEdges][2];
    int   edgePointsArray[nvertices];


    int arrayIdx           = 0;
    int edgeArrayIdx       = 0;
    int edgePointsArrayIdx = 0;
    ctype cCoords[3][3];
    
    for (i=0; i<numTriangles; i++) {
        const DomainTriangle<ctype>& cT = par->triangles(i);

        std::vector<int> newIdxlocal(cT.nodes.size());
        
        newIdxlocal[cT.cornerNode(0)] = 0;
        newIdxlocal[cT.cornerNode(1)] = 1;
        newIdxlocal[cT.cornerNode(2)] = 2;
        int localArrayIdx = 3;
        for (size_t cN=0; cN<cT.nodes.size(); cN++) {
            if (cT.nodes[cN].isINTERSECTION_NODE()){
              for(k = 0; k < 2; k++)
                    domainPositions[arrayIdx][k] = (cT.nodes[cN].domainPos())[k];

                nodeNumber[arrayIdx]     = cT.nodes[cN].getNodeNumber();
    //            triId[arrayIdx] = i;
                newIdxlocal[cN] = localArrayIdx;
                arrayIdx++;
                localArrayIdx++;
            }
        }

        for (size_t cN=0; cN<cT.nodes.size(); cN++) {

            if (cT.nodes[cN].isTOUCHING_NODE()){
                for(k = 0; k < 2; k++)
                    domainPositions[arrayIdx][k] = (cT.nodes[cN].domainPos())[k];
                
                nodeNumber[arrayIdx]     = cT.nodes[cN].getNodeNumber();
      //          triId[arrayIdx] = i;
                newIdxlocal[cN] = localArrayIdx;
                arrayIdx++;
                localArrayIdx++;
            }
        }

        for (size_t cN=0; cN<cT.nodes.size(); cN++) {

            if (cT.nodes[cN].isINTERIOR_NODE()){
               for(k = 0; k < 2; k++)
                    domainPositions[arrayIdx][k] = (cT.nodes[cN].domainPos())[k];
                
                nodeNumber[arrayIdx]     = cT.nodes[cN].getNodeNumber();
        //        triId[arrayIdx] = i;
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
                parameterEdgeArrayLocal[edgeArrayIdx][0] = newIdxlocal[cE.from()];
                parameterEdgeArrayLocal[edgeArrayIdx][1] = newIdxlocal[cE.to()];
                edgeArrayIdx++;
            }
        }

        // the edgePoints for this triangle
        for (j=0; j<3; j++){
            for (size_t k=1; k<cT.edgePoints[j].size()-1; k++)
                edgePointsArray[edgePointsArrayIdx++] = newIdxlocal[cT.edgePoints[j][k]];
        }
    }              
    
    //6) 'Nodes'
    // barycentric coordinates on the respective triangle and 
    // x, y, and z-coordinates of the image position.
    //NodeData                 
    dimz[0] = numNodes;
    dimz[1] = 2;
    writeFloatDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dimz, &status, "/LocalNodePosition", domainPositions);

    //7)NodeNumbers
    dims[0] = numNodes;
    writeIntDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dims, &status, "/NodeNumber", nodeNumber);

    //8) 'ParameterEdges'
    //connection array(in global index) of parameter edges
    dimz[0] = numParamEdges;
    dimz[1] = 2;
    writeIntDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dimz, &status, "/LocalParamEdge", parameterEdgeArrayLocal);    

    //9) 'EdgePoints'
    dims[0] = edgePointsArrayIdx;
    writeIntDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dims, &status, "/EdgePoints", edgePointsArray);

    //Supportive data
    //params
    int psurfaceparams[4];
    psurfaceparams[0] = numVertices;
    psurfaceparams[1] = numNodes;
    psurfaceparams[2] = numTriangles;
    psurfaceparams[3] = numParamEdges;

    dims[0] = 4;
    writeIntDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dims, &status, "/Params", psurfaceparams);

    //close the file
    status = H5Fclose(file_id);
    return; 
  };

  template<class ctype,int dim>
  ///writhe the xdmf file which store the structure information of hdf5 file.
  void psurface::Hdf5IO<ctype,dim>::writeXdmf(const std::string&  xdf_filename, const std::string&  hdf_filename)
  {
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
      fprintf(xmf, "<DataItem ItemType=\"HyperSlab\" Dimensions=\"%d %d\" Type=\"HyperSlab\">\n",numTriangles,1);

        fprintf(xmf, "<DataItem  Dimensions=\"3 2\" Format=\"XML\">\n");
        fprintf(xmf, "%d %d\n",1, 4);
        fprintf(xmf, "%d %d\n",1, 11);
        fprintf(xmf, "%d %d\n", numTriangles,1);
        fprintf(xmf, "</DataItem>\n"); 

        fprintf(xmf, " <DataItem Dimensions=\"%d %d\" NumberType=\"int\" Format=\"HDF\">\n", numTriangles, 11);
        fprintf(xmf, "%s:/numNodesAndEdgesArray\n", hdf_filename.c_str());
        fprintf(xmf, "</DataItem>\n");

      fprintf(xmf, "</DataItem>\n");
    
    fprintf(xmf, "</Attribute>\n");

    fprintf(xmf, "</Grid>\n");

    fprintf(xmf, "<Grid Name=\"paramedge\" GridType=\"Uniform\">\n");
    fprintf(xmf, "<Topology TopologyType=\"Mixed\" NumberOfElements=\"%d\">\n", numParamEdges);
          fprintf(xmf, "<DataItem Dimensions = \"%d %d\" NumberType=\"int\" Format=\"HDF\">\n",  numParamEdges, 4);
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
    return;
  };

  ///Create PSurface which only have BaseGridVertex and BaseTriangle
/*  template<class ctype,int dim>
  void psurface::Hdf5IO<ctype,dim>::initBaseCasePSurface(Surface* surf, const std::string&  filename)
  {
    PSurfaceFactory<2,ctype> factory(par);
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

    hid_t file;
    hid_t datatype, dataset;
    hid_t filespace;
    hid_t       memspace;
    hsize_t dims[1];
    hsize_t dimz[2];
    herr_t status,status_n;
    int rank;
    int i,j,k;
    file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    //read params
    int *psurfaceparams;
    readIntDataFromFile(&file, &dataset, &filespace, &memspace, dims,"/Params", psurfaceparams);
    numVertices = psurfaceparams[0];
    numNodes = psurfaceparams[1];
    numTriangles = psurfaceparams[2];
    numParamEdges = psurfaceparams[3];
    nvertices = numVertices + numNodes;
    ncells = numTriangles + numParamEdges;
    hdf_close(dataset, filespace, memspace);    
    
    
    ///insert vertex
    StaticVector<ctype,3> newVertex;
    for(int i = 0; i < numVertices; i++)
    {
      for(int j = 0; j < 3; j++) newVertex[j] = (baseGridVertexCoordsArray[i])[j];
      factory.insertVertex(newVertex);
    }

    ///insert image node position
    ctype *coords;
    readFloatDataFromFile(&file, &dataset, &filespace, &memspace, dimz,"/iPos",coords);
    hdf_close(dataset, filespace, memspace);

    psurf->iPos.resize(dimz[0]/3);
    for (size_t i=0; i< dimz[0]/3; i++)
      for (int j=0; j<3; j++)
        psurf->iPos[i][j] = coords[3*i + j];

    ///insert trianlges and the plain graph onto it.
    int edgeCounter=0, edgePointCounter=0;
    int nodeArrayIdx = 0;
    int *tri;
    readIntDataFromFile(&file, &dataset, &filespace, &memspace, dimz,"/BaseTri", tri);
    hdf_close(dataset, filespace, memspace);
    for (int i=0; i<numTriangles; i++)
    {
          std::tr1::array<unsigned int, 3> triangleVertices = { tri[4*i + 1], tri[4*i + 2], tri[4*i + 3]};
          int newTriIdx = factory.insertSimplex(triangleVertices);
          psurf->triangles(newTriIdx).patch = 0;
    }

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
    return;
  } 
 */

  ///Creat Psurface
  template<class ctype,int dim>
  void psurface::Hdf5IO<ctype,dim>::initCompletePSurface(Surface* surf, const std::string&  filename)
  {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      hid_t file;
      hid_t datatype, dataset;
      hid_t filespace;
      hid_t       memspace;
      hsize_t dims[1];
      hsize_t dimz[2];
      herr_t status,status_n;
      int rank;
      int i,j,k;
      file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

      //read params
      int *psurfaceparams;
      readIntDataFromFile(&file, &dataset, &filespace, &memspace, dims,"/Params", psurfaceparams);
      numVertices = psurfaceparams[0];
      numNodes = psurfaceparams[1];
      numTriangles = psurfaceparams[2];
      numParamEdges = psurfaceparams[3];
      nvertices = numVertices + numNodes;
      ncells = numTriangles + numParamEdges;
      hdf_close(dataset, filespace, memspace);

      //read xyz coordinate of vertices and insert vertex into psurface
        
          ctype *tricoords;  
          readFloatDataFromFile(&file, &dataset, &filespace, &memspace, dimz, "/BaseCoords", tricoords);
          hdf_close(dataset, filespace, memspace);
/*          StaticVector<ctype,3> newVertex;
            baseGridVertexCoordsArray.resize(3*numVertices);      
          for(j = 0; j < numVertices ;j++){
              for(i = 0; i < 3; i++)
                  newVertex[j] = tricoords[3*j + i];
          }*/
      
      
      //triangle
      int *baseGridTriArray;
      readIntDataFromFile(&file, &dataset, &filespace, &memspace, dimz,"/BaseTri", baseGridTriArray);
      hdf_close(dataset, filespace, memspace);      
/*      baseGridTriArray.resize(numTriangles);
      for(j = 0; j < numTriangles;j++)
      {
        for(i = 0; i < 3; i++)
          (baseGridTriArray[j])[i] = tri[4*j+ 1 + i];
      } */
      //Parameter Edge
      int *parameterEdgeArrayLocal;
      readIntDataFromFile(&file, &dataset, &filespace, &memspace, dimz,"/LocalParamEdge", parameterEdgeArrayLocal);
      hdf_close(dataset, filespace, memspace);
/*      parameterEdgeArrayLocal.resize(numParamEdges);
      for(j = 0; j < dimz[0];j++)
      {
         for(i = 0; i < 2; i++)
             (parameterEdgeArrayLocal[j])[i] = lparam[2*j + i];
      }*/
      //nodeNumber
      int *nodeNumber;
      readIntDataFromFile(&file, &dataset, &filespace, &memspace, dims,"/NodeNumber", nodeNumber);
      hdf_close(dataset, filespace, memspace);
/*      nodeNumber.resize(numNodes);
      for(i = 0; i < numNodes;i++) nodeNumber[i] = nodenumber[i];*/

      //iPos
      ctype *coords;
      readFloatDataFromFile(&file, &dataset, &filespace, &memspace, dimz,"/iPos",coords);
      hdf_close(dataset, filespace, memspace);
      int iPos_size = dimz[0];
/*      iPos.resize(dimz[0]);
      for(j = 0; j < dimz[0]; j++)
      {
        (iPos[j])[0] = coords[3*j];
        (iPos[j])[1] = coords[3*j+1];
        (iPos[j])[2] = coords[3*j+2];
      }*/

      //nodeNodesAndEdgesArray
      int *numNodesAndEdgesArray;
      readIntDataFromFile(&file, &dataset, &filespace, &memspace, dimz, "/numNodesAndEdgesArray", numNodesAndEdgesArray);
      hdf_close(dataset, filespace, memspace);      
/*      numNodesAndEdgesArray.resize(11*numTriangles);
      for(i = 0; i < 11*numTriangles;i++)
          numNodesAndEdgesArray[i] = nodearray[i];*/

      //local position of nodes on triangle
      ctype *domainPositions;
      readFloatDataFromFile(&file, &dataset, &filespace, &memspace, dimz, "/LocalNodePosition",domainPositions);
      hdf_close(dataset, filespace, memspace);      
/*      domainPositions.resize(numNodes);
      for(i = 0; i < numNodes;i++)
      {
          for(j = 0; j < 2;j++)
              (domainPositions[i])[j] = localpos[2*i + j];
      }*/
      
      //edgepointsarray
      int *edgePointsArray;
      readIntDataFromFile(&file, &dataset, &filespace, &memspace, dims, "/EdgePoints", edgePointsArray);
      hdf_close(dataset, filespace, memspace);      
/*      edgePointsArray.resize(dims[0]);
      for(i = 0; i < dims[0];i++)
      {
         edgePointsArray[i] = edgep[i];
      }*/

      //patches
      int *patch;
      readIntDataFromFile(&file, &dataset, &filespace, &memspace, dims, "/patches_vec", patch);
      hdf_close(dataset, filespace, memspace);
      int patch_size = dims[0]/3;
/*      patches.resize(dims[0]/3);
      for(i = 0; i < dims[0]/3;i++)
      {
        (patches[i]).innerRegion = patch[3*i];
        (patches[i]).outerRegion = patch[3*i + 1];
        (patches[i]).boundaryId = patch[3*i + 2];
      }*/
      H5Fclose(file);
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////      
    /// Create PSurface factory
    PSurfaceFactory<2,ctype> factory(par);
    ///(Assume) Target surface already exists
    factory.setTargetSurface(surf);

    //set param for psurface
    AmiraMesh am;
    am.parameters.set("ContentType", "Parametrization");
    am.parameters.remove(am.parameters[0]);
    par->getPaths(am.parameters);

    //patches
    PSurface<2,float>::Patch Patch; //= new PSurface<2, ctype>::Patch;
    for (i=0; i< patch_size; i++)
    {
      Patch.innerRegion = patch[3*i];
      Patch.outerRegion = patch[3*i + 1];  
      Patch.boundaryId = patch[3*i + 2];
      par->patches.push_back(Patch);
    }

    ///insert vertex
    StaticVector<ctype,3> newVertex;
    for(i = 0; i < numVertices; i++)
    {
       for(int j = 0; j < 3; j++) newVertex[j] = tricoords[3*i + j]; 
       factory.insertVertex(newVertex);
    }

    ///insert image node position
    par->iPos.resize(iPos_size);
    for ( i=0; i< iPos_size; i++)
      for (j=0; j<3; j++)
        par->iPos[i][j] = coords[3*i + j];

    ///insert trianlges and the plain graph onto it.
    int edgeCounter=0, edgePointCounter=0;
    int nodeArrayIdx = 0;

    for (i = 0; i< numTriangles; i++){
        std::tr1::array<unsigned int, 3> triangleVertices = {baseGridTriArray[4*i + 1],baseGridTriArray[4*i + 2],baseGridTriArray[4*i + 3]};
        int newTriIdx = factory.insertSimplex(triangleVertices);
        par->triangles(newTriIdx).patch = numNodesAndEdgesArray[11*i+4];
    /// get the parametrization on this triangle
    int numIntersectionNodes = numNodesAndEdgesArray[11*i+0];
    int numTouchingNodes     = numNodesAndEdgesArray[11*i+1];
    int numInteriorNodes     = numNodesAndEdgesArray[11*i+2];
    int numParamEdges        = numNodesAndEdgesArray[11*i+3];

    ///nodes
    par->triangles(newTriIdx).nodes.resize(numIntersectionNodes + numTouchingNodes + numInteriorNodes + 3);
    int nodenumber;
    // three corner nodes
    StaticVector<ctype,2> domainPos(1, 0);
    nodenumber =  numNodesAndEdgesArray[11*i+8];
    par->triangles(newTriIdx).nodes[0].setValue(domainPos, nodenumber, Node<ctype>::CORNER_NODE);
    par->triangles(newTriIdx).nodes[0].makeCornerNode(0, nodenumber);

    domainPos = StaticVector<ctype,2>(0, 1);
    nodenumber = numNodesAndEdgesArray[11*i+9];
    par->triangles(newTriIdx).nodes[1].setValue(domainPos, nodenumber, Node<ctype>::CORNER_NODE);
    par->triangles(newTriIdx).nodes[1].makeCornerNode(1, nodenumber);

    domainPos = StaticVector<ctype,2>(0, 0);
    nodenumber = numNodesAndEdgesArray[11*i+10];
    par->triangles(newTriIdx).nodes[2].setValue(domainPos, nodenumber, Node<ctype>::CORNER_NODE);
    par->triangles(newTriIdx).nodes[2].makeCornerNode(2, nodenumber);

     int nodeCounter = 3;
    ///the intersection nodes
    for (j=0; j<numIntersectionNodes; j++, nodeCounter++, nodeArrayIdx++){
      domainPos = StaticVector<ctype,2>(domainPositions[2*nodeArrayIdx], domainPositions[2*nodeArrayIdx + 1]);
      nodenumber = nodeNumber[nodeArrayIdx];
      par->triangles(newTriIdx).nodes[nodeCounter].setValue(domainPos, nodenumber, Node<ctype>::INTERSECTION_NODE);

    }

    // the touching nodes
    for (j=0; j<numTouchingNodes; j++, nodeCounter++, nodeArrayIdx++){
      domainPos = StaticVector<ctype,2>(domainPositions[2*nodeArrayIdx], domainPositions[2*nodeArrayIdx + 1]);
      int nodenumber    = nodeNumber[nodeArrayIdx];
      par->triangles(newTriIdx).nodes[nodeCounter].setValue(domainPos, nodenumber, Node<ctype>::TOUCHING_NODE);
    }

    // the interior nodes
    for (j=0; j<numInteriorNodes; j++, nodeCounter++, nodeArrayIdx++){
      domainPos = StaticVector<ctype,2>(domainPositions[2*nodeArrayIdx], domainPositions[2*nodeArrayIdx+1]);
      int nodenumber    = nodeNumber[nodeArrayIdx];
      par->triangles(newTriIdx).nodes[nodeCounter].setValue(domainPos, nodenumber, Node<ctype>::INTERIOR_NODE);
    }

    /// the parameterEdges
    for(j = 0; j < numParamEdges;j++, edgeCounter++)
    {
      par->triangles(newTriIdx).addEdge(parameterEdgeArrayLocal[2*edgeCounter],parameterEdgeArrayLocal[2*edgeCounter + 1]);
    }
    /// the edgePoints arrays on each triangular edge
    for (j=0; j<3; j++){
      par->triangles(newTriIdx).edgePoints[j].resize(numNodesAndEdgesArray[11*i+5+j] + 2);
      par->triangles(newTriIdx).edgePoints[j][0]     = j;
      par->triangles(newTriIdx).edgePoints[j].back() = (j+1)%3;


      for ( k = 0; k<numNodesAndEdgesArray[11*i+5+j]; k++){
          par->triangles(newTriIdx).edgePoints[j][k+1] = edgePointsArray[edgePointCounter];
          edgePointCounter++;
      }
    }
    }
    par->hasUpToDatePointLocationStructure = false;
    par->setupOriginalSurface();
    return;
  };

  ///initialize PsurfaceConvert from the psurface object
  template<class ctype,int dim>
  psurface::Hdf5IO<ctype,dim>::Hdf5IO(PSurface<dim,ctype>* psurface)
  {
    par = psurface;
  }

template class Hdf5IO<float,2>;


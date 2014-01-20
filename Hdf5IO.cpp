#include "config.h"

#include <vector>
#include <string.h>
#include <hdf5.h>
#include <memory>
#include <tr1/memory>
#include <stdexcept>

#include "StaticVector.h"
#include "Domains.h"
#include "PSurface.h"
#include "TargetSurface.h"
#include "PSurfaceFactory.h"
#include "Hdf5IO.h"

using namespace psurface;

//Writes one dimensional int type data array to hdf5 file
void writeIntDataToFile(hid_t* file_id, hid_t* dataset_id, hid_t* dataspace_id, hid_t* datatype, hsize_t* dims, herr_t* status, const char*  name, int address[])
{
    *dataspace_id = H5Screate_simple(1, dims, NULL);
    *dataset_id = H5Dcreate(*file_id, name, H5T_NATIVE_INT, *dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    *status = H5Dwrite(*dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, address);
    *status = H5Dclose(*dataset_id);
    *status = H5Sclose(*dataspace_id);
}

//Writes two dimensional int type data array to hdf5 file
void writeIntDataToFile(hid_t* file_id, hid_t* dataset_id, hid_t* dataspace_id, hid_t* datatype, hsize_t* dims, herr_t* status, const char*  name, int address[][2])
{
    *dataspace_id = H5Screate_simple(2, dims, NULL);
    *dataset_id = H5Dcreate(*file_id, name, H5T_NATIVE_INT, *dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    *status = H5Dwrite(*dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, address);
    *status = H5Dclose(*dataset_id);
    *status = H5Sclose(*dataspace_id);
}

//Writes two dimensional int type data array to hdf5 file
void writeIntDataToFile(hid_t* file_id, hid_t* dataset_id, hid_t* dataspace_id, hid_t* datatype, hsize_t* dims, herr_t* status, const char*  name, int address[][4])
{
    *dataspace_id = H5Screate_simple(2, dims, NULL);
    *dataset_id = H5Dcreate(*file_id, name, H5T_NATIVE_INT, *dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    *status = H5Dwrite(*dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, address);
    *status = H5Dclose(*dataset_id);
    *status = H5Sclose(*dataspace_id);
}

//Writes two dimensional int type data array to hdf5 file
void writeIntDataToFile(hid_t* file_id, hid_t* dataset_id, hid_t* dataspace_id, hid_t* datatype, hsize_t* dims, herr_t* status, const char*  name, int address[][11])
{
    *dataspace_id = H5Screate_simple(2, dims, NULL);
    *dataset_id = H5Dcreate(*file_id, name, H5T_NATIVE_INT, *dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    *status = H5Dwrite(*dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, address);
    *status = H5Dclose(*dataset_id);
    *status = H5Sclose(*dataspace_id);
}

//Writes one dimensional float type data array to hdf5 file
void writeFloatDataToFile(hid_t* file_id, hid_t* dataset_id, hid_t* dataspace_id, hid_t* datatype, hsize_t* dims, herr_t* status, const char* name, float address[])
  {
    *dataspace_id = H5Screate_simple(1, dims, NULL);
    *dataset_id = H5Dcreate(*file_id, name, H5T_NATIVE_FLOAT, *dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    *status = H5Dwrite(*dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, address);
    *status = H5Dclose(*dataset_id);
    *status = H5Sclose(*dataspace_id);
  }

//Writes two dimensional float type data array to hdf5 file
void writeFloatDataToFile(hid_t* file_id, hid_t* dataset_id, hid_t* dataspace_id, hid_t* datatype, hsize_t* dims, herr_t* status, const char* name, float address[][2])
  {
    *dataspace_id = H5Screate_simple(2, dims, NULL);
    *dataspace_id = H5Screate_simple(2, dims, NULL);
    *dataset_id = H5Dcreate(*file_id, name, H5T_NATIVE_FLOAT, *dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    *status = H5Dwrite(*dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, address);
    *status = H5Dclose(*dataset_id);
  }

//Writes two dimensional float type data array to hdf5 file
void writeFloatDataToFile(hid_t* file_id, hid_t* dataset_id, hid_t* dataspace_id, hid_t* datatype, hsize_t* dims, herr_t* status, const char* name, float address[][3])
  {
    *dataspace_id = H5Screate_simple(2, dims, NULL);
    *dataspace_id = H5Screate_simple(2, dims, NULL);
    *dataset_id = H5Dcreate(*file_id, name, H5T_NATIVE_FLOAT, *dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    *status = H5Dwrite(*dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, address);
    *status = H5Dclose(*dataset_id);
  }

//Writes one dimensional double type data array to hdf5 file
  void writeDoubletDataToFile(hid_t* file_id, hid_t* dataset_id, hid_t* dataspace_id, hid_t* datatype, hsize_t* dims, herr_t* status, const char* name, double address[])
  {
    *dataspace_id = H5Screate_simple(1, dims, NULL);
    *dataset_id = H5Dcreate(*file_id, name, H5T_NATIVE_DOUBLE, *dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    *status = H5Dwrite(*dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, address);
    *status = H5Dclose(*dataset_id);
    *status = H5Sclose(*dataspace_id);
  }

//Writes two dimensional double type data array to hdf5 file
  void writeDoubleDataToFile(hid_t* file_id, hid_t* dataset_id, hid_t* dataspace_id, hid_t* datatype, hsize_t* dims, herr_t* status, const char* name, double address[][2])
  {
    *dataspace_id = H5Screate_simple(2, dims, NULL);
    *dataset_id = H5Dcreate(*file_id, name, H5T_NATIVE_DOUBLE, *dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    *status = H5Dwrite(*dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, address);
    *status = H5Dclose(*dataset_id);
    *status = H5Sclose(*dataspace_id);
  }

//Writes two dimensional double type data array to hdf5 file
  void writeDoubleDataToFile(hid_t* file_id, hid_t* dataset_id, hid_t* dataspace_id, hid_t* datatype, hsize_t* dims, herr_t* status, const char* name, double address[][3])
  {
    *dataspace_id = H5Screate_simple(2, dims, NULL);
    *dataset_id = H5Dcreate(*file_id, name, H5T_NATIVE_DOUBLE, *dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    *status = H5Dwrite(*dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, address);
    *status = H5Dclose(*dataset_id);
    *status = H5Sclose(*dataspace_id);
  }

//Reads int type data array from hdf5 file to one dimensional int array
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

//Reads float type data array from hdf5 file to one dimensional float array
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

//Reads double type data array from hdf5 file to one dimensional double array
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
  psurface::PSurface<2, ctype>* psurface::Hdf5IO<ctype,dim>::read(const std::string& filename)
  {
      PSurface<2,float>* p = new PSurface<2,float>;
      Surface* surf = new Surface;
      Hdf5IO<float,dim>* hdf5io = new Hdf5IO<float,dim>(p);
      hdf5io->initCompletePSurface(surf, filename);
      delete hdf5io;
      return p;
  }

  template<class ctype,int dim>
  void psurface::Hdf5IO<ctype,dim>::createHdfAndXdmf(const std::string&  xdf_filename, const std::string&  hdf_filename, bool base)
  {
    if(!base)
    //Creates hdf5 file that contain all data needed to display it in paraview.
    {
      writeHdf5Data(hdf_filename);
      writeXdmf(xdf_filename, hdf_filename);
    }
    else
    // create hdf5 file that only contain necessary data of the psurface.
      writeBaseHdf5Data(hdf_filename);
  }

  template<class ctype,int dim>
  ///write the data array into hdf5 data structure
  void psurface::Hdf5IO<ctype,dim>::writeHdf5Data(const std::string& filename)
  {
    hid_t     file_id;
    file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    hid_t     dataset_id, dataspace_id,datatype;
    hsize_t   dims1d[1];
    hsize_t   dims2d[2];
    herr_t    status;

    numVertices  = par->getNumVertices();
    numTriangles = par->getNumTriangles();

    //2) 'BaseGridVertexCoords'
    ctype basecoords[numVertices][3];
    for (int i = 0; i < numVertices; i++)
    {
      basecoords[i][0] = par->vertices(i)[0];
      basecoords[i][1] = par->vertices(i)[1];
      basecoords[i][2] = par->vertices(i)[2];
    }
    dims2d[0] = numVertices;
    dims2d[1] = 3;
    writeFloatDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dims2d, &status, "/BaseCoords", basecoords);

    //3) 'BaseGridTriangles'
    int base_tri[numTriangles][4];
    for (int i = 0; i < numTriangles;i++)
    {
      base_tri[i][0] = 4; //topology number of triangle in xdmf
      base_tri[i][1] = par->triangles(i).vertices[0];
      base_tri[i][2] = par->triangles(i).vertices[1];
      base_tri[i][3] = par->triangles(i).vertices[2];
    }
    dims2d[0] = numTriangles;
    dims2d[1] = 4;
    writeIntDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dims2d, &status, "/BaseTri", base_tri);

    //4) NodePositions(x, y, and z-coordinates of the image position).
    //ipos
    ctype ipos[par->iPos.size()][3];
    for (size_t i = 0; i < par->iPos.size(); i++)
    {
        ipos[i][0] = par->iPos[i][0];
        ipos[i][1] = par->iPos[i][1];
        ipos[i][2] = par->iPos[i][2];
    }
    dims2d[0] = par->iPos.size();
    dims2d[1] = 3;
    writeFloatDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dims2d, &status, "/iPos", ipos);

    //5) 'NumNodesAndParameterEdgesPerTriangle'
    int num_nodes_and_edges_array[numTriangles][11];
    numNodes      = 0;
    numParamEdges = 0;
    numEdgePoints = 0;

    for (int i = 0; i < numTriangles;i++)
    {
        const DomainTriangle<ctype>& cT = par->triangles(i);

        int numIntersectionNodes;
        int numTouchingNodes;
        int numInteriorNodes;

        cT.countNodes(numIntersectionNodes, numTouchingNodes, numInteriorNodes);

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
        numNodes += 3;  // corner nodes

        numEdgePoints += cT.edgePoints[0].size() + cT.edgePoints[1].size() + cT.edgePoints[2].size() - 6;

        numParamEdges += numEdges;
    }

    //Write numNodesAndEdgesArray into hdf5
    dims2d[0] = numTriangles;
    dims2d[1] = 11;
    writeIntDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dims2d, &status, "/numNodesAndEdgesArray", num_nodes_and_edges_array);

    ncells = numTriangles + numParamEdges;
    nvertices = numVertices + numNodes;


    ctype nodePositions[numNodes][3];
    ctype domainPositions[numNodes][2];
    ctype imagePos[nvertices][3];

    int   nodeNumber[numNodes];
    int   nodeType[nvertices];
    int   parameterEdgeArrayLocal[numParamEdges][2];
    int   parameterEdgeArray[numParamEdges][4];
    int   edgePointsArray[nvertices];

    for (int i = 0; i < numVertices; i++)
        nodeType[i] = Node<ctype>::CORNER_NODE;

    int arrayIdx           = 0;
    int edgeArrayIdx       = 0;
    int edgePointsArrayIdx = 0;

    for (int i=0; i<numTriangles; i++) {
        const DomainTriangle<ctype>& cT = par->triangles(i);

        std::vector<int> newIdx(cT.nodes.size());
        std::vector<int> newIdxlocal(cT.nodes.size());
        StaticVector<ctype,3> cCoords[3];
        for (int j = 0; j < 3; j++)
            cCoords[j] = par->vertices(par->triangles(i).vertices[j]);

        int localArrayIdx = 0;

        //Corner Node
        for (size_t cN=0; cN<cT.nodes.size(); cN++) {

            for (int k = 0; k < 2; k++)
                domainPositions[arrayIdx][k] = (cT.nodes[cN].domainPos())[k];
            for (int k = 0; k < 3; k++)
                nodePositions[arrayIdx][k] = cCoords[0][k]*cT.nodes[cN].domainPos()[0]
                                              +cCoords[1][k]*cT.nodes[cN].domainPos()[1]
                                              +cCoords[2][k]*(1 - cT.nodes[cN].domainPos()[0] - cT.nodes[cN].domainPos()[1]);

            nodeNumber[arrayIdx]     = cT.nodes[cN].getNodeNumber();
            nodeType[arrayIdx + numVertices] = cT.nodes[cN].type;
            for (int k = 0; k < 3; k++)
                imagePos[arrayIdx+ numVertices][k] = par->imagePos(i,cN)[k];
            newIdx[cN] = arrayIdx;
            newIdxlocal[cN] = localArrayIdx;
            arrayIdx++;
            localArrayIdx++;

        }

        // the parameterEdges for this triangle
        typename PlaneParam<ctype>::UndirectedEdgeIterator cE;
        for (cE = cT.firstUndirectedEdge(); cE.isValid(); ++cE){
            if(cE.isRegularEdge())
            {
                //Store ParameterEdge in global index of end points
                parameterEdgeArray[edgeArrayIdx][0] = 2;
                parameterEdgeArray[edgeArrayIdx][1] = 2;
                parameterEdgeArray[edgeArrayIdx][2] = newIdx[cE.from()] + numVertices;
                parameterEdgeArray[edgeArrayIdx][3] = newIdx[cE.to()] + numVertices;

                //Store ParameterEdge in local index of end points
                parameterEdgeArrayLocal[edgeArrayIdx][0] = newIdxlocal[cE.from()];
                parameterEdgeArrayLocal[edgeArrayIdx][1] = newIdxlocal[cE.to()];
                edgeArrayIdx++;
            }
        }

        // the edgePoints for this triangle
        for (int j=0; j<3; j++){
            for (size_t k=1; k<cT.edgePoints[j].size()-1; k++)
                edgePointsArray[edgePointsArrayIdx++] = newIdxlocal[cT.edgePoints[j][k]];
        }
    }

    //6) 'Nodes'
    // barycentric coordinates on the respective triangle and
    // x, y, and z-coordinates of the image position.
    //nodes on plane surface of triangle
    dims2d[0] = numNodes;
    dims2d[1] = 3;
    writeFloatDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dims2d, &status, "/NodeCoords", nodePositions);

    //NodeData
    dims2d[0] = numNodes;
    dims2d[1] = 2;
    writeFloatDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dims2d, &status, "/LocalNodePosition", domainPositions);

    //image position
    dims2d[0] = nvertices;
    dims2d[1] = 3;
    writeFloatDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dims2d, &status, "/ImagePosition", imagePos);

    //7)NodeNumbers
    dims1d[0] = numNodes;
    writeIntDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dims1d, &status, "/NodeNumber", nodeNumber);

    //8) 'ParameterEdges'
    //connection array(in global index) of parameter edges
    dims2d[0] = numParamEdges;
    dims2d[1] = 2;
    writeIntDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dims2d, &status, "/LocalParamEdge", parameterEdgeArrayLocal);

    //param edge
    dims2d[0] = numParamEdges;
    dims2d[1] = 4;
    writeIntDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dims2d, &status, "/ParamEdge", parameterEdgeArray);

    //9) 'EdgePoints'
    dims1d[0] = edgePointsArrayIdx;
    writeIntDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dims1d, &status, "/EdgePoints", edgePointsArray);

    //Supportive data
    //params
    int psurfaceparams[4];
    psurfaceparams[0] = numVertices;
    psurfaceparams[1] = numNodes;
    psurfaceparams[2] = numTriangles;
    psurfaceparams[3] = numParamEdges;

    dims1d[0] = 4;
    writeIntDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dims1d, &status, "/Params", psurfaceparams);

    //nodetype
    dims1d[0] = nvertices;
    writeIntDataToFile(&file_id, &dataset_id, &dataspace_id, &datatype, dims1d, &status, "/NodeType", nodeType);

    //close the file
    status = H5Fclose(file_id);
  };

  template<class ctype,int dim>
  ///write the data array into hdf5 data structure(This function store exactly the same data as amiramesh does)
  void psurface::Hdf5IO<ctype,dim>::writeBaseHdf5Data(const std::string& filename)
  {
    hid_t     file_id;
    file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    hid_t     dataset_id, dataspace_id,datatype;
    hsize_t   dims[1];
    hsize_t   dimz[2];
    herr_t    status;

    int i, j, k;
    numVertices  = par->getNumVertices();
    numTriangles = par->getNumTriangles();

    //2) 'BaseGridVertexCoords'
    ctype basecoords[numVertices][3];
    for(i = 0; i < numVertices; i++)
    {
      basecoords[i][0] = par->vertices(i)[0];
      basecoords[i][1] = par->vertices(i)[1];
      basecoords[i][2] = par->vertices(i)[2];
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
        ipos[i][0] = par->iPos[i][0];
        ipos[i][1] = par->iPos[i][1];
        ipos[i][2] = par->iPos[i][2];
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
                //Store ParameterEdge in local index of end points
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
  };

  template<class ctype,int dim>
  void psurface::Hdf5IO<ctype,dim>::initCompletePSurface(Surface* surf, const std::string&  filename)
  {
      ///////////////////////////////////////////////////////////////////////////////////////////////
      //Read parametrization data from hdf5 file
      hid_t datatype, dataset;
      hid_t filespace;
      hid_t       memspace;
      hsize_t dims[1];
      hsize_t dimz[2];
      herr_t status,status_n;
      int rank;
      int i,j,k;
      hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
      if (file<0)
        throw std::runtime_error("Couldn't open file '" + filename + "' for reading!");


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

      //triangle
      int *baseGridTriArray;
      readIntDataFromFile(&file, &dataset, &filespace, &memspace, dimz,"/BaseTri", baseGridTriArray);
      hdf_close(dataset, filespace, memspace);

      //Parameter Edge
      int *parameterEdgeArrayLocal;
      readIntDataFromFile(&file, &dataset, &filespace, &memspace, dimz,"/LocalParamEdge", parameterEdgeArrayLocal);
      hdf_close(dataset, filespace, memspace);

      //nodeNumber
      int *nodeNumber;
      readIntDataFromFile(&file, &dataset, &filespace, &memspace, dims,"/NodeNumber", nodeNumber);
      hdf_close(dataset, filespace, memspace);

      //iPos
      ctype *coords;
      readFloatDataFromFile(&file, &dataset, &filespace, &memspace, dimz,"/iPos",coords);
      hdf_close(dataset, filespace, memspace);
      int iPos_size = dimz[0];

      //nodeNodesAndEdgesArray
      int *numNodesAndEdgesArray;
      readIntDataFromFile(&file, &dataset, &filespace, &memspace, dimz, "/numNodesAndEdgesArray", numNodesAndEdgesArray);
      hdf_close(dataset, filespace, memspace);

      //local position of nodes on triangle
      ctype *domainPositions;
      readFloatDataFromFile(&file, &dataset, &filespace, &memspace, dimz, "/LocalNodePosition",domainPositions);
      hdf_close(dataset, filespace, memspace);

      //edgepointsarray
      int *edgePointsArray;
      readIntDataFromFile(&file, &dataset, &filespace, &memspace, dims, "/EdgePoints", edgePointsArray);
      hdf_close(dataset, filespace, memspace);

      H5Fclose(file);

      //////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Create psurface object from parametrization data
      // Create PSurface factory
      PSurfaceFactory<2,ctype> factory(par);
      //(Assume) Target surface already exists
      factory.setTargetSurface(surf);

      //insert vertex
      StaticVector<ctype,3> newVertex;
      for(i = 0; i < numVertices; i++)
      {
          for(int j = 0; j < 3; j++)
              newVertex[j] = tricoords[3*i + j];
          factory.insertVertex(newVertex);
      }

      //insert image node position
      par->iPos.resize(iPos_size);
      for ( i=0; i< iPos_size; i++)
          for (j=0; j<3; j++)
              par->iPos[i][j] = coords[3*i + j];

      //insert triangles and the plain graph onto it.
      int edgeCounter=0, edgePointCounter=0;
      int nodeArrayIdx = 0;

      for (i = 0; i< numTriangles; i++){
          std::tr1::array<unsigned int, 3> triangleVertices = {static_cast<unsigned int>(baseGridTriArray[4*i + 1]),
                                                               static_cast<unsigned int>(baseGridTriArray[4*i + 2]),
                                                               static_cast<unsigned int>(baseGridTriArray[4*i + 3])};

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

          //the intersection nodes
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

        // the parameterEdges
        for(j = 0; j < numParamEdges;j++, edgeCounter++)
        {
            par->triangles(newTriIdx).addEdge(parameterEdgeArrayLocal[2*edgeCounter],parameterEdgeArrayLocal[2*edgeCounter + 1]);
        }

        // the edgePoints arrays on each triangular edge
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
  };

  //initialize PsurfaceConvert from the psurface object
  template<class ctype,int dim>
  psurface::Hdf5IO<ctype,dim>::Hdf5IO(PSurface<2,ctype>* psurface)
  {
    par = psurface;
  }

//   Explicit template instantiations.
namespace psurface {
  template class Hdf5IO<float,1>;
  template class Hdf5IO<float,2>;
}

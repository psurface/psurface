//the code to write psurface object into vtk file
//vtuwrite
//#ifndef HDF5IO_H
//#define HDF5IO_H
#include <vector>
#include "StaticVector.h"
#include "Domains.h"
#include <stdio.h>
#include <stdlib.h>
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
class gmsh{
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
/*  bool readGmsh(File *file)
  {
      number_of_real_vertices = 0;
      element_count = 0;

      // process header
      double version_number;
      int file_type, data_size;

      readfile(file,1,"%s\n",buf);
      if (strcmp(buf,"$MeshFormat")!=0)
        printf("expected $MeshFormat in first line\n");
      readfile(file,3,"%lg %d %d\n",&version_number,&file_type,&data_size);
      if( (version_number < 2.0) || (version_number > 2.2) )
        printf("can only read Gmsh version 2 files\n");
      readfile(file,1,"%s\n",buf);
      if (strcmp(buf,"$EndMeshFormat")!=0)
        printf("expected $EndMeshFormat\n");

      // node section
      int number_of_nodes;
      readfile(file,1,"%s\n",buf);
      if (strcmp(buf,"$Nodes")!=0)
        printf("expected $Nodes\n");
      readfile(file,1,"%d\n",&number_of_nodes);

      // read nodes
      int id;
      double x[ 3 ];
      for( int i = 1; i <= number_of_nodes; ++i )
      {
          readfile(file,4, "%d %lg %lg %lg\n", &id, &x[ 0 ], &x[ 1 ], &x[ 2 ] );
          if( id != i )
            DUNE_THROW( Dune::IOError, "Expected id " << i << "(got id " << id << "." );

          // just store node position
          for( int j = 0; j < dimWorld; ++j )
            nodes[ i ][ j ] = x[ j ];
      }
      readfile(file,1,"%s\n",buf);
      if (strcmp(buf,"$EndNodes")!=0)
        DUNE_THROW(Dune::IOError, "expected $EndNodes");
      

      // element section
      readfile(file,1,"%s\n",buf);
      if (strcmp(buf,"$Elements")!=0)
        DUNE_THROW(Dune::IOError, "expected $Elements");
      int number_of_elements;
      readfile(file,1,"%d\n",&number_of_elements);
      if (verbose) std::cout << "file contains " << number_of_elements << " elements" << std::endl;

      //=========================================
      // Pass 1: Renumber needed vertices
      //=========================================

      long section_element_offset = ftell(file);
      std::map<int,unsigned int> renumber;
      for (int i=1; i<=number_of_elements; i++)
      {
        int id, elm_type, number_of_tags;
        readfile(file,3,"%d %d %d ",&id,&elm_type,&number_of_tags);
        for (int k=1; k<=number_of_tags; k++)
        {
          int blub;
          readfile(file,1,"%d ",&blub);
          // k == 1: physical entity (not used here)
          // k == 2: elementary entity (not used here either)
          // if version_number < 2.2:
          //   k == 3: mesh partition 0
          // else
          //   k == 3: number of mesh partitions
          //   k => 4: mesh partition k-4
        }
        pass1HandleElement(file, elm_type, renumber, nodes);
      }
      if (verbose) std::cout << "number of real vertices = " << number_of_real_vertices << std::endl;
      if (verbose) std::cout << "number of boundary elements = " << boundary_element_count << std::endl;
      if (verbose) std::cout << "number of elements = " << element_count << std::endl;
      readfile(file,1,"%s\n",buf);
      if (strcmp(buf,"$EndElements")!=0)
        DUNE_THROW(Dune::IOError, "expected $EndElements");
      boundary_id_to_physical_entity.resize(boundary_element_count);
      element_index_to_physical_entity.resize(element_count);

    
  }
*/
  void readfile(FILE * file, int cnt, const char * format,
                  void* t1, void* t2 = 0, void* t3 = 0, void* t4 = 0,
                  void* t5 = 0, void* t6 = 0, void* t7 = 0, void* t8 = 0,
                  void* t9 = 0, void* t10 = 0)
  {
      off_t pos = ftello(file);
      int c = fscanf(file, format, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10);
      if (c != cnt)
          printf("error in readfile\n");
  }

  bool writeGmsh(const char* gmshFileName)
  {
      FILE *gm = 0;
      gm = fopen(gmshFileName, "w");
      fprintf(gm,"$MeshFormat\n");
      fprintf(gm, "%1.1g %d %d\n",2.2, 0, 8);
      fprintf(gm, "$EndMeshFormat\n");
      fprintf(gm, "$Nodes\n");
      fprintf(gm, "%d\n", nvertices);
      for( int i = 0; i < baseGridVertexCoordsArray.size(); i++)
          fprintf(gm, "%d %f %f %f\n", i, (baseGridVertexCoordsArray[i])[0], (baseGridVertexCoordsArray[i])[1],(baseGridVertexCoordsArray[i])[2]);
      for( int i = 0; i <  nodePositions.size(); i++)
          fprintf(gm, "%d %f %f %f\n", i + baseGridVertexCoordsArray.size(), (nodePositions[i])[0], (nodePositions[i])[1],(nodePositions[i])[2]);
     fprintf(gm, "$EndNodes\n");
     fprintf(gm, "$Elements\n");
     fprintf(gm, "%d\n", ncells);
     for(int i = 0; i < baseGridTriArray.size(); i++)
        fprintf(gm, "%d %d %d %d %d %d %d %d %d\n", i, 2, 3, 1, 1, 1, (baseGridTriArray[i])[0], (baseGridTriArray[i])[1], (baseGridTriArray[i])[2]);
     for(int i = 0; i < parameterEdgeArray.size(); i++)
         fprintf(gm, "%d %d %d %d %d %d %d %d\n", i + parameterEdgeArray.size(), 1, 3, 1, 1, 1, (parameterEdgeArray[i])[0], (parameterEdgeArray[i])[1]);
     fprintf(gm,"$EndElements\n");
     fprintf(gm,"$NodeData\n");
     //string-tag
     fprintf(gm,"1\n");
     fprintf(gm, "\"NodeType\"\n");
     //real tag
     fprintf(gm, "1\n");
     fprintf(gm, "%f\n", 0.0);
     //integer tag
     fprintf(gm, "%d\n", 3); // interger tags
     fprintf(gm, "%d\n", 0); // the time step (0; time steps always start at 0)
     fprintf(gm, "%d\n", 1); // 1-component (scalar) field
     //node type
     fprintf(gm, "%d\n", nvertices);
     for( int i = 0; i < nodeNumber.size(); i++)
        fprintf(gm,"%d %d\n",i, nodeNumber[i]);
     //nodeNumber
     fprintf(gm, "%d\n", nvertices);
     for( int i = 0; i < baseGridVertexCoordsArray.size(); i++)
        fprintf(gm, "%d %d\n", i, 0);
     for( int i = 0; i <  nodePositions.size(); i++)
        fprintf(gm, "%d %d\n",i, nodeNumber[i]);
     //local position
     fprintf(gm, "%d\n", nvertices);
     for( int i = 0; i < baseGridVertexCoordsArray.size(); i++)
        fprintf(gm, "%d %f\n", i, 0);
     for( int i = 0; i <  nodePositions.size(); i++)
        fprintf(gm, "%d %f\n",i, nodeNumber[i]);             
     //parameter_edge_array_local
     //ElementData
     //iPos
     fprintf(gm, "$ElementData\n");
     //string-tag
     fprintf(gm,"1\n");
     fprintf(gm, "\"A scalar view\"\n");
     //real tag
     fprintf(gm, "1\n");
     fprintf(gm, "%f\n", 0.0);
     //integer tag
     fprintf(gm, "%d\n", 3); // interger tags
     fprintf(gm, "%d\n", 0); // the time step (0; time steps always start at 0)
     fprintf(gm, "%d\n", 1); // 1-component (scalar) field
     //patch
     fprintf(gm, "%d\n", patches.size());
     for(int i = 0; i < patches.size(); i++)
         fprintf(gm, "%d %d\n", i, patches[i]);
     //edgepoints array
     //numNodesandEdges array
     fprintf(gm, "$EndElementData\n");
  }
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



/*  bool initFromHDF5(PSurface<2,ctype>* psurf, Surface* surf, const char* filename)
  {
    readHdf5Data(filename);
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
  };*/
};
}

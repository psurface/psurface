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
#include "hdf5IO.h"
#ifdef PSURFACE_STANDALONE
#include "TargetSurface.h"
#else
#include <hxsurface/Surface.h>
#endif
#include <stdexcept>

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
  /**@name data structure to store psurface nodes, trianles and edges on origrinal surface*/
  //@{
  ///the position of corner nodes in global coordinate
  std::vector<StaticVector<ctype,3> > baseGridVertexCoordsArray;

  ///triangles
  std::vector<StaticVector<int, 3> > baseGridTriArray;

  ///edges
  std::vector<StaticVector<int,2> > parameterEdgeArray;
  //@}

  int numVertices;
  int numTriangles;

  bool readGmsh(FILE* file)
  {
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

  void skipline(FILE * file)
  {
    int c;
    do {
      c = fgetc(file);
    } while(c != '\n' && c != EOF);
  }

  bool initFromGmsh(PSurface<2,ctype>* psurf, Surface* surf, const char* filename)
  {
    FILE* file = fopen(filename,"r");
    bool i = readGmsh(file);
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
  };

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

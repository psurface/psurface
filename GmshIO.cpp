#include <vector>
#include <string.h>
#include <stdexcept>
#include "StaticVector.h"
#include "Domains.h"
#include "PSurface.h"
#include "PSurfaceFactory.h"
#include "GmshIO.h"


using namespace psurface;
  //initialize PsurfaceConvert from the psurface object
  template<class ctype,int dim>
  psurface::GmshIO<ctype,dim>::GmshIO(PSurface<dim,ctype>* psurface)
  {
    par = psurface;
  }


  //read psurface_convert from Gmsh file
  template<class ctype,int dim>
  void psurface::GmshIO<ctype,dim>::readGmsh(Surface* surf, const std::string&  filename)
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

        StaticVector<int, 3> elementDofs(3);

        readfile(file, 3, "%d %d %d", &(elementDofs[0]),&(elementDofs[1]),&(elementDofs[2]));
        skipline(file);

        triArray.push_back(elementDofs);
      }

      //remove vertices that are not corners of a triangle
      std::vector<bool> nodeInTri(number_of_nodes);
      std::vector<int>  newNodeIndex(number_of_nodes);

      std::fill(nodeInTri.begin(), nodeInTri.end(), false);

      for(int i = 0; i < triArray.size(); i++)
      {
          nodeInTri[(triArray[i])[0] - 1] = true;
          nodeInTri[(triArray[i])[1] - 1] = true;
          nodeInTri[(triArray[i])[2] - 1] = true;
      }
      int newIndx = 0;
      for(int i = 0; i < number_of_nodes; i++)
          newNodeIndex[i] = (nodeInTri[i]) ? newIndx++ : -1;

      //creat parace based on the base triangles
      PSurfaceFactory<2,ctype> factory(par);
      factory.setTargetSurface(surf);

      //patches
      PSurface<2,float>::Patch Patch; //= new PSurface<2, ctype>::Patch;
      Patch.innerRegion = 0;
      Patch.outerRegion = 1;
      Patch.boundaryId =  2;
      par->patches.push_back(Patch);
      Patch.innerRegion = 0;
      Patch.outerRegion = 1;
      Patch.boundaryId =  1;
      par->patches.push_back(Patch);

      ///insert vertex
      int numVertices = 0;
      for(int i = 0; i < number_of_nodes; i++)
      {
          if(nodeInTri[i])
          {
            factory.insertVertex(coordsArray[i]);
            numVertices ++;
          }
      }

      ///insert image node position
      for (int i=0; i< number_of_nodes; i++)
      {
          if(nodeInTri[i])
              par->iPos.push_back(coordsArray[i]);
      }
      ///insert triangles and the plane graph on them
      int edgeCounter=0, edgePointCounter=0;
      int nodeArrayIdx = 0;
      for (int i=0; i< triArray.size(); i++){
          std::tr1::array<unsigned int, 3> triangleVertices = {newNodeIndex[(triArray[i])[0] - 1], newNodeIndex[(triArray[i])[1] - 1], newNodeIndex[(triArray[i])[2] - 1]};
          int newTriIdx = factory.insertSimplex(triangleVertices);
          par->triangles(newTriIdx).patch = 0;

          /// get the parametrization on this triangle
          ///nodes
          par->triangles(newTriIdx).nodes.resize(3);
          int nodenumber;
          // three corner nodes
          StaticVector<ctype,2> domainPos(1, 0);
          nodenumber =  0;
          par->triangles(newTriIdx).nodes[0].setValue(domainPos, nodenumber, Node<ctype>::CORNER_NODE);
          par->triangles(newTriIdx).nodes[0].makeCornerNode(0, nodenumber);

          domainPos = StaticVector<ctype,2>(0, 1);
          nodenumber = 1;
          par->triangles(newTriIdx).nodes[1].setValue(domainPos, nodenumber, Node<ctype>::CORNER_NODE);
          par->triangles(newTriIdx).nodes[1].makeCornerNode(1, nodenumber);

          domainPos = StaticVector<ctype,2>(0, 0);
          nodenumber = 2;
          par->triangles(newTriIdx).nodes[2].setValue(domainPos, nodenumber, Node<ctype>::CORNER_NODE);
          par->triangles(newTriIdx).nodes[2].makeCornerNode(2, nodenumber);

          /// the parameterEdges
          for(int j = 0; j < 3; j++)
              par->triangles(newTriIdx).addEdge(j, (j+1)%3);

          /// the edgePoints arrays on each triangular edge
          for (int j=0; j<3; j++){
              par->triangles(newTriIdx).edgePoints[j].resize(2);
              par->triangles(newTriIdx).edgePoints[j][0]     = j;
              par->triangles(newTriIdx).edgePoints[j].back() = (j+1)%3;
          }
      }

      par->hasUpToDatePointLocationStructure = false;
      par->setupOriginalSurface();
  }

  template<class ctype,int dim>
  void psurface::GmshIO<ctype,dim>::skipline(FILE * file)
  {
    int c;
    do {
      c = fgetc(file);
    } while(c != '\n' && c != EOF);
  }

  template<class ctype,int dim>
  void psurface::GmshIO<ctype,dim>::readfile(FILE * file, int cnt, const char* format,
                  void* t1, void* t2 = 0, void* t3 = 0, void* t4 = 0,
                  void* t5 = 0, void* t6 = 0, void* t7 = 0, void* t8 = 0,
                  void* t9 = 0, void* t10 = 0)
  {
      off_t pos = ftello(file);
      int c = fscanf(file, format, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10);
      if (c != cnt)
          throw(std::runtime_error("error in readfile\n"));
  }

//   Explicit template instantiations.
template class GmshIO<float,2>;


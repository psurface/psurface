#include <vector>
#include <string.h>
#include <stdexcept>
#include "StaticVector.h"
#include "Domains.h"
#include "PSurface.h"
#include "PSurfaceFactory.h"
#include "GmshIO.h"


#ifdef PSURFACE_STANDALONE
#include "TargetSurface.h"
#else
#include "hxsurface/Surface.h"
#endif

  //read psurface_convert from Gmsh file
  template<class ctype,int dim>
  psurface::PSurface<dim, ctype>* psurface::GmshIO<ctype,dim>::readGmsh(const std::string& filename)
  {
      PSurface<dim, ctype>* par = new PSurface<dim, ctype>;
      par->surface = new Surface;

      FILE* file = fopen(filename.c_str(),"r");
      if (not file)
          throw(std::runtime_error("Could open file '" + filename + "' for reading!"));

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

      bool patches = false;
      std::vector<int> patchNums;

      for (int i=1; i<=number_of_elements; i++)
      {
        int id, elm_type, number_of_tags;
        readfile(file,3,"%d %d %d ",&id,&elm_type,&number_of_tags);

        // handle tags
        int new_patch;

        // find the patch this element is belonging to
        if (number_of_tags > 0) {
          // if this is not the first entry and we did not find any patches before, which patches are the elements before this one belonging to ?
          if (1 != i && !patches)
            throw std::runtime_error("found elements with and without tags\n");

          // we found a patch, thous we are using patches
          patches = true;

          // first tag is the patch number
          readfile(file, 1, "%d ", &new_patch);

          // other tags will be ignored
          for (int k=2; k<=number_of_tags; ++k) {
              int dummy;
              readfile(file, 1, "%d ", &dummy);
          }
        } else {
          // if we did not any find tags though we found some before, which patch is this element belonging to ?
          if (patches)
            throw std::runtime_error("found elements with and without tags\n");

          // if no patch numbers are provided, assign them all to the first
          new_patch = 0;
        }

        // store new patch
        patchNums.push_back(new_patch);

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

      fclose(file); // we got everything that we need

      //remove vertices that are not corners of a triangle
      std::vector<bool> nodeInTri(number_of_nodes);
      std::vector<int>  newNodeIndex(number_of_nodes);

      std::fill(nodeInTri.begin(), nodeInTri.end(), false);

      for(size_t i = 0; i < triArray.size(); i++)
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
      factory.setTargetSurface(par->surface);

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
      for (size_t i=0; i<triArray.size(); i++){

          std::tr1::array<int, 3> vertexIdx;

          for (int j=0; j<3; j++)
              vertexIdx[j] = newNodeIndex[triArray[i][j]-1];

          int newTriangle = par->createSpaceForTriangle(vertexIdx[0], vertexIdx[1], vertexIdx[2]);

          par->triangles(newTriangle).makeOneTriangle(vertexIdx[0], vertexIdx[1], vertexIdx[2]);

          par->triangles(newTriangle).patch = patchNums[i];

          par->integrateTriangle(newTriangle);
      }

      par->hasUpToDatePointLocationStructure = false;
      par->setupOriginalSurface();

      return par;
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
                                             void* t1, void* t2, void* t3, void* t4,
                                             void* t5 , void* t6, void* t7, void* t8,
                                             void* t9 , void* t10 )
  {
      int c = fscanf(file, format, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10);
      if (c != cnt)
          throw(std::runtime_error("error in readfile\n"));
  }

//   Explicit template instantiations.
namespace psurface {
  template class GmshIO<float,2>;
}

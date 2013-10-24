#define  HAVE_AMIRAMESH
#include <vector>
#include <string.h>
#include <memory>
#include <tr1/memory>
#include <fstream>

#include "StaticVector.h"
#include "Domains.h"
#include "PSurface.h"
#include "PSurfaceFactory.h"
#include "VtkIO.h"
#include "vtuwriter.hh"


  template<class ctype,int dim>
  psurface::VTKIO<ctype,dim>::VTKIO(PSurface<dim,ctype>* psurface)
  {
    par = psurface;

    numVertices  = par->getNumVertices();
    numTriangles = par->getNumTriangles();

    numNodes      = 0;
    numParamEdges = 0;

    for (int i=0; i<numTriangles; i++) {
        numNodes      += par->triangles(i).nodes.size();
        numParamEdges += par->triangles(i).getNumRegularEdges();
    }

    ncells = numTriangles + numParamEdges;
    nvertices = numVertices + numNodes;
    /////////////////////////////////////////////////////////////////////

    //plane graph on each base grid triangle, saved as a list of nodes and a list of edges.
    int arrayIdx           = 0;
    int edgeArrayIdx       = 0;

    nodeType.resize(numVertices + numNodes);
    nodePositions.resize(numNodes);
    parameterEdgeArray.resize(numParamEdges);
    for (int i = 0; i < numVertices; i++)
        nodeType[i] = Node<ctype>::CORNER_NODE;

    for (int i=0; i<numTriangles; i++) {
        const DomainTriangle<ctype>& cT = par->triangles(i);

        std::vector<int> newIdx(cT.nodes.size());

        // Copy triangle vertex coordinates, for easier-to-write access
        StaticVector<ctype,3> cCoords[3];
        for(int j = 0; j < 3; j++)
            cCoords[j] = par->vertices(par->triangles(i).vertices[j]);

        for (size_t cN=0; cN<cT.nodes.size(); cN++) {

                for (int k = 0; k < 3; k++)
                (nodePositions[arrayIdx])[k] = cCoords[0][k]*cT.nodes[cN].domainPos()[0]
                                              +cCoords[1][k]*cT.nodes[cN].domainPos()[1]
                                              +cCoords[2][k]*(1 - cT.nodes[cN].domainPos()[0] - cT.nodes[cN].domainPos()[1]);
                nodeType[arrayIdx+ numVertices] = cT.nodes[cN].type;
                newIdx[cN] = arrayIdx;
                arrayIdx++;

        }

        // the parameterEdges for this triangle
        typename PlaneParam<ctype>::UndirectedEdgeIterator cE;
        for (cE = cT.firstUndirectedEdge(); cE.isValid(); ++cE){
            if(cE.isRegularEdge())
            {
                parameterEdgeArray[edgeArrayIdx][0] = newIdx[cE.from()] + numVertices;
                parameterEdgeArray[edgeArrayIdx][1] = newIdx[cE.to()] + numVertices;
                edgeArrayIdx++;
            }
        }
    }
  };

  //write the psurface into vtu file
  template<class ctype,int dim>
  void psurface::VTKIO<ctype,dim>::createVTU(std::string filename, bool basegrid)
  {
    std::ofstream file;
    file.open(filename.c_str());
    if (! file.is_open())
      std::cout << filename << "does not exist!" << std::endl;
    writeDataFile(file, basegrid);
    file.close();
  }

  //write data file to stream
  template<class ctype,int dim>
  void psurface::VTKIO<ctype,dim>::writeDataFile(std::ostream& s, bool basegrid)
  {
    VTK::FileType fileType = VTK::unstructuredGrid;

    ///output vtu type
    VTK::OutputType outputtype = VTK::ascii;
    VTK::VTUWriter writer(s, outputtype,fileType);//Most inportant structure used here

    if(basegrid)
      writer.beginMain(numTriangles, numVertices);
    else
      writer.beginMain(numTriangles + numParamEdges, numVertices + numNodes);

    writeAllData(writer, basegrid);
    writer.endMain();
  }

  //write the data section in vtu
  template<class ctype,int dim>
  void psurface::VTKIO<ctype,dim>::writeAllData(VTK::VTUWriter& writer, bool basegrid) {
    //PointData
    writePointData(writer,basegrid);
    // Points
    writeGridPoints(writer,basegrid);
    // Cells
    writeGridCells(writer,basegrid);
  }

  // write point data
  template<class ctype,int dim>
  void psurface::VTKIO<ctype,dim>::writePointData(VTK::VTUWriter& writer, bool basegrid)
  {
      std::string scalars = "nodetype";
      std::string vectors= "";
      int numpoints;
      if(basegrid)
          numpoints = numVertices;
      else
          numpoints = nvertices;

      writer.beginPointData(scalars, vectors);

      {
            std::tr1::shared_ptr<VTK::DataArrayWriter<ctype> > p
            (writer.makeArrayWriter<ctype>(scalars, 1, numpoints));
            for( int i = 0; i < numpoints;i++)
                p->write(nodeType[i]);
      }

      writer.endPointData();
  }

  // write the positions of vertices
  template<class ctype,int dim>
  void psurface::VTKIO<ctype,dim>::writeGridPoints(VTK::VTUWriter& writer, bool basegrid)
  {
      int numpoints;
      if(basegrid)
          numpoints = numVertices;
      else
          numpoints = nvertices;

      writer.beginPoints();
      {
            std::tr1::shared_ptr<VTK::DataArrayWriter<ctype> > p
            (writer.makeArrayWriter<ctype>("Coordinates", 3, numpoints));
            if(!p->writeIsNoop()) {
                  for(int i = 0; i < numVertices; i++)
                        for(int l = 0; l < 3; l++)
                            p->write((par->vertices(i))[l]);

                  if(!basegrid)
                  {
                        for(int i = 0; i < numNodes; i++)
                              for(int l = 0; l < 3; l++)
                                    p->write((nodePositions[i])[l]);
                  }
            }
      }
      writer.endPoints();
  }

  // write the connectivity array
  template<class ctype,int dim>
  void psurface::VTKIO<ctype,dim>::writeGridCells(VTK::VTUWriter& writer, bool basegrid)
  {
      int num_cell = (basegrid) ? numTriangles : ncells;

      writer.beginCells();
      // connectivity
      {
          std::tr1::shared_ptr<VTK::DataArrayWriter<int> > p1
          (writer.makeArrayWriter<int>("connectivity", 1, 3*numTriangles));
          if(!p1->writeIsNoop())
          {
              for(int i = 0; i < numTriangles; i++)
                  for( int l = 0; l < 3; l++)
                      p1->write(par->triangles(i).vertices[l]);

              if(!basegrid)
              {
                  for (size_t i = 0; i <  parameterEdgeArray.size(); i++)
                      for(int l = 0; l < 2; l++)
                          p1->write((parameterEdgeArray[i])[l]);
              }
          }
      }

      // offsets
      {
          std::tr1::shared_ptr<VTK::DataArrayWriter<int> > p2
          (writer.makeArrayWriter<int>("offsets", 1, numTriangles));
          if(!p2->writeIsNoop()) {
              int offset = 0;
              for(int i = 0; i < numTriangles; i++)
              {
                  offset += 3;
                  p2->write(offset);
              }
              if(!basegrid)
              {
                  for (size_t i = 0; i < parameterEdgeArray.size(); i++)
                  {
                      offset += 2;
                      p2->write(offset);
                  }
              }
          }
      }

      // types
      {
          std::tr1::shared_ptr<VTK::DataArrayWriter<unsigned char> > p3
          (writer.makeArrayWriter<unsigned char>("types", 1, num_cell));
          if(!p3->writeIsNoop())
          {
              for(int i = 0; i < numTriangles;i++)
              p3->write(5); //vtktype of triangle
              if(!basegrid)
              {
                  for (size_t i = 0; i < parameterEdgeArray.size(); i++)
                      p3->write(3);//vtktype of edges
              }
          }
      }

      writer.endCells();
  }

//   Explicit template instantiations.
namespace psurface {
  template class VTKIO<float,2>;
}

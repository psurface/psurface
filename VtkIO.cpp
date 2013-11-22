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

    /////////////////////////////////////////////////////////////////////

    //plane graph on each base grid triangle, saved as a list of nodes and a list of edges.
    int arrayIdx           = 0;
    int edgeArrayIdx       = 0;

    nodeType.resize(numNodes);
    nodePositions.resize(numNodes);
    parameterEdgeArray.resize(numParamEdges);

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
                nodeType[arrayIdx] = cT.nodes[cN].type;
                newIdx[cN] = arrayIdx;
                arrayIdx++;

        }

        // the parameterEdges for this triangle
        typename PlaneParam<ctype>::UndirectedEdgeIterator cE;
        for (cE = cT.firstUndirectedEdge(); cE.isValid(); ++cE){
            if(cE.isRegularEdge())
            {
                parameterEdgeArray[edgeArrayIdx][0] = newIdx[cE.from()];
                parameterEdgeArray[edgeArrayIdx][1] = newIdx[cE.to()];
                edgeArrayIdx++;
            }
        }
    }
  };

  //write the psurface into vtu file
  template<class ctype,int dim>
  void psurface::VTKIO<ctype,dim>::createVTU(const std::string& element_filename, const std::string& graph_filename)
  {
    std::ofstream element_file;
    element_file.open(element_filename.c_str());
    if (! element_file.is_open())
      std::cout << "Could not create " << element_filename << std::endl;
    writeElementDataFile(element_file);
    element_file.close();

    if (!graph_filename.empty()) {
      std::ofstream graph_file;
      graph_file.open(graph_filename.c_str());
      if (! graph_file.is_open())
        std::cout << "Could not create " << graph_filename << std::endl;
      writeGraphDataFile(graph_file);
      graph_file.close();
    }
  }

  //write data file to stream
  template<class ctype,int dim>
  void psurface::VTKIO<ctype,dim>::writeElementDataFile(std::ostream& s)
  {
    VTK::FileType fileType = VTK::unstructuredGrid;
    VTK::OutputType outputtype = VTK::ascii;
    VTK::VTUWriter writer(s, outputtype, fileType);

    writer.beginMain(numTriangles, numVertices);

    // Points
    writeElementGridPoints(writer);
    // Cells
    writeElementGridCells(writer);
    // Cell data
    writeElementGridCellData(writer);

    writer.endMain();
  }

  //write data file to stream
  template<class ctype,int dim>
  void psurface::VTKIO<ctype,dim>::writeGraphDataFile(std::ostream& s)
  {
    VTK::FileType fileType = VTK::unstructuredGrid;
    VTK::OutputType outputtype = VTK::ascii;
    VTK::VTUWriter writer(s, outputtype, fileType);

    writer.beginMain(numParamEdges, numNodes);

    // Write nodes types
    writeGraphNodeTypes(writer);
    // Points
    writeGraphGridPoints(writer);
    // Cells
    writeGraphGridCells(writer);

    writer.endMain();
  }

  template<class ctype,int dim>
  void psurface::VTKIO<ctype,dim>::writeGraphNodeTypes(VTK::VTUWriter& writer)
  {
      std::string scalars = "nodetype";
      std::string vectors = "";

      writer.beginPointData(scalars, vectors);
      {
            std::tr1::shared_ptr<VTK::DataArrayWriter<ctype> > p
            (writer.makeArrayWriter<ctype>(scalars, 1, numNodes));
            for (int i = 0; i < numNodes; i++)
                p->write(nodeType[i]);
      } // p needs to go out of scope before we call endPointData()

      writer.endPointData();
  }

  // write the positions of vertices
  template<class ctype,int dim>
  void psurface::VTKIO<ctype,dim>::writeElementGridPoints(VTK::VTUWriter& writer)
  {
      writer.beginPoints();
      {
            std::tr1::shared_ptr<VTK::DataArrayWriter<ctype> > p
            (writer.makeArrayWriter<ctype>("Coordinates", 3, numVertices));
            if(!p->writeIsNoop()) {
                  for(int i = 0; i < numVertices; i++)
                        for(int l = 0; l < 3; l++)
                            p->write((par->vertices(i))[l]);
            }
      }
      writer.endPoints();
  }

  template<class ctype,int dim>
  void psurface::VTKIO<ctype,dim>::writeGraphGridPoints(VTK::VTUWriter& writer)
  {
      writer.beginPoints();
      {
            std::tr1::shared_ptr<VTK::DataArrayWriter<ctype> > p
            (writer.makeArrayWriter<ctype>("Coordinates", 3, numNodes));
            if(!p->writeIsNoop()) {
              for(int i = 0; i < numNodes; i++)
                for(int l = 0; l < 3; l++)
                  p->write((nodePositions[i])[l]);
            }
      }
      writer.endPoints();
  }

  // write the connectivity array
  template<class ctype,int dim>
  void psurface::VTKIO<ctype,dim>::writeElementGridCells(VTK::VTUWriter& writer)
  {
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
          }
      }

      // types
      {
          std::tr1::shared_ptr<VTK::DataArrayWriter<unsigned char> > p3
          (writer.makeArrayWriter<unsigned char>("types", 1, numTriangles));
          if(!p3->writeIsNoop())
          {
              for(int i = 0; i < numTriangles;i++)
              p3->write(5); //vtktype of triangle
          }
      }

      writer.endCells();
  }

  template<class ctype,int dim>
  void psurface::VTKIO<ctype,dim>::writeGraphGridCells(VTK::VTUWriter& writer)
  {
      writer.beginCells();
      // connectivity
      {
          std::tr1::shared_ptr<VTK::DataArrayWriter<int> > p1
          (writer.makeArrayWriter<int>("connectivity", 1, 2*numParamEdges));
          if(!p1->writeIsNoop())
          {
            for (size_t i = 0; i < parameterEdgeArray.size(); i++)
              for(int l = 0; l < 2; l++)
                p1->write((parameterEdgeArray[i])[l]);
          }
      }

      // offsets
      {
          std::tr1::shared_ptr<VTK::DataArrayWriter<int> > p2
          (writer.makeArrayWriter<int>("offsets", 1, numParamEdges));
          if(!p2->writeIsNoop()) {
              int offset = 0;
              for (size_t i = 0; i < parameterEdgeArray.size(); i++)
                {
                  offset += 2;
                  p2->write(offset);
                }
          }
      }

      // types
      {
          std::tr1::shared_ptr<VTK::DataArrayWriter<unsigned char> > p3
          (writer.makeArrayWriter<unsigned char>("types", 1, numParamEdges));
          if(!p3->writeIsNoop())
          {
            for (size_t i = 0; i < parameterEdgeArray.size(); i++)
              p3->write(3);//vtktype of edges
          }
      }

      writer.endCells();
  }

  template<class ctype,int dim>
  void psurface::VTKIO<ctype,dim>::writeElementGridCellData(VTK::VTUWriter& writer)
  {
    writer.beginCellData();

    // patch numbers
    {
      std::tr1::shared_ptr<VTK::DataArrayWriter<int> > p
        (writer.makeArrayWriter<int>("Patch", 1, numTriangles));
      if(!p->writeIsNoop()) {
        for(int i = 0; i < numTriangles; i++)
          p->write(par->triangles(i).patch);
      }
    }

    writer.endCellData();
  }


//   Explicit template instantiations.
namespace psurface {
  template class VTKIO<float,2>;
}

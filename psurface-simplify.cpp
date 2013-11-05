#include "config.h"

#include <iostream>
#include <sstream>
#include <fstream>

#include <stdexcept>
#include <memory>

#include <algorithm>
#include <list>

#include <vector>
#include <string>

#include <getopt.h>

#include "TargetSurface.h"
#include "AmiraMeshIO.h"
#include "PSurface.h"
#include "Hdf5IO.h"
#include "GmshIO.h"
#include "VtkIO.h"

#if defined HAVE_AMIRAMESH
#include <amiramesh/AmiraMesh.h>
#endif

#include "MultiDimOctree.h"
#include "EdgeIntersectionFunctor.h"
#include "QualityRequest.h"
#include "HxParamToolBox.h"

#include "VertexHeap.h"
#include "Triangulator.h"


using namespace std;
using namespace psurface;


/** \brief List of the file types that we support */
enum FileType
{
  AMIRA,
  HDF5,
  VTU,
  GMSH
};

////////////////////////////////////////////////////////////////////////////////
//// Routines for printing instructions and file handling
////////////////////////////////////////////////////////////////////////////////

void print_usage() {
  QualityRequest req;

  cerr << "Usage:" << endl
       << "   psurface-simplify -i <inputfilename> -o <outputfilename> (-n <node number> | -c <number of nodes>)" << endl
       << endl
       << "where " << endl
       << "-n removes a specific node with given node number (quality request will be ignored)" << endl
       << "-c removes given number of points" << endl
       << endl
       << "Optional arguments:" << endl
       << "-t x : set dihedral angle threshold to x         (default: " << req.dihedralAngleThreshold << ")" << endl
       << "-l x : set allowed max edge length to x          (default: " << req.maxEdgeLength          << ")" << endl
       << "-r x : set importance of aspect ratio to x       (default: " << req.aspectRatio            << ")" << endl
       << "-d x : set importance of Hausdorff distance to x (default: " << req.hausdorffDistance      << ")" << endl
       << "-b x : 1 to output just the basegrid             (default: " << 1                          << ")" << endl
       << "-s   : do not allow self intersection            (default: " << req.intersections          << ")" << endl
       << endl;
}

bool hasExtension(const string& input, const string& ext) {
  return input.find(ext) == input.length() - ext.length();
}

FileType filetypeOf(const string& filename) {
  FileType type;

  if(hasExtension(filename, ".am") || hasExtension(filename,".par"))
    type = AMIRA;
  else if(hasExtension(filename,".h5"))
    type = HDF5;
  else if(hasExtension(filename,".vtu"))
    type = VTU;
  else if(hasExtension(filename,".msh"))
    type = GMSH;
  else
    throw runtime_error(string("File ") + filename + " has unknown extension.");

  return type;
}


////////////////////////////////////////////////////////////////////////////////
//// Routines for removing multiple points according to a QualityRequest.
////////////////////////////////////////////////////////////////////////////////

void calcError(int vertex, const QualityRequest& quality, VertexHeap::ErrorValue& error,
               PSurface<2, float>* par, MultiDimOctree<Edge, EdgeIntersectionFunctor, float, 3>& edgetree) {
  int featureEdgeA, featureEdgeB;

  std::vector<std::vector<int> > halfStarVertices;
  std::vector<int>               fullStarVertices;
  std::vector<std::vector<int> > halfStarTris;
  std::vector<int>               fullStarTris;


  const int featureStatus = psurface::ParamToolBox::computeFeatureStatus(par, vertex, featureEdgeA, featureEdgeB);

  // Remove point according to feature status.
  if (psurface::ParamToolBox::FEATURE_POINT == featureStatus)
    error.block();
  else if (psurface::ParamToolBox::REGULAR_POINT == featureStatus) {
    error.unblock();

    // Any two edges will do here:
    featureEdgeA = par->vertices(vertex).edges[0];
    featureEdgeB = par->vertices(vertex).edges[1];

    // Finds the two halfstars that make up the full star.
    std::vector<int> patches;
    if (!psurface::ParamToolBox::findAllHalfStars(vertex, featureEdgeA, featureEdgeB, halfStarVertices, halfStarTris, patches, par)) {
      error.block();
      return;
    }

    psurface::ParamToolBox::makeFullStarOutOfHalfStars(halfStarVertices[0], halfStarTris[0],
                                                       halfStarVertices[1], halfStarTris[1],
                                                       fullStarVertices, fullStarTris);

    // Simulate the retriangulation to obtain its error.
    Triangulator::estimateStarError(fullStarVertices, vertex, quality, fullStarTris, error, edgetree, par);

    // Error is counter per halfstar.
    error.value /= 2;
  } else {
    error.unblock();
    error.value = 0;

    std::vector<int> patches;
    if (!psurface::ParamToolBox::findAllHalfStars(vertex, featureEdgeA, featureEdgeB,
                                                  halfStarVertices, halfStarTris, patches, par)) {
      error.block();
      return;
    }

    // Simulate the retriangulation to obtain its error.
    for (int i = 0; i < halfStarVertices.size(); ++i){
      if (halfStarTris[i].size()>1){
        VertexHeap::ErrorValue qualityValue;

        psurface::Triangulator::estimateHalfStarError(halfStarVertices[i], vertex,
                                                      quality, halfStarTris[i], qualityValue,
                                                      edgetree, par);

        if (qualityValue.isBlocked()){
          error.block();
          break;
        }

        error.value += qualityValue.value;
      }
    }

    // Error is counted per halfstar.
    error.value /= halfStarVertices.size();
  }
}

// Returns number of points removed.
int removePoint(int vertex, const psurface::QualityRequest& quality,
                PSurface<2, float>* par, MultiDimOctree<Edge, EdgeIntersectionFunctor, float, 3>* pedgetree) {
    int featureEdgeA, featureEdgeB;

    const int featureStatus = psurface::ParamToolBox::computeFeatureStatus(par, vertex, featureEdgeA, featureEdgeB);

    // Remove point according to feature status.
    if (psurface::ParamToolBox::FEATURE_POINT == featureStatus)
      return 0;
    else if (psurface::ParamToolBox::REGULAR_POINT == featureStatus) {
      if (!psurface::ParamToolBox::removeRegularPoint(par, vertex, quality, pedgetree))
        return 0;
    } else {
      if (!psurface::ParamToolBox::removeFeatureLinePoint(par, vertex, quality, featureStatus, featureEdgeA, featureEdgeB, pedgetree))
        return 0;
    }

    return 1;
}

void updateErrors(int vertex, const psurface::QualityRequest& quality,
                  PSurface<2, float>* par, MultiDimOctree<Edge, EdgeIntersectionFunctor, float, 3>& edgetree, VertexHeap& vertexHeap) {
    vector<int> neighbors = par->getNeighbors(vertex);

    for (int k = 0; k < neighbors.size(); ++k) {
      VertexHeap::ErrorValue error = vertexHeap.getError(neighbors[k]);

      calcError(neighbors[k], quality, error, par, edgetree);
      vertexHeap.reposition(neighbors[k], error);
    }
}

// Returns number of points removed.
int removeNumberOfPoints (int n, psurface::QualityRequest& req, PSurface<2, float>* par) {
  //// Setup certain objects

  // Setup quality request.
  req.normalize();

  // Remove triangular closure on each triangle.
  par->removeExtraEdges();

  // Setup octree for interesection tests.
  EdgeIntersectionFunctor ef(&(par->vertices(0)));
  MultiDimOctree<Edge, EdgeIntersectionFunctor, float, 3> edgetree;

  // Actually fill the tree with data only if we check for intersections.
  if (req.intersections) {
    Box<float, 3> box;
    par->getBoundingBox(box);

    // Careful: Storing a POINTER to the EdgeIntersectionIterator here !
    edgetree.init(box, &ef);

    // Simply insert every edge without specific order.
    for (int k = 0; k < par->getNumEdges(); ++k)
      edgetree.insert(&(par->edges(k)));
  }


  // Calculate the error that the removal of a certain point would introduce according to QualityRequest and save them in a heap.
  typedef vector<VertexHeap::ErrorValue> ErrorContainer;

  const int numVertices = par->getNumVertices();
  ErrorContainer error(numVertices);
  VertexHeap vertexHeap;

  for (int k = 0; k < numVertices; ++k)
    calcError(k, req, error[k], par, edgetree);

  vertexHeap.buildHeap(error);
  error.resize(0);

  //// Finally actually remove the points.
  int removedPoints = 0;
  while (removedPoints < n) {
    // Check whether there are still points available for removal.
    if ((-1 == vertexHeap.getMin()) or vertexHeap.isBlockedMin()) {
      cerr << "Could not find another point to remove." << endl;
      break;
    }

    // Get the index of the vertex to remove next.
    int index = vertexHeap.extractMin();

    // Now really remove a point.
    if (removePoint(index, req, par, &edgetree)) {
      updateErrors(index, req, par, edgetree, vertexHeap);
      ++removedPoints;
    } else {
      VertexHeap::ErrorValue oldErr = vertexHeap.getMinErrorStatus();

      oldErr.block();
      vertexHeap.insert(index, oldErr);
    }
  }

  //// Tidy up.
  if (req.intersections)
    edgetree.clear();

  par->garbageCollection();

  if (removedPoints > 0) {
    par->hasUpToDatePointLocationStructure = false;
    par->createPointLocationStructure();
  }

  return removedPoints;
}


////////////////////////////////////////////////////////////////////////////////
//// Main
////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) try {
  ////// Parse arguments.
  string input, output;
  int base = 1;
  QualityRequest req;

  bool nodeCount = false, nodeNumber = false;
  int n;

  int opt;

  while ((opt = getopt(argc, argv, ":i:o:n:b:t:l:r:d:s:c:")) != EOF) {
    switch (opt) {
    case 'i':
      input = optarg;
      break;
    case 'o':
      output = optarg;
      break;
    case 'n':
      if (nodeCount or nodeNumber)
        throw runtime_error("Specified already a node or a number of nodes to be removed.");

      stringstream(optarg) >> n;
      nodeNumber = true;
      break;
    case 'c':
      if (nodeCount or nodeNumber)
        throw runtime_error("Specified already a node or a number of nodes to be removed.");

      stringstream(optarg) >> n;
      nodeCount = true;
      break;
    case 'b':
      stringstream(optarg) >> base;
      break;
    case 't':
      req.smallDihedralAngles = true;
      stringstream(optarg) >> req.dihedralAngleThreshold;
      break;
    case 'l':
      stringstream(optarg) >> req.maxEdgeLength;
      break;
    case 'r':
      stringstream(optarg) >> req.aspectRatio;
      break;
    case 'd':
      stringstream(optarg) >> req.hausdorffDistance;
      break;
    case 's':
      req.intersections = true;
      break;
    default:
      print_usage();
      throw runtime_error("Tried to set invalid flag.");
    }
  }

  ////// Check arguments.
  // Got input and output filenames ?
  if (input.empty() or output.empty()) {
    print_usage();
    throw runtime_error("Input or output file not specified.");
  }

  // Got a node argument ?
  if (false == nodeNumber and false == nodeCount) {
    print_usage();
    throw runtime_error("Specified neither a node nor a number of nodes to be removed.");
  }

  // Check Filetype.
  FileType inputType = filetypeOf(input),
    outputType = filetypeOf(output);

  ////// Read input file.
  auto_ptr<PSurface<2,float> > par(new PSurface<2,float>);
  auto_ptr<Surface> surf(new Surface);

  switch(inputType) {
  case HDF5:
    {
#if HAVE_HDF5
      auto_ptr<Hdf5IO<float,2> > pconvert(new Hdf5IO<float,2>(par.get()));
      pconvert->initCompletePSurface(surf.get(), input);
#else
      std::cerr << "You have given an hdf5 input file, but psurface-simplify" << std::endl;
      std::cerr << "has been compiled without hdf5 support!" << std::endl;
      throw runtime_error("No hdf5 support.");
#endif
    }
    break;

  case GMSH:
    {
      par = auto_ptr<PSurface<2, float> >(GmshIO<float,2>::readGmsh(input));
    }
    break;

  case AMIRA:
    {
#if defined HAVE_AMIRAMESH
      AmiraMesh* am = AmiraMesh::read(input.c_str());
      PSURFACE_API AmiraMeshIO<float> amIO;
      if( !amIO.initFromAmiraMesh(par.get(), am, input.c_str(), surf.get()))
        throw runtime_error("Unable to initiate psurface from amiramesh file!");
#else
      std::cerr << "You have given an amira input file, but psurface-simplify" << std::endl;
      std::cerr << "has been compiled without amira support!" << std::endl;
      throw runtime_error("No amira support.");
#endif

    }
    break;
  default:
    {
      throw runtime_error("Unknown input type.");
    }
  };


  ////// Remove points.
  int ret;

  if (true == nodeCount)
    ret = removeNumberOfPoints(n, req, par.get());
  else // true == nodeNumber
    ret = removePoint(n, req, par.get(), NULL);

  // Print number of nodes removed.
  cout << ret << endl;


  ////// Write output file.
  switch(outputType) {
  case HDF5:
    {
#if HAVE_HDF5
      string xdmffile(output);
      xdmffile.erase (xdmffile.end() - 3, xdmffile.end());
      xdmffile.append(".xdmf");
      auto_ptr<Hdf5IO<float,2> > pn(new Hdf5IO<float,2>(par.get()));
      pn->createHdfAndXdmf(xdmffile, output, base);
#else
      cerr << "You have given an hdf5 output file, but psurface-simplify" << endl;
      cerr << "has been compiled without hdf5 support!" << endl;
      throw runtime_error("No hdf5 support.");
#endif
    }
    break;

  case VTU:
    {
      auto_ptr<VTKIO<float,2> > pn(new VTKIO<float,2>(par.get()));
      pn->createVTU(output.c_str(), base);
    }
    break;

  case AMIRA:
    {
#if defined HAVE_AMIRAMESH
      PSURFACE_API AmiraMeshIO<float> amIO;
      amIO.writeAmiraMesh(par.get(), output.c_str());
#else
      std::cerr << "You have given an amira output file, but psurface-simplify" << std::endl;
      std::cerr << "has been compiled without amira support!" << std::endl;
      throw runtime_error("No amira support.");
#endif
    }
    break;
  default:
    throw runtime_error("Unknown output type.");
  };

  return 0;
 } catch (const exception& e) {
  cerr << "ERROR: " << e.what() << endl;

  return -1;
 }

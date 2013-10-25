#include "config.h"

#include <iostream>
#include <sstream>
#include <fstream>

#include <stdexcept>
#include <memory>

#include <vector>
#include <string>

#include <getopt.h>

#ifdef PSURFACE_STANDALONE
#include "TargetSurface.h"
#else
#include "hxsurface/Surface.h"
#endif

#include "AmiraMeshIO.h"
#include "psurfaceAPI.h"
#include "PSurface.h"
#include "Hdf5IO.h"
#include "GmshIO.h"
#include "VtkIO.h"

#if defined HAVE_AMIRAMESH || !defined PSURFACE_STANDALONE
#include <amiramesh/AmiraMesh.h>
#endif

#include "MultiDimOctree.h"
#include "EdgeIntersectionFunctor.h"
#include "QualityRequest.h"
#include "HxParamToolBox.h"

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

void print_usage() {
  cerr << "Usage: psurface-simplify -i <inputfilename> -o <outputfilename> -n <nodenumber> -b <1/0>" << endl;
}

// The following two functions should actually be part of a class like "MeshFile".
bool hasExtension(const string& input, const string& ext) {
  return input.find(ext) == input.length() - ext.length();
}

// Could be solved more elegantly.
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

int main(int argc, char **argv) try {
  //// Parse arguments.
  if (argc < 8) {
    print_usage();
    return 1;
  }

  string input, output;
  int base = true;
  int index;

  int opt;

  // Invokes member function `int operator ()(void);'
  while ((opt = getopt(argc, argv, ":i:o:n:b:")) != EOF) {
    switch (opt) {
    case 'i':
      input = optarg;
      break;
    case 'o':
      output = optarg;
      break;
    case 'n':
      stringstream(optarg) >> index;
      break;
    case 'b':
      stringstream(optarg) >> base;
      break;
    default:
      print_usage();
      return 1;
    }
  }


  //// Check Filetype.
  FileType inputType = filetypeOf(input),
    outputType = filetypeOf(output);

  auto_ptr<PSurface<2,float> > par(new PSurface<2,float>);
  auto_ptr<Surface> surf(new Surface);
  //par->surface = new Surface;

  //// Read input file.
  switch(inputType) {
  case HDF5:
    {
#if HAVE_HDF5
      auto_ptr<Hdf5IO<float,2> > pconvert(new Hdf5IO<float,2>(par.get()));
      pconvert->initCompletePSurface(surf.get(), input);
#else
      std::cerr << "You have given an hdf5 input file, but psurface-simplify" << std::endl;
      std::cerr << " has been compiled without hdf5 support!" << std::endl;
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
      AmiraMesh* am = AmiraMesh::read(input.c_str());
      PSURFACE_API AmiraMeshIO<float> amIO;
      if( !amIO.initFromAmiraMesh(par.get(), am, input.c_str(), surf.get()))
    throw runtime_error("Unable to initiate psurface from amiramesh file!");
    }
    break;
  default:
    {
      throw runtime_error("Unkown input type.");
    }
  };

  //// Remove point.
  cout << "Removing node " << index << "." << endl;

  Box<float, 3> box;
  par->getBoundingBox(box);
  EdgeIntersectionFunctor ef(&(par->vertices(0)));
  MultiDimOctree<Edge, EdgeIntersectionFunctor, float, 3> edgebox(box, &ef);

  QualityRequest req;
  //req.intersections = true;
  //req.smallDihedralAngles = true;
  //req.paths = false;

  ParamToolBox::removeRegularPoint(par.get(), index, req, NULL); //&edgebox);
  par->garbageCollection();

  cout << "Node removed." << endl;

  //// Write output file.
  switch(outputType) {
  case HDF5:
    {
#if HAVE_HDF5
      string xdmffile(output);
      xdmffile.erase (xdmffile.end() - 3, xdmffile.end());
      xdmffile.append(".xdmf");
      auto_ptr<Hdf5IO<float,2> > pn(new Hdf5IO<float,2>(par.get()));
      pn->creatHdfAndXdmf(xdmffile, output, base);
#else
        cerr << "You have given an hdf5 input file, but psurface-simplify" << endl;
        cerr << " has been compiled without hdf5 support!" << endl;
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
      PSURFACE_API AmiraMeshIO<float> amIO;
      amIO.writeAmiraMesh(par.get(), output.c_str());
    }
    break;
  default:
    throw runtime_error("Unknown output type.");
  };

  return 0;
 } catch (const exception& e) {
  cout << e.what() << endl;
 }

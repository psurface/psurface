#include "config.h"

#include <iostream>
#include <sstream>

#include <memory>
#include <stdexcept>

#include <vector>
#include <string>

#include <getopt.h>

#include "TargetSurface.h"
#include "AmiraMeshIO.h"
#include "PSurface.h"
#include "Hdf5IO.h"
#include "VtkIO.h"

#if defined HAVE_AMIRAMESH
#include <amiramesh/AmiraMesh.h>
#endif

#include "PSurfaceSmoother.h"


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
  cerr << "Usage:" << endl
       << "   psurface-smooth -i <inputfilename> -o <outputfilename> -n <number of smoothing cycles>" << endl
       << "                   -k (optional) set if patches need NOT to be preserved (default: not set)"
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
    throw runtime_error("Gmsh files are not supported here.");
  else
    throw runtime_error(string("File ") + filename + " has unknown extension.");

  return type;
}



////////////////////////////////////////////////////////////////////////////////
//// Main
////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) try {
  ////// Parse arguments.
  string input, output;
  int n = -1;
  bool keepPatches = true;

  int opt;

  while ((opt = getopt(argc, argv, ":i:o:n:k")) != EOF) {
    switch (opt) {
    case 'i':
      input = optarg;
      break;
    case 'o':
      output = optarg;
      break;
    case 'n':
      stringstream(optarg) >> n;
      break;
    case 'k':
      keepPatches = false;
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

  // Got a valid number of smoothing cycles ?
  if (0 > n) {
    print_usage();
    throw runtime_error("Number of smoothing cycles not set or negative.");
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
      cerr << "You have given an hdf5 input file, but psurface-simplify" << endl;
      cerr << "has been compiled without hdf5 support!" << endl;
      throw runtime_error("No hdf5 support.");
#endif
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
      cerr << "You have given an amira input file, but psurface-simplify" << endl;
      cerr << "has been compiled without amira support!" << endl;
      throw runtime_error("No amira support.");
#endif

    }
    break;
  default:
    {
      throw runtime_error("Unknown input type.");
    }
  };


  ////// Smooth.
  std::vector<unsigned int> nodeStack;

  for (int i = 0; i < n; ++i)
    for (int k = 0; k < par->getNumEdges(); ++k)
      PSurfaceSmoother<float>::applyEdgeRelaxation(par.get(), k, keepPatches, nodeStack);


  ////// Write output file.
  switch(outputType) {
  case HDF5:
    {
#if HAVE_HDF5
      string xdmffile(output);
      xdmffile.erase (xdmffile.end() - 3, xdmffile.end());
      xdmffile.append(".xdmf");
      auto_ptr<Hdf5IO<float,2> > pn(new Hdf5IO<float,2>(par.get()));
      pn->createHdfAndXdmf(xdmffile, output, false);
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
      pn->createVTU(output.c_str(), output.substr(0, output.length()-string(".vtu").length()) + "-graph.vtu");
    }
    break;

  case AMIRA:
    {
#if defined HAVE_AMIRAMESH
      PSURFACE_API AmiraMeshIO<float> amIO;
      amIO.writeAmiraMesh(par.get(), output.c_str());
#else
      cerr << "You have given an amira output file, but psurface-simplify" << endl;
      cerr << "has been compiled without amira support!" << endl;
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

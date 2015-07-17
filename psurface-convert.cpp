#include "config.h"

#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <iostream>
#include <string.h>
#include <fstream>
#include <memory>

#ifdef PSURFACE_STANDALONE
#include "TargetSurface.h"
#else
#include "hxsurface/Surface.h"
#endif

#include "AmiraMeshIO.h"
#include "PSurface.h"
#include "Hdf5IO.h"
#include "GmshIO.h"
#include "VtkIO.h"

#if defined HAVE_AMIRAMESH || !defined PSURFACE_STANDALONE
#include <amiramesh/AmiraMesh.h>
#endif

/** \brief List of the file types that we support */
enum FileTypes
{
  AMIRA,
  HDF5,
  VTU,
  GMSH
};


using namespace psurface;
int main(int argc, char **argv)
{
    if (argc < 4) {
      fprintf(stderr, "Usage: psurface_convert -i inputname -o outputname (-t type) \n");
      fprintf(stderr, "Input file type could be amiramesh(*.am) , hdf5(*.h5) or gmsh(*.msh).\n");
      fprintf(stderr, "Output file type could be amiramesh(*.am) , hdf5(*.h5) or vtu(*.vtu).\n");
      fprintf(stderr, "type could be b(basegrid) or r(readable hdf5 file).\n -t b means that the output should only have base grid trianlge(This option is used when the output type is vtu type.\n -t r means that we get readable output hdf5 type data(This option is used when the output type is hdf5).\n");
      exit(0);
    }

    //use get opt to deal with the argv
    char *input, *output, *type;
    bool basegrid = 0, basehdf5 = 1;
    int opt=0;
    int i=0;
    const char* optstring=":i:o:t:";
    const int num=3;

    while((opt=getopt(argc,argv,optstring)) != -1)
    {
        switch(opt)
        {
        case 'i':
             input = optarg;
             break;
        case 'o':
             output = optarg;
             break;
        case 't':
             type = optarg;
             break;
        case ':':
            printf("the option needs a value\n");
            break;
        case '?':
            printf("unknown option '%c'\n",optopt);
            break;
        }
    }

    for(i=0;optind<argc;i++,optind++)
    {
        if(i<num)
            printf("argument:%s/n",argv[optind]);
        else
            printf("excess argument:%s/n",argv[optind]);
    }

    FileTypes inputType, outputType;
    if(strstr(input,".am") != NULL || strstr(input,".par") != NULL )
        inputType = AMIRA;
    else if(strstr(input,".h5") != NULL)
        inputType = HDF5;
    else if(strstr(input,".vtu") != NULL)
        inputType = VTU;
    else if(strstr(input,".msh") != NULL)
        inputType = GMSH;
    else
      printf(" could not tell the input type by file extension\n");

    if(strstr(output,".am") != NULL)
        outputType = AMIRA;
    else if(strstr(output,".h5") != NULL)
    {
        outputType = HDF5;
        if( type != NULL && *type == 'r')  basehdf5 = 0;
    }
    else if(strstr(output,".vtu") != NULL)
    {
        outputType = VTU;
        if( type != NULL && *type == 'b')  basegrid = 1;
    }
    else if(strstr(output,".msh") != NULL)
        outputType = GMSH;
    else
      printf(" could not tell the output type by file extension\n");

  PSurface<2,float>* par;

  switch(inputType)
  {
    case HDF5:
      {
#if HAVE_HDF5
        par = Hdf5IO<float,2>::read(input);
#else
        std::cerr << "You have given an hdf5 input file, but psurface-convert" << std::endl;
        std::cerr << "has been compiled without hdf5 support!" << std::endl;
        exit(1);
#endif
      }
      break;

    case GMSH:
      {
        par = GmshIO<float,2>::readGmsh(input);
      }
      break;

    case AMIRA:
      {
#if HAVE_AMIRAMESH
        AmiraMesh* am = AmiraMesh::read(input);
        AmiraMeshIO<float> amIO;
        Surface* surf = new Surface;
        par = new PSurface<2,float>;

        if( !amIO.initFromAmiraMesh(par,am,input, surf))
        {
          printf("unable to initiate psurface from amira mesh file!\n");
          return 0;
        }
#else
        std::cerr << "You have given an AmiraMesh input file, but psurface-convert" << std::endl;
        std::cerr << "has been compiled without AmiraMesh support!" << std::endl;
        exit(1);
#endif
      }
      break;
    default:
    {
//     throw(std::runtime_error("unkown input type\n"));
     return 0;
    }
  };

  std::string str(output);
  switch(outputType)
   {
    case HDF5:
      {
#if HAVE_HDF5
        std::string xdmffile(output);
        xdmffile.erase (xdmffile.end() - 3, xdmffile.end());
        xdmffile.append(".xdmf");
        Hdf5IO<float,2>* pn = new Hdf5IO<float,2>(par);
        pn->createHdfAndXdmf(xdmffile, output,basehdf5);
#else
        std::cerr << "You have given an hdf5 output file, but psurface-convert" << std::endl;
        std::cerr << "has been compiled without hdf5 support!" << std::endl;
        exit(1);
#endif
      }
      break;

    case VTU:
      {
        VTKIO<float,2> pn(par);
        pn.createVTU(str.c_str(), basegrid ? "" : (str.substr(0, str.length()-std::string(".vtu").length())) + "-graph.vtu");
      }
      break;

    case AMIRA:
      {
#if HAVE_AMIRAMESH
          PSURFACE_API AmiraMeshIO<float> amIO;
          amIO.writeAmiraMesh(par, output);
#else
        std::cerr << "You have given an AmiraMesh output file, but psurface-convert" << std::endl;
        std::cerr << "has been compiled without AmiraMesh support!" << std::endl;
        exit(1);
#endif
      }
      break;
    default:
       printf("unknown output type\n");
//     throw(std::runtime_error("unkown output type\n"));
   };

  return 0;
}


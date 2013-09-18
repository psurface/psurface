#define  HAVE_AMIRAMESH
#include "PSurface.h"
#include "psurfaceAPI.h"
#ifdef PSURFACE_STANDALONE
#include "TargetSurface.h"
#endif
#if defined HAVE_AMIRAMESH || !defined PSURFACE_STANDALONE
#include <amiramesh/AmiraMesh.h>
#endif
#include "AmiraMeshIO.h"
#include "psurface_convert.h"
#include <iostream>
#include <stdio.h>  
#include <unistd.h>
#include <string.h>
namespace psurface
{
  enum FileTypes{
  AMIRA,
  HDF5,
  VTU,
  GMSH};
}
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
    bool basegrid = 0, readablehdf = 0;
    int opt=0;  
    int i=0;
    const char *optstring=":i:o:t:";  
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
            printf("the option need a value\n");  
            break;  
        case '?':  
            printf("unknow option%c\n",optopt);  
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
      printf(" could not tell the input type by file extention\n");
      
    if(strstr(output,".am") != NULL)
        outputType = AMIRA;
    else if(strstr(output,".h5") != NULL)
    {
        outputType = HDF5; 
        if( type != NULL && *type == 'r')  readablehdf = 1;
    }
    else if(strstr(output,".vtu") != NULL)
    {
        outputType = VTU;
        if( type != NULL && *type == 'b')  basegrid = 1;
    }
    else if(strstr(output,".msh") != NULL)
        outputType = GMSH;
    else
      printf(" could not tell the output type by file extention\n");
  
  PSurface<2,float>* par = new PSurface<2,float>;
  Surface* surf = new Surface;

  switch(inputType)
  {
    case HDF5:
      {
        PsurfaceConvert<float,2>* pconvert = new PsurfaceConvert<float,2>(input, 1);
        if(!pconvert->initPsurface(par, surf, 0))
        {
          printf("unable to initiate psurface from hdf5 file!\n");
          return 0;
        }
      }
      break;

    case GMSH:
      {
        PsurfaceConvert<float,2>* pconvert = new PsurfaceConvert<float,2>(input, 0);      
        if( !pconvert->initPsurface(par, surf, 1))
        {
          printf("unable to initiate psurface from GMSH file!\n");
          return 0;
        }

      }
      break;

    case AMIRA:
      {
        AmiraMesh* am = AmiraMesh::read(input);
        PSURFACE_API AmiraMeshIO<float> amIO;    
        if( !amIO.initFromAmiraMesh(par,am,input, surf))
        {
          printf("unable to initiate psurface from amira mesh file!\n");
          return 0;
        }
      }
      break;
    default:
    {
     throw(std::runtime_error("unkown input type\n"));
     return 0;
    }
  };

   switch(outputType)
   {
    case HDF5:
      {
      char xdmffile[100];
      strcpy(xdmffile,output);
      strcat(xdmffile, ".xdmf");
      PsurfaceConvert<float,2>* pn = new PsurfaceConvert<float,2>(par);      
      pn->creatHdfAndXdmf(xdmffile, output,readablehdf);
      }
      break;

    case VTU:
      {
        PsurfaceConvert<float,2>* pn = new PsurfaceConvert<float,2>(par);
        pn->creatVTU(output,basegrid);
      }
      break;

    case AMIRA:
      {
      PSURFACE_API AmiraMeshIO<float> amIO;
      amIO.writeAmiraMesh(par, output);
      }
      break;
    default:
     throw(std::runtime_error("unkown output type\n"));
   };

  return 0;
}

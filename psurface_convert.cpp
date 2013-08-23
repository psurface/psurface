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
namespace psurface{
  enum FileTypes{
  AMIRA,
  HDF5,
  VTU,
  GMSH};
}
using namespace psurface;
int main(int argc, char **argv)
{
    if (argc < 2) {
	fprintf(stderr, "Usage: psurface_convert -i inputname -o outputname ...\n");
	fprintf(stderr, "Input file type could be amiramesh(*.am) , hdf5(*.h5) or gmsh(*.msh).\n");
	fprintf(stderr, "Output file type could be amiramesh(*.am) , hdf5(*.h5) or vtu(*.vtu).\n");
	exit(0);
    }
    //use get opt to deal with the argv
    char *input, *output;
    int opt=0;  
    int i=0;  
    const char *optstring=":i:o:";  
    const int num=2;  
  
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
        case ':':  
            printf("the option need a value/n");  
            break;  
        case '?':  
            printf("unknow option：%c/n",optopt);  
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
    
    printf("input = %s ouput = %s\n", input, output);

    FileTypes inputType, outputType;    
    if(strstr(input,".am") != NULL)
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
        outputType = HDF5;
    else if(strstr(output,".vtu") != NULL)
        outputType = VTU;
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
      PsurfaceConvert<float,2>* pconvert = new PsurfaceConvert<float,2>(input, 0);
      pconvert->initPsurface(par, surf, 0);
      }
      break;

    case GMSH:
      {
      PsurfaceConvert<float,2>* pconvert = new PsurfaceConvert<float,2>(input, 1);      
      pconvert->initPsurface(par, surf, 1);
      }
      break;

    case AMIRA:
      {
        AmiraMesh* am = AmiraMesh::read(input);
        PSURFACE_API AmiraMeshIO<float> amIO;    
        amIO.initFromAmiraMesh(par,am,input, surf);
      }
      break;
    default:
     throw(std::runtime_error("unkown input type\n"));
  };

   switch(outputType)
   {
    case HDF5:
      {
      char *xdmffile;
      strcpy (xdmffile,output);
      strcat(xdmffile, ".xdmf");
      PsurfaceConvert<float,2>* pn = new PsurfaceConvert<float,2>(par);      
      pn->creatHdfAndXdmf(xdmffile, output);
      }
      break;

    case GMSH:
      {
        PsurfaceConvert<float,2>* pn = new PsurfaceConvert<float,2>(par);
        pn->creatVTU(output);
      }
      break;

    case AMIRA:
      {
      PSURFACE_API AmiraMeshIO<float> amIO;    
      amIO.writeAmiraMesh(par, "output");
      }
      break;
    default:
     throw(std::runtime_error("unkown output type\n"));
   };

  return 0;
/*  //psurface_convert -i hdf5 -o amiramesh
  {
    PsurfaceConvert pconvert = new PsurfaceConvert<2,float>(inputfile, 0);
    initPsurface(par, surf, 0);
    PSURFACE_API AmiraMeshIO<float> amIO;    
    amIO.writeAmiraMesh(par1, "outputfile");
  }
  //psurface_convert -i hdf5 -o vtu
  {
    PsurfaceConvert pconvert = new PsurfaceConvert<2,float>(inputfile, 0);
    initPsurface(par, surf, 0);
    PSurfaceConvert pn = new PSurfaceConvert<2,float>(par);
    pn.creatVTU(outputfile);
  }
  //psurface_convert -i gmsh -o hdf5
  {
    PsurfaceConvert pconvert = new PsurfaceConvert<2,float>(inputfile, 1);
    initPsurface(par, surf, 1);
    PSurfaceConvert pn = new PSurfaceConvert<2,float>(par);
    pn.creatHdfAndXdmf(xdf_filename, hdf_filename);
  }
  //psurface_convert -i gmsh -o amiramsh
  {
    PsurfaceConvert pconvert = new PsurfaceConvert<2,float>(inputfile, 1);
    initPsurface(par, surf, 1);
    PSURFACE_API AmiraMeshIO<float> amIO;    
    amIO.writeAmiraMesh(par1, "outputfile");
  }
  //psurface_convert -i gmsh -o vtu
  {
    PsurfaceConvert pconvert = new PsurfaceConvert<2,float>(inputfile, 1);
    initPsurface(par, surf, 1);
    PSurfaceConvert pn = new PSurfaceConvert<2,float>(par);
    pn.creatVTU(outputfile);
  }

  //psurface_convert -i amiramesh -o hdf5
  {
    AmiraMesh* am = AmiraMesh::read(inputfile);
    if (!am.initFromAmiraMesh(par,am,inputfile, surf)) {
      printf("error in getting psurface from sphere.par.am");
    }
    else
    {
      PSurfaceConvert pn = new PSurfaceConvert<2,float>(par);
      pn.creatHdfAndXdmf(xdf_filename, hdf_filename);
    }     
  }
  //psurface_convert -i amiramesh -o vtu
  {
    AmiraMesh* am = AmiraMesh::read(inputfile);
    if (!am.initFromAmiraMesh(par,am,inputfile, surf)) {
      printf("error in getting psurface from sphere.par.am");
    }
    else
    {
      PSurfaceConvert pn = new PSurfaceConvert<2,float>(par);
      pn.creatVTU(outputfile);
    } 
  }*/
}

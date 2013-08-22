#define  HAVE_AMIRAMESH
#define  PSURFACE_STANDALONE

#include <vector>
#include "StaticVector.h"
#include "Domains.h"
#include "PSurface.h"
#include "PSurfaceFactory.h"
#include "psurfaceAPI.h"

#ifdef PSURFACE_STANDALONE
#include "TargetSurface.h"
#else
#include <hxsurface/Surface.h>
#endif

#if defined HAVE_AMIRAMESH || !defined PSURFACE_STANDALONE
#include <amiramesh/AmiraMesh.h>
#endif
#include "AmiraMeshIO.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <hdf5.h>
#include "hdf5IO.h"
using namespace psurface;

int main(int argc, char *argv[])
{
//    char* amfile="fibula_proximal_15mm.par";
//    char* amfile="femur_distal_15mm.par";
    char* amfile = "tibia_proximal_15mm.par";
    char* hdffile = "psurface.h5";
    char* outhdffile = "psurface2.h5";
    if(argc>1)
      amfile = argv[1];
    PSURFACE_API AmiraMeshIO<float> amIO;
    AmiraMesh* am = AmiraMesh::read(amfile);
    PSurface<2,float>* par1 = new PSurface<2,float>;
    Surface* surf1 = new Surface;
    if (!amIO.initFromAmiraMesh(par1,am,amfile, surf1)) {
    printf("error in getting psurface ");
    }
    else
    {
      Hdf5IO<float,2> hw;
      hw.par = par1;
      hw.update();
      hw.writeHdf5Data(hdffile);
      hw.writeXdmf("t1.xmf",hdffile);
    }

    amIO.writeAmiraMesh(par1, "amira.am");


    Hdf5IO<float,2> hx;
    PSurface<2,float>* par2 = new PSurface<2,float>;
    Surface* surf2 = new Surface;
    hx.initFromHDF5(par2,surf2, hdffile);

    PSURFACE_API AmiraMeshIO<float> amIO2;
    amIO2.writeAmiraMesh(par2, "hdf5.am");
    return 0; 
  }

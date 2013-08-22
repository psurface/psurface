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
#include <psurface/AmiraMeshIO.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gmshIO.h" 
using namespace psurface;

int main(int argc, char *argv[])
{
//    char* amfile="fibula_proximal_15mm.par";
//    char* amfile="femur_distal_15mm.par";
    char* amfile = "tibia_proximal_15mm.par";
//    char* amfile = "sphere.par.am";
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
      gmsh<float,2> gm;
      gm.par = par1;
      gm.update();
      gm.writeGmsh("test.msh");
    }
    return 0; 
  }

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
#include "hdf5IO.h"
using namespace psurface;

int main(int argc, char *argv[])
{
    //get psurface from msh file
    char* mshfile = "curved2d.msh";
    gmsh<float,2> gm;
    PSurface<2,float>* par = new PSurface<2,float>;
    Surface* surf = new Surface;
    gm.initFromGmsh(par, surf, mshfile);

    //write it to hdf and xdmf file
    char* hdffile = "psurface.h5";
    Hdf5IO<float,2> hw;
    hw.par = par;
    hw.update();
    hw.writeHdf5Data(hdffile);
    hw.writeXdmf("t1.xmf",hdffile);

    //write it to amiramesh
    PSURFACE_API AmiraMeshIO<float> amIO2;
    amIO2.writeAmiraMesh(par, "curved2d.am");

    return 0; 
  }

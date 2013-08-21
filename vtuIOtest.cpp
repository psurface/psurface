//the code to write psurface object into vtk file
//vtuwrite
#ifdef HAVE_CONFIG_H
#include "config.h"     
#endif
#define  HAVE_AMIRAMESH
#include <vector>
#include "StaticVector.h"
#include "Domains.h"
#include "PSurface.h"
#include "PSurfaceFactory.h"
#include "psurfaceAPI.h"
//#ifdef PSURFACE_STANDALONE
#include "TargetSurface.h"
//#endif
//#if defined HAVE_AMIRAMESH || !defined PSURFACE_STANDALONE
#include <amiramesh/AmiraMesh.h>
//#endif
//vtu file writen by christin
#include "dune/vtuwriter.hh"
#include "psurface/AmiraMeshIO.h"
#include <iostream>
#include "vtuIO.h"
using namespace Dune;
int main(int argc, char **argv)
{
  psurface::PSURFACE_API AmiraMeshIO<double> amIO,amIO1,amIO2;
  const char *file1, *file2,*file3;
  file1  = "cube.par.am";
  file2  = "sphere.par.am";
  file3  = "trice.par.am";
  AmiraMesh* am1 = AmiraMesh::read(file1);
  psurface::PSurface<2,double>* par1 = new psurface::PSurface<2,double>;
  psurface::Surface* surf1 = new psurface::Surface;
  if (!amIO.initFromAmiraMesh(par1,am1,file1, surf1)) {
  printf("error in getting psurface from cube.par.am");
  }
  else
  {
    VtkPWriter<double,2> vw;
    vw.update(par1);
    vw.writeUnstructureMesh("cube.vtu");
  }

  AmiraMesh* am2 = AmiraMesh::read(file2);
  psurface::PSurface<2,double>* par2 = new psurface::PSurface<2,double>;
  psurface::Surface* surf2 = new psurface::Surface;
  if (!amIO1.initFromAmiraMesh(par2,am2,file2, surf2)) {
  printf("error in getting psurface from sphere.par.am");
  }
  else
  {
    VtkPWriter<double,2> vw;
    vw.update(par1);
    vw.writeUnstructureMesh("sphere.vtu");
  }

  AmiraMesh* am3 = AmiraMesh::read(file3);
  psurface::PSurface<2,double>* par3 = new psurface::PSurface<2,double>;
  psurface::Surface* surf3 = new psurface::Surface;
  if (!amIO2.initFromAmiraMesh(par3,am3,file3, surf3)) {
  printf("error in getting psurface from trice.par.am");
  }
  else
  {
    VtkPWriter<double,2> vw;
    vw.update(par1);
    vw.writeUnstructureMesh("trice.vtu");
  }

  return 0;
}

#ifndef PSURFACE_CONVERT_NEW_H
#define PSURFACE_CONVERT_NEW_H
#include <memory>
#include <iostream>
#include <string.h>
#include <fstream>

#define  HAVE_AMIRAMESH
#ifdef PSURFACE_STANDALONE
#include "TargetSurface.h"
#else
#include "hxsurface/Surface.h"
#endif

#if defined HAVE_AMIRAMESH || !defined PSURFACE_STANDALONE
#include <amiramesh/AmiraMesh.h>
#endif
#include "AmiraMeshIO.h"

#ifdef PSURFACE_STANDALONE
namespace psurface { class Surface; }
#else
class Surface;
#endif

namespace psurface{
  
  enum NodeTypes{INTERIOR_NODE,
  INTERSECTION_NODE,
  CORNER_NODE,
  TOUCHING_NODE,
  GHOST_NODE};

  enum FileTypes{
  AMIRA,
  HDF5,
  VTU,
  GMSH};

}
#endif

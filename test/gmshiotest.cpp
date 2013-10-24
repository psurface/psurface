#include "config.h"

#include <exception>
#include <iostream>
#include <memory>
#include <string>
#include <vector>


#ifdef PSURFACE_STANDALONE
#include "TargetSurface.h"
#else
#include "hxsurface/Surface.h"
#endif

#include "PSurface.h"
#include "GmshIO.h"


using namespace std;
using namespace psurface;

typedef vector<string> StringVector;


int main(int argc, char* argv[])
{
  //// Meshes to test.
  const string basepath("examplefiles/");
  StringVector input(2);
  input[0] = "tricube-anticlockwise.msh";
  input[1] = "tricube-clockwise.msh";

  // Create new PSurface.
  auto_ptr<PSurface<2,float> > par(new PSurface<2,float>);
  par->surface = new Surface;

  // Read in meshes and check for consistency.
  try {
    for (StringVector::const_iterator it = input.begin(); it != input.end(); ++it) {
      // Read mesh.
      const string filename = basepath + *it;
      cout << "Testing using " << filename << endl;

      auto_ptr<GmshIO<float,2> > pIn(new GmshIO<float,2>(par.get()));
      pIn->readGmsh(par->surface, basepath + *it);

      // Check for consistency.
      par->checkConsistency((basepath + *it).c_str());
    }
  } catch (const exception& e) {
    cout << e.what() << endl;

    return 1;
  }

  return 0;
}

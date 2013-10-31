#include "config.h"

#include "fenv.h"

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

#include "MultiDimOctree.h"
#include "EdgeIntersectionFunctor.h"
#include "QualityRequest.h"
#include "HxParamToolBox.h"


using namespace std;
using namespace psurface;

typedef vector<string> StringVector;


int main(int argc, char* argv[])
{
  feenableexcept(FE_INVALID);

  //// Meshes to test.
  const string basepath("examplefiles/");
  StringVector input(2);
  input[0] = "tricube-anticlockwise.msh";
  input[1] = "tricube-clockwise.msh";

  // Read in a mesh, remove a node and check for consistency.
  try {
    for (StringVector::const_iterator it = input.begin(); it != input.end(); ++it) {
      const string filename = basepath + *it;

      cout << "Testing using " << filename << endl;

      // Remove the first eight nodes individually.
      for (size_t index = 0; index < 8; ++index) {
        // Read mesh.
        auto_ptr<PSurface<2,float> > par(GmshIO<float,2>::readGmsh(filename));

        // Remove node.
        cout << "   Removing node " << index << "." << endl;

        Box<float, 3> box;
        par->getBoundingBox(box);
        EdgeIntersectionFunctor ef(&(par->vertices(0)));
        MultiDimOctree<Edge, EdgeIntersectionFunctor, float, 3> edgebox(box, &ef);

        QualityRequest req;
        //req.intersections = true;
        //req.smallDihedralAngles = true;
        //req.paths = false;

        ParamToolBox::removeRegularPoint(par.get(), index, req, &edgebox);

        cout << "   Node removed." << endl;

        // Check consistency.
        par->checkConsistency(filename.c_str());
      }

      cout << "Done." << endl;
    }
  } catch (const exception& e) {
    cout << e.what() << endl;

    return 1;
  }

  return 0;
}

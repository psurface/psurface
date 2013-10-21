#ifndef NORMAL_PROJECTOR_H
#define NORMAL_PROJECTOR_H

#include "StaticVector.h"
#include "PSurface.h"
#include "PathVertex.h"

#include "psurfaceAPI.h"

#include <vector>

#ifdef PSURFACE_STANDALONE
namespace psurface { class Surface; }
#else
class Surface;
#endif

namespace psurface {

template <int dim, class ctype>
class PSurfaceFactory;
class NodeBundle;
class GlobalNodeIdx;
template <int dimworld, class ctype>
struct DirectionFunction;

/** \brief Construct a PSurface object by projecting one surface in normal direction onto another

\tparam ctype The type used for coordinates
*/
template <class ctype>
class PSURFACE_API NormalProjector {
public:

    NormalProjector(PSurface<2,ctype>* psurface)
        : psurface_(psurface)
    {}

    void project(const Surface* targetSurface,
                 const DirectionFunction<3,ctype>* domainDirection,
                 const DirectionFunction<3,ctype>* targetDirection
                    );

protected:

    void computeDiscreteDomainDirections(const DirectionFunction<3,ctype>* direction,
                                         std::vector<StaticVector<ctype,3> >& normals);

    void computeDiscreteTargetDirections(const Surface* targetSurface,
                                         const DirectionFunction<3,ctype>* direction,
                                         std::vector<StaticVector<ctype,3> >& normals);
                     
    void setupEdgePointArrays();

    /** \brief Insert a target edge using the vertices stored the edgePath vector. */
    void insertEdge(PSurfaceFactory<2,ctype>& factory,
                    int from, int to,
                    std::vector<PathVertex<ctype> >& edgePath);

    /** \brief Insert a segment on the edge path of a target edge. */
    void insertEdgeSegment(PSurfaceFactory<2,ctype>& factory,
                                int from, int to,
                                std::vector<PathVertex<ctype> >& edgePath);

    // ///////////////////////////////////////////////////////////////////////
    //   Methods needed to test whether an edge can be projected completely
    // ///////////////////////////////////////////////////////////////////////

    /** \brief Check if we can insert a target edge and store the path on the domain surface. */
    bool edgeCanBeInserted(const std::vector<StaticVector<ctype,3> >& normals,
                           int from, 
                           int to,
                           const std::vector<NodeBundle>& projectedTo,
                           std::vector<PathVertex<ctype> >& edgePath);

    bool testInsertEdgeFromInteriorNode(const std::vector<StaticVector<ctype,3> >& normals, 
                                        int from, int to, ctype& lambda,
                                        NodeBundle& curr,
                                        typename Node<ctype>::NodeType& currType, int& currTri,
                                        int& enteringEdge, std::vector<PathVertex<ctype> >& edgePath,
                                        int offset);

    bool testInsertEdgeFromIntersectionNode(const std::vector<StaticVector<ctype,3> >& normals, 
                                            int from, int to, ctype& lambda,
                                            NodeBundle& curr,
                                            typename Node<ctype>::NodeType& currType, int& currTri,
                                            int& enteringEdge,std::vector<PathVertex<ctype> >& edgePath,
                                            int offset);

    bool testInsertEdgeFromTouchingNode(const std::vector<StaticVector<ctype,3> >& normals, 
                                        int from, int to, ctype& lambda,
                                        NodeBundle& curr,
                                        typename Node<ctype>::NodeType& currType, int& currTri,
                                        int& enteringEdge, std::vector<PathVertex<ctype> >& edgePath,
                                        int offset);

    bool testInsertEdgeFromCornerNode(const std::vector<StaticVector<ctype,3> >& normals, 
                                      int from, int to, ctype& lambda,
                                      NodeBundle& curr, 
                                      typename Node<ctype>::NodeType& currType, int& currTri,
                                      int& enteringEdge, std::vector<PathVertex<ctype> >& edgePath,
                                      int offset);

    void addCornerNodeBundle(int v, 
                             int nN
                             );

    /** \brief Return true if the two NodeBundles share at least one triangle */
    bool onSameTriangle(const NodeBundle& a, const NodeBundle& others) const;

    /** \brief Return true if the NodeBundle contains a Node on the given triangle */
    bool onSameTriangle(const int& tri, const NodeBundle& others) const;

    int getCommonTri(const NodeBundle& a, const NodeBundle& b);

    std::vector<int> getCommonTris(const NodeBundle& a, const NodeBundle& b);

    NodeIdx getCornerNode(const DomainTriangle<ctype>& cT, int corner);

    /** \brief Get the type of the nodes in a bundle (it must be the same for all) */
    typename Node<ctype>::NodeType type(const NodeBundle& b) const;

    /** basically by solving the nonlinear system of equations 
     * \f$ F(x) := x_0 (p_0 - p_1) + x_1 (p_1 - p_2) + x_2 x_0(n_0 - n_2) + x_2 x_1 (n_1 - n_2)
     * + x_2 n_2 + p_2 - p = 0\f$ using standard Newton iteration.
     */
    bool computeInverseNormalProjection(const StaticVector<ctype,3>& p0, const StaticVector<ctype,3>& p1, const StaticVector<ctype,3>& p2,
                                        const StaticVector<ctype,3>& n0, const StaticVector<ctype,3>& n1, const StaticVector<ctype,3>& n2,
                                        const StaticVector<ctype,3>& target, StaticVector<ctype,3>& x);

      
    bool edgeIntersectsNormalFan(const StaticVector<ctype,3>& p0, const StaticVector<ctype,3>& p1,
                                 const StaticVector<ctype,3>& q0, const StaticVector<ctype,3>& q1,
                                 const StaticVector<ctype,3>& n0, const StaticVector<ctype,3>& n1,
                                 StaticVector<ctype,3>& x);
    
    /** The case of parallel ray and triangle is not considered an intersection
     * no matter whether it is or not.
     */
    bool rayIntersectsTriangle(const StaticVector<ctype,3>& basePoint, 
                               const StaticVector<ctype,3>& direction,
                               const StaticVector<ctype,3>& a, const StaticVector<ctype,3>& b, const StaticVector<ctype,3>& c,
                               StaticVector<ctype,2>& localCoords,
                               ctype& normalDist,
                               ctype eps);

    // ///////////////////////////////////////////////////////////////
    //   A few static methods for the 1d-in-2d case.
    // ///////////////////////////////////////////////////////////////
public:
    static bool computeInverseNormalProjection(const StaticVector<ctype,2>& p0,
                                               const StaticVector<ctype,2>& p1,
                                               const StaticVector<ctype,2>& n0,
                                               const StaticVector<ctype,2>& n1,
                                               const StaticVector<ctype,2>& q,
                                               ctype& local);
    
    static bool normalProjection(const StaticVector<ctype,2>& base,
                                 const StaticVector<ctype,2>& direction,
                                 int& bestSegment,
                                 ctype& rangeLocalPosition,
                                 const std::vector<std::tr1::array<int,2> >& targetSegments,
                                 const std::vector<std::tr1::array<ctype, 2> >& coords);
    
    static bool rayIntersectsLine(const StaticVector<ctype, 2>& basePoint, 
                                  const StaticVector<ctype, 2>& direction,
                                  const StaticVector<ctype, 2>& a, 
                                  const StaticVector<ctype, 2>& b, 
                                  ctype& distance, ctype& targetLocal);

protected:
    // /////////////////////////////////////
    // Data members
    // /////////////////////////////////////

    PSurface<2,ctype>* psurface_;

};

} // namespace psurface

#endif

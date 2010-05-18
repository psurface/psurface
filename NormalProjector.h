#ifndef NORMAL_PROJECTOR_H
#define NORMAL_PROJECTOR_H

#include <psurface/StaticVector.h>
#include "PSurface.h"

#include <vector>

class NodeBundle;
class Surface;
class GlobalNodeIdx;
class ContactBoundary;


/** \brief Construct a PSurface object by projecting one surface in normal direction onto another

\tparam ctype The type used for coordinates
*/
template <class ctype>
class NormalProjector {
public:

    NormalProjector(PSurface<2,ctype>* psurface)
        : psurface_(psurface)
    {}

    void project(const ContactBoundary& contactPatch,
                    void (*directions)(const double* pos, double* dir)
                    );

protected:

    void setupEdgePointArrays();

    void insertEdge(const std::vector<StaticVector<double,3> >& normals,
                    int from, 
                    int to,
                    const std::vector<NodeBundle>& projectedTo
                    );

    void insertEdgeFromInteriorNode(const std::vector<StaticVector<double,3> >& normals, 
                                    int from, int to, double& lambda,
                                    const std::vector<NodeBundle>& projectedTo,
                                    NodeBundle& curr, int& enteringEdge);

    void insertEdgeFromIntersectionNode(const std::vector<StaticVector<double,3> >& normals, 
                                        int from, int to, double& lambda,
                                        const std::vector<NodeBundle>& projectedTo,
                                        NodeBundle& curr, int& enteringEdge);

    void insertEdgeFromTouchingNode(const std::vector<StaticVector<double,3> >& normals, 
                                    int from, int to, double& lambda,
                                    const std::vector<NodeBundle>& projectedTo,
                                    NodeBundle& curr, int& enteringTri);

    void insertEdgeFromCornerNode(const std::vector<StaticVector<double,3> >& normals, 
                                  int from, int to, double& lambda,
                                    const std::vector<NodeBundle>& projectedTo,
                                    NodeBundle& curr, int& enteringEdge);


    // ///////////////////////////////////////////////////////////////////////
    //   Methods needed to test whether an edge can be projected completely
    // ///////////////////////////////////////////////////////////////////////

    bool edgeCanBeInserted(const std::vector<StaticVector<double,3> >& normals,
                           int from, 
                           int to,
                           const std::vector<NodeBundle>& projectedTo
                           );

    bool testInsertEdgeFromInteriorNode(const std::vector<StaticVector<double,3> >& normals, 
                                        int from, int to, double& lambda,
                                        NodeBundle& curr,
                                        typename Node<ctype>::NodeType& currType, int& currTri,
                                        int& enteringEdge);

    bool testInsertEdgeFromIntersectionNode(const std::vector<StaticVector<double,3> >& normals, 
                                            int from, int to, double& lambda,
                                            NodeBundle& curr,
                                            typename Node<ctype>::NodeType& currType, int& currTri,
                                            int& enteringEdge);

    bool testInsertEdgeFromTouchingNode(const std::vector<StaticVector<double,3> >& normals, 
                                        int from, int to, double& lambda,
                                        NodeBundle& curr,
                                        typename Node<ctype>::NodeType& currType, int& currTri,
                                        int& enteringEdge);

    bool testInsertEdgeFromCornerNode(const std::vector<StaticVector<double,3> >& normals, 
                                      int from, int to, double& lambda,
                                      NodeBundle& curr, 
                                      typename Node<ctype>::NodeType& currType, int& currTri,
                                      int& enteringEdge);

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
                                        const StaticVector<double,3>& n0, const StaticVector<double,3>& n1, const StaticVector<double,3>& n2,
                                        const StaticVector<ctype,3>& target, StaticVector<double,3>& x);

      
    bool edgeIntersectsNormalFan(const StaticVector<ctype,3>& p0, const StaticVector<ctype,3>& p1,
                                 const StaticVector<ctype,3>& q0, const StaticVector<ctype,3>& q1,
                                 const StaticVector<double,3>& n0, const StaticVector<double,3>& n1,
                                 StaticVector<double,3>& x);
    
    /** The case of parallel ray and triangle is not considered an intersection
     * no matter whether it is or not.
     */
    bool rayIntersectsTriangle(const StaticVector<ctype,3>& basePoint, 
                               const StaticVector<ctype,3>& direction,
                               const StaticVector<ctype,3>& a, const StaticVector<ctype,3>& b, const StaticVector<ctype,3>& c,
                               StaticVector<ctype,2>& localCoords,
                               ctype& normalDist,
                               double eps);

    // ///////////////////////////////////////////////////////////////
    //   A few static methods for the 1d-in-2d case.
    // ///////////////////////////////////////////////////////////////
public:
    static bool computeInverseNormalProjection(const StaticVector<double,2>& p0,
                                               const StaticVector<double,2>& p1,
                                               const StaticVector<double,2>& n0,
                                               const StaticVector<double,2>& n1,
                                               const StaticVector<double,2>& q,
                                               double& local);
    
    static bool normalProjection(const StaticVector<double,2>& base,
                                 const StaticVector<double,2>& direction,
                                 int& bestSegment,
                                 double& rangeLocalPosition,
                                 const std::vector<std::tr1::array<int,2> >& targetSegments,
                                 const std::vector<std::tr1::array<double, 2> >& coords);
    
    static bool rayIntersectsLine(const StaticVector<double, 2>& basePoint, 
                                  const StaticVector<double, 2>& direction,
                                  const StaticVector<double, 2>& a, 
                                  const StaticVector<double, 2>& b, 
                                  double& distance, double& targetLocal);

protected:
    // /////////////////////////////////////
    // Data members
    // /////////////////////////////////////

    PSurface<2,ctype>* psurface_;

    std::vector<StaticVector<double, 3> > targetNormals;

};

#endif

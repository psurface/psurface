#ifndef NORMAL_PROJECTOR_H
#define NORMAL_PROJECTOR_H

#include "McVec3d.h"
#include <psurface/Parametrization.h>

#include <vector>

class NodeBundle;
class McVec3f;
class McVec2f;
class McVec2d;
class Surface;
class ContactBoundary;
class GlobalNodeIdx;


class NormalProjector {
public:

    void handleSide(Parametrization* par, 
                    const ContactBoundary& contactPatch,
                    void (*obsDirections)(const double* pos, double* dir)
                    );

    void setupEdgePointArrays(Parametrization* par);

    void insertEdge(Parametrization* par,  
                    const std::vector<McVec3d>& normals,
                    int from, 
                    int to,
                    const std::vector<NodeBundle>& projectedTo
                    );

    void insertEdgeFromInteriorNode(Parametrization* par, 
                                    const std::vector<McVec3d>& normals, 
                                    int from, int to, double& lambda,
                                    const std::vector<NodeBundle>& projectedTo,
                                    NodeBundle& curr, int& enteringEdge);

    void insertEdgeFromIntersectionNode(Parametrization* par, 
                                        const std::vector<McVec3d>& normals, 
                                        int from, int to, double& lambda,
                                        const std::vector<NodeBundle>& projectedTo,
                                        NodeBundle& curr, int& enteringEdge);

    void insertEdgeFromTouchingNode(Parametrization* par, 
                                    const std::vector<McVec3d>& normals, 
                                    int from, int to, double& lambda,
                                    const std::vector<NodeBundle>& projectedTo,
                                    NodeBundle& curr, int& enteringTri);

    void insertEdgeFromCornerNode(Parametrization* par, 
                                  const std::vector<McVec3d>& normals, 
                                  int from, int to, double& lambda,
                                    const std::vector<NodeBundle>& projectedTo,
                                    NodeBundle& curr, int& enteringEdge);


    // ///////////////////////////////////////////////////////////////////////
    //   Methods needed to test whether an edge can be projected completely
    // ///////////////////////////////////////////////////////////////////////
    bool edgeCanBeInserted(const Parametrization* par,  
                           const std::vector<McVec3d>& normals,
                           int from, 
                           int to,
                           const std::vector<NodeBundle>& projectedTo
                           );

    bool testInsertEdgeFromInteriorNode(const Parametrization* par, 
                                        const std::vector<McVec3d>& normals, 
                                        int from, int to, double& lambda,
                                        const std::vector<NodeBundle>& projectedTo,
                                        Node::NodeType& currType, TriangleIdx& currTri,
                                        int& enteringEdge);

    bool testInsertEdgeFromIntersectionNode(const Parametrization* par, 
                                            const std::vector<McVec3d>& normals, 
                                            int from, int to, double& lambda,
                                            const std::vector<NodeBundle>& projectedTo,
                                            Node::NodeType& currType, TriangleIdx& currTri,
                                            int& enteringEdge);

    bool testInsertEdgeFromTouchingNode(const Parametrization* par, 
                                        const std::vector<McVec3d>& normals, 
                                        int from, int to, double& lambda,
                                        const std::vector<NodeBundle>& projectedTo,
                                        const NodeBundle& curr,
                                        Node::NodeType& currType, TriangleIdx& currTri,
                                        int& enteringEdge);

    bool testInsertEdgeFromCornerNode(const Parametrization* par, 
                                      const std::vector<McVec3d>& normals, 
                                      int from, int to, double& lambda,
                                      const std::vector<NodeBundle>& projectedTo,
                                      const NodeBundle& curr, 
                                      Node::NodeType& currType, TriangleIdx& currTri,
                                      int& enteringEdge);

    void insertGhostNodeAtVertex(Parametrization* par, 
                                 int v, 
                                 int targetTri, 
                                 const McVec2d& localTargetCoords
                                 );

    void addCornerNodeBundle(Parametrization* par, 
                             int v, 
                             int nN
                             );

    /** \brief Return true if the two NodeBundles share at least one triangle */
    bool onSameTriangle(const NodeBundle& a, const NodeBundle& others) const;

    /** \brief Return true if the NodeBundle contains a Node on the given triangle */
    bool onSameTriangle(const TriangleIdx& tri, const NodeBundle& others) const;

    TriangleIdx getCommonTri(const NodeBundle& a, const NodeBundle& b);

    McSmallArray<TriangleIdx, 2> getCommonTris(const NodeBundle& a, const NodeBundle& b);

    NodeIdx getCornerNode(const DomainTriangle& cT, int corner);

    /** basically by solving the nonlinear system of equations 
     * \f$ F(x) := x_0 (p_0 - p_1) + x_1 (p_1 - p_2) + x_2 x_0(n_0 - n_2) + x_2 x_1 (n_1 - n_2)
     * + x_2 n_2 + p_2 - p = 0\f$ using standard Newton iteration.
     */
    bool computeInverseNormalProjection(const McVec3f& p0, const McVec3f& p1, const McVec3f& p2,
                                        const McVec3d& n0, const McVec3d& n1, const McVec3d& n2,
                                        const McVec3f& target, McVec3d& x);

      
    bool edgeIntersectsNormalFan(const McVec3f& p0, const McVec3f& p1,
                                 const McVec3f& q0, const McVec3f& q1,
                                 const McVec3d& n0, const McVec3d& n1,
                                 McVec3d& x);
    
    /** The case of parallel ray and triangle is not considered an intersection
     * no matter whether it is or not.
     */
    bool rayIntersectsTriangle(const McVec3d& basePoint, 
                               const McVec3d& direction,
                               const McVec3f& a, const McVec3f& b, const McVec3f& c,
                               McVec2d& localCoords,
                               double& normalDist,
                               double eps);

    /** \todo Können wir nicht die entsprechende Routine in Surface verwenden? */
    void computeVertexNormals();

    // /////////////////////////////////////
    // Data members

    const ContactBoundary* contactBoundary[2];
    std::vector<McVec3d> targetNormals;

};

#endif

#ifndef NORMAL_PROJECTOR_H
#define NORMAL_PROJECTOR_H

#include <psurface/StaticVector.h>
#include "PSurface.h"

#include <vector>

class NodeBundle;
class Surface;
class ContactBoundary;
class GlobalNodeIdx;


class NormalProjector {
public:

    void handleSide(PSurface<2,float>* par, 
                    const ContactBoundary& contactPatch,
                    void (*obsDirections)(const double* pos, double* dir)
                    );

    void setupEdgePointArrays(PSurface<2,float>* par);

    void insertEdge(PSurface<2,float>* par,  
                    const std::vector<StaticVector<double,3> >& normals,
                    int from, 
                    int to,
                    const std::vector<NodeBundle>& projectedTo
                    );

    void insertEdgeFromInteriorNode(PSurface<2,float>* par, 
                                    const std::vector<StaticVector<double,3> >& normals, 
                                    int from, int to, double& lambda,
                                    const std::vector<NodeBundle>& projectedTo,
                                    NodeBundle& curr, int& enteringEdge);

    void insertEdgeFromIntersectionNode(PSurface<2,float>* par, 
                                        const std::vector<StaticVector<double,3> >& normals, 
                                        int from, int to, double& lambda,
                                        const std::vector<NodeBundle>& projectedTo,
                                        NodeBundle& curr, int& enteringEdge);

    void insertEdgeFromTouchingNode(PSurface<2,float>* par, 
                                    const std::vector<StaticVector<double,3> >& normals, 
                                    int from, int to, double& lambda,
                                    const std::vector<NodeBundle>& projectedTo,
                                    NodeBundle& curr, int& enteringTri);

    void insertEdgeFromCornerNode(PSurface<2,float>* par, 
                                  const std::vector<StaticVector<double,3> >& normals, 
                                  int from, int to, double& lambda,
                                    const std::vector<NodeBundle>& projectedTo,
                                    NodeBundle& curr, int& enteringEdge);


    // ///////////////////////////////////////////////////////////////////////
    //   Methods needed to test whether an edge can be projected completely
    // ///////////////////////////////////////////////////////////////////////
    bool edgeCanBeInserted(const PSurface<2,float>* par,  
                           const std::vector<StaticVector<double,3> >& normals,
                           int from, 
                           int to,
                           const std::vector<NodeBundle>& projectedTo
                           );

    bool testInsertEdgeFromInteriorNode(const PSurface<2,float>* par, 
                                        const std::vector<StaticVector<double,3> >& normals, 
                                        int from, int to, double& lambda,
                                        const std::vector<NodeBundle>& projectedTo,
                                        Node<float>::NodeType& currType, int& currTri,
                                        int& enteringEdge);

    bool testInsertEdgeFromIntersectionNode(const PSurface<2,float>* par, 
                                            const std::vector<StaticVector<double,3> >& normals, 
                                            int from, int to, double& lambda,
                                            const std::vector<NodeBundle>& projectedTo,
                                            Node<float>::NodeType& currType, int& currTri,
                                            int& enteringEdge);

    bool testInsertEdgeFromTouchingNode(const PSurface<2,float>* par, 
                                        const std::vector<StaticVector<double,3> >& normals, 
                                        int from, int to, double& lambda,
                                        const std::vector<NodeBundle>& projectedTo,
                                        const NodeBundle& curr,
                                        Node<float>::NodeType& currType, int& currTri,
                                        int& enteringEdge);

    bool testInsertEdgeFromCornerNode(const PSurface<2,float>* par, 
                                      const std::vector<StaticVector<double,3> >& normals, 
                                      int from, int to, double& lambda,
                                      const std::vector<NodeBundle>& projectedTo,
                                      const NodeBundle& curr, 
                                      Node<float>::NodeType& currType, int& currTri,
                                      int& enteringEdge);

    void insertGhostNodeAtVertex(PSurface<2,float>* par, 
                                 int v, 
                                 int targetTri, 
                                 const StaticVector<double,2>& localTargetCoords
                                 );

    void addCornerNodeBundle(PSurface<2,float>* par, 
                             int v, 
                             int nN
                             );

    /** \brief Return true if the two NodeBundles share at least one triangle */
    bool onSameTriangle(const NodeBundle& a, const NodeBundle& others) const;

    /** \brief Return true if the NodeBundle contains a Node on the given triangle */
    bool onSameTriangle(const int& tri, const NodeBundle& others) const;

    int getCommonTri(const NodeBundle& a, const NodeBundle& b);

    std::vector<int> getCommonTris(const NodeBundle& a, const NodeBundle& b);

    NodeIdx getCornerNode(const DomainTriangle& cT, int corner);

    /** basically by solving the nonlinear system of equations 
     * \f$ F(x) := x_0 (p_0 - p_1) + x_1 (p_1 - p_2) + x_2 x_0(n_0 - n_2) + x_2 x_1 (n_1 - n_2)
     * + x_2 n_2 + p_2 - p = 0\f$ using standard Newton iteration.
     */
    bool computeInverseNormalProjection(const StaticVector<float,3>& p0, const StaticVector<float,3>& p1, const StaticVector<float,3>& p2,
                                        const StaticVector<double,3>& n0, const StaticVector<double,3>& n1, const StaticVector<double,3>& n2,
                                        const StaticVector<float,3>& target, StaticVector<double,3>& x);

      
    bool edgeIntersectsNormalFan(const StaticVector<float,3>& p0, const StaticVector<float,3>& p1,
                                 const StaticVector<float,3>& q0, const StaticVector<float,3>& q1,
                                 const StaticVector<double,3>& n0, const StaticVector<double,3>& n1,
                                 StaticVector<double,3>& x);
    
    /** The case of parallel ray and triangle is not considered an intersection
     * no matter whether it is or not.
     */
    bool rayIntersectsTriangle(const StaticVector<double,3>& basePoint, 
                               const StaticVector<double,3>& direction,
                               const StaticVector<float,3>& a, const StaticVector<float,3>& b, const StaticVector<float,3>& c,
                               StaticVector<double,2>& localCoords,
                               double& normalDist,
                               double eps);

    /** \todo Können wir nicht die entsprechende Routine in Surface verwenden? */
    void computeVertexNormals();

    // /////////////////////////////////////
    // Data members

    const ContactBoundary* contactBoundary[2];
    std::vector<StaticVector<double, 3> > targetNormals;

};

#endif

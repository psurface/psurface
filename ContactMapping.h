#ifndef CONTACT_MAPPING_HH
#define CONTACT_MAPPING_HH

#include <vector>
#include <iostream>
#include <tr1/array>

#include <psurface/StaticVector.h>
#include <psurface/contact.h>

template <int dim>
class ContactMapping {};

template <>
class ContactMapping<2>
{
public:
    void build(const std::vector<std::tr1::array<double,2> >& coords1,  ///< The vertices of the first surface as \f$x_0 ,y_0 ,z_0, x_1, y_1, z_1 ...\f$
               const std::vector<std::tr1::array<int,2> >& tri1,       ///< The triangles of the first surface
               const std::vector<std::tr1::array<double,2> >& coords2,  ///< The vertices of the second surface
               const std::vector<std::tr1::array<int,2> >& tri2,
               float epsilon,   ///< The estimate maximum deformation for the contact oracle
               void (*obsDirections)(const double* pos, double* dir)
               );

    void getOverlaps(std::vector<IntPrimitive>& overlaps);

private:

    bool isCompletelyCovered(int i) const;

    static bool computeInverseNormalProjection(const StaticVector<double,2>& p0,
                                               const StaticVector<double,2>& p1,
                                               const StaticVector<double,2>& n0,
                                               const StaticVector<double,2>& n1,
                                               const StaticVector<double,2>& q,
                                               double& local);

    bool normalProjection(const StaticVector<double,2>& base,
                          const StaticVector<double,2>& direction,
                          int& bestSegment,
                          double& rangeLocalPosition,
                          const std::vector<std::tr1::array<int,2> >& targetSegments,
                          const std::vector<std::tr1::array<double, 2> >& coords) const;

    bool rayIntersectsLine(const StaticVector<double, 2>& basePoint, 
                           const StaticVector<double, 2>& direction,
                           const StaticVector<double, 2>& a, 
                           const StaticVector<double, 2>& b, 
                           double& distance, double& targetLocal) const;

    // /////////////////////////////////////////////

    class Node {
    public:

        Node() {}

        Node(double dLP, double rLP, bool nOV, bool nOTV, int rangeSegment0, int rangeSegment1)
            : domainLocalPosition(dLP), rangeLocalPosition(rLP),
              isNodeOnVertex(nOV), isNodeOnTargetVertex(nOTV),
              rightRangeSegment(-1)
        {
            rangeSegments[0] = rangeSegment0;
            rangeSegments[1] = rangeSegment1;
        }

        friend
        std::ostream& operator<< (std::ostream& s, const Node& node)
        {
            s << node.domainLocalPosition << ",   " << node.rangeLocalPosition << ",   "
              << ((node.isNodeOnVertex) ? "true" : "false") << "  "
              << ((node.isNodeOnTargetVertex) ? "true" : "false") << "  --  "
              << "rangeSegments: " << node.rangeSegments[0] << "  " << node.rangeSegments[1] 
              << " -- " << node.rightRangeSegment << std::endl;
            return s;
        }

        double domainLocalPosition;

        double rangeLocalPosition;

        bool isNodeOnVertex;
        
        bool isNodeOnTargetVertex;

        int rangeSegments[2];

        int rightRangeSegment;
    };

    class DomainSegment {
    public:
        std::vector<Node> nodes;

        int points[2];

        int neighbor[2];
    };

    std::vector<StaticVector<double, 2> > vertices;

    std::vector<DomainSegment> domainSegments;

    std::vector<StaticVector<double, 2> > domainNormals;

    std::vector<StaticVector<double, 2> > targetVertices;

    std::vector<StaticVector<double, 2> > targetNormals;

};


template <>
class ContactMapping<3>
{
public:	

    ~ContactMapping() {
        // delete the contact surface after use
        deleteContactSurface();
    }

    void build(const std::vector<std::tr1::array<double,3> >& coords1,  ///< The vertices of the first surface as \f$x_0 ,y_0 ,z_0, x_1, y_1, z_1 ...\f$
               const std::vector<std::tr1::array<int,3> >& tri1,       ///< The triangles of the first surface
               const std::vector<std::tr1::array<double,3> >& coords2,  ///< The vertices of the second surface
               const std::vector<std::tr1::array<int,3> >& tri2,       ///< The triangles of the second surface
               float epsilon,   ///< The estimate maximum deformation for the contact oracle
               void (*obsDirections)(const double* pos, double* dir)
               ) {
        buildContactMapping(coords1, tri1, coords2, tri2,
                            epsilon,
                            obsDirections);
    }

    void getOverlaps(std::vector<IntPrimitive>& overlaps) {
        getMergedGrid(overlaps);
    }

};

#endif

#ifndef CONTACT_MAPPING_HH
#define CONTACT_MAPPING_HH

#include <vector>

#include <dune/common/fvector.hh>

#include <contact/contact.h>

template <int dim>
class ContactMapping {};

template <>
class ContactMapping<2>
{
public:
    void build(const std::vector<double>& coords1,  ///< The vertices of the first surface as \f$x_0 ,y_0 ,z_0, x_1, y_1, z_1 ...\f$
               const std::vector<int>& tri1,       ///< The triangles of the first surface
               const std::vector<double>& coords2,  ///< The vertices of the second surface
               const std::vector<int>& tri2,
               float epsilon,   ///< The estimate maximum deformation for the contact oracle
               void (*obsDirections)(const double* pos, double* dir)
               );

    void getOverlaps(std::vector<IntPrimitive>& overlaps);

private:

    bool isCompletelyCovered(int i) const;

    static bool computeInverseNormalProjection(const Dune::FieldVector<double,2>& p0,
                                               const Dune::FieldVector<double,2>& p1,
                                               const Dune::FieldVector<double,2>& n0,
                                               const Dune::FieldVector<double,2>& n1,
                                               const Dune::FieldVector<double,2>& q,
                                               double& local);

    bool normalProjection(const Dune::FieldVector<double,2>& base,
                          const Dune::FieldVector<double,2>& direction,
                          int& bestSegment,
                          double& rangeLocalPosition,
                          const std::vector<int>& targetSegments,
                          const std::vector<double>& coords) const;

    bool rayIntersectsLine(const Dune::FieldVector<double, 2>& basePoint, 
                           const Dune::FieldVector<double, 2>& direction,
                           const Dune::FieldVector<double, 2>& a, 
                           const Dune::FieldVector<double, 2>& b, 
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

    std::vector<Dune::FieldVector<double, 2> > vertices;

    std::vector<DomainSegment> domainSegments;

    std::vector<Dune::FieldVector<double, 2> > domainNormals;

    std::vector<Dune::FieldVector<double, 2> > targetVertices;

    std::vector<Dune::FieldVector<double, 2> > targetNormals;

};


template <>
class ContactMapping<3>
{
public:	

    ~ContactMapping() {
        // delete the contact surface after use
        deleteContactSurface();
    }

    void build(const std::vector<double>& coords1,  ///< The vertices of the first surface as \f$x_0 ,y_0 ,z_0, x_1, y_1, z_1 ...\f$
               const std::vector<int>& tri1,       ///< The triangles of the first surface
               const std::vector<double>& coords2,  ///< The vertices of the second surface
               const std::vector<int>& tri2,       ///< The triangles of the second surface
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

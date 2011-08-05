#ifndef POINT_INTERSECTION_FUNCTOR_H
#define POINT_INTERSECTION_FUNCTOR_H

// Check for VC9 / VS2008 with installed feature pack.
#if defined(_MSC_VER) && (_MSC_VER>=1500)
    // Dummy-include to define _CPPLIB_VER.
    #include <vector>

    #if defined(_CPPLIB_VER) && _CPPLIB_VER>=505
        #include <array>
    #else
        #error Please install the Visual Studio 2008 SP1 for TR1 support.
    #endif
#else
    #include <tr1/array>
#endif

#include "SurfaceParts.h"

/** \brief Functor class needed to insert McEdge objects into a MultiDimOctree
 */
struct EdgeIntersectionFunctor
{
    class MyMcBox2f 
        : std::tr1::array<float,4>
    {
    public:
        MyMcBox2f(float a, float b, float c, float d)
        {
            (*this)[0] = a;
            (*this)[1] = b;
            (*this)[2] = c;
            (*this)[3] = d;
        }

        /// Returns TRUE if intersection of given point and Box is not empty
        bool intersect(const StaticVector<float,2> &pt) const {
            return (pt[0]>=(*this)[0] && pt[0]<=(*this)[1] 
                    && pt[1]>=(*this)[2] && pt[1]<=(*this)[3]);
        }

    };

    /** \brief Constructor */
    EdgeIntersectionFunctor(const McVertex<float>* vertices)
        : vertices_(vertices)
    {}

    bool operator()(const std::tr1::array<float,3>& lower,
                    const std::tr1::array<float,3>& upper, const McEdge& item) const {

        // The edge intersects the box if one of the endpoints is within the box
        if (Box<float,3>(lower,upper).contains(vertices_[item.from]) 
            || Box<float,3>(lower,upper).contains(vertices_[item.to]))
            return true;

        // Otherwise, the edge intersects the box if it intersects one of the six faces
        return (intersectsXYPatch( MyMcBox2f(lower[0], upper[0], lower[1], upper[1]), lower[2], &item) ||
                intersectsXYPatch( MyMcBox2f(lower[0], upper[0], lower[1], upper[1]), upper[2], &item) ||
                intersectsXZPatch( MyMcBox2f(lower[0], upper[0], lower[2], upper[2]), lower[1], &item) ||
                intersectsXZPatch( MyMcBox2f(lower[0], upper[0], lower[2], upper[2]), upper[1], &item) ||
                intersectsYZPatch( MyMcBox2f(lower[1], upper[1], lower[2], upper[2]), lower[0], &item) ||
                intersectsYZPatch( MyMcBox2f(lower[1], upper[1], lower[2], upper[2]), upper[0], &item));
    }

protected:

    bool intersectsXYPatch(const MyMcBox2f& rect, float z, const McEdge* item) const {

        const StaticVector<float,3>& f = vertices_[item->from];
        const StaticVector<float,3>& t = vertices_[item->to];

        if ( (f[2] < z && t[2] < z) || (f[2] > z && t[2] > z))
            return false;

        float lambda = (z - f[2]) / (t[2] - f[2]);

        StaticVector<float,2> intersection;
        intersection[0] = f[0] + lambda*(t[0] - f[0]);
        intersection[1] = f[1] + lambda*(t[1] - f[1]);

        return rect.intersect(intersection);
    }



   bool intersectsXZPatch(const MyMcBox2f& rect, float y, const McEdge* item) const {

        const StaticVector<float,3>& f = vertices_[item->from];
        const StaticVector<float,3>& t = vertices_[item->to];

        if ( (f[1] < y && t[1] <y) || (f[1] > y && t[1] > y))
            return false;

        float lambda = (y - f[1]) / (t[1] - f[1]);

        StaticVector<float,2> intersection;
        intersection[0] = f[0] + lambda*(t[0] - f[0]);
        intersection[1] = f[2] + lambda*(t[2] - f[2]);

        return rect.intersect(intersection);
    }

    bool intersectsYZPatch(const MyMcBox2f& rect, float x, const McEdge* item) const {

        const StaticVector<float,3>& f = vertices_[item->from];
        const StaticVector<float,3>& t = vertices_[item->to];

        if ( (f[0] < x && t[0] <x) || (f[0] > x && t[0] > x))
            return false;

        float lambda = (x - f[0]) / (t[0] - f[0]);

        StaticVector<float,2> intersection;
        intersection[0] = f[1] + lambda*(t[1] - f[1]);
        intersection[1] = f[2] + lambda*(t[2] - f[2]);

        return rect.intersect(intersection);
    }

    //const std::vector<McVertex>& vertices_;
    const McVertex<float>* vertices_;

};


#endif

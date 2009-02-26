#ifndef CONTACT_MCVEC3F
#define CONTACT_MCVEC3F

#include <mclib/McVec3f.h>

// octree-insertable
class MyMcVec3f : public McVec3f {
public:
    /// Assignment.
    MyMcVec3f& operator=(const McVec3f& v) {
        x=v.x; y=v.y; z=v.z; return *this;
    }
    
    /// for insertion into McOctree
    bool intersect(const McBox3f& box, void* userData) const {
        return (box.getMin().x <= x) && (x <= box.getMax().x)
            && (box.getMin().y <= y) && (y <= box.getMax().y)
            && (box.getMin().z <= z) && (z <= box.getMax().z);
    }
};

struct MyMcVec3fIntersector
{
    bool operator()(const Box<std::tr1::array<float,3>, 3>& box, const MyMcVec3f* item) const {
        return (box.lower()[0] <= item->x) && (item->x <= box.upper()[0])
            && (box.lower()[1] <= item->y) && (item->y <= box.upper()[1])
            && (box.lower()[2] <= item->z) && (item->z <= box.upper()[2]);
    }

    bool operator()(const std::tr1::array<float,3>& lower, 
                    const std::tr1::array<float,3>& upper, const MyMcVec3f& item) const {
        return (lower[0] <= item.x) && (item.x <= upper[0])
            && (lower[1] <= item.y) && (item.y <= upper[1])
            && (lower[2] <= item.z) && (item.z <= upper[2]);
    }

};


#endif

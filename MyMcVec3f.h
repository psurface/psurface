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


#endif

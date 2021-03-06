#ifndef GLOBAL_NODE_IDX_H
#define GLOBAL_NODE_IDX_H

namespace psurface {

/** This class represents a global index for parametrization nodes on the base grid.
 * It consists of the base grid triangle number and the node index on that triangle.
 */
class GlobalNodeIdx {

public:
    //! Default constructor
    GlobalNodeIdx() : tri(-1), idx(-1) {}

    //! Construct from integers
    GlobalNodeIdx(int t, int i) {
        tri = t;
        idx = i;
    }

    //! Set index
    void setValue(int t, int i) {
        tri = t;
        idx = i;
    }

    //! Check if the index is valid
    bool isValid() const { return tri >= 0;}

    //! Print index
    void print() const {
        printf("tri: %d   idx: %d\n", tri, idx);
    }

    //! Base grid triangle index
    int tri;
    //! Local index of the node
    int idx;

};

} // namespace psurface

#endif

#ifndef GLOBAL_NODE_IDX_H
#define GLOBAL_NODE_IDX_H

namespace psurface {

/** This class represents a global index for parametrization nodes on the base grid.
 * It consists of the base grid triangle number and the node index on that triangle.
 */
class GlobalNodeIdx {

public:

    GlobalNodeIdx() : tri(-1), idx(-1) {}

    GlobalNodeIdx(int t, int i) {
        tri = t;
        idx = i;
    }

    void setValue(int t, int i) {
        tri = t;
        idx = i;
    }

    bool isValid() const { return tri >= 0;}

    void print() const {
        printf("tri: %d   idx: %d\n", tri, idx);
    }

    int tri;

    int idx;

};

} // namespace psurface

#endif

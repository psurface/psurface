#ifndef NODE_BUNDLE_H
#define NODE_BUNDLE_H

#include <McSmallArray.h>
#include <psurface/GlobalNodeIdx.h>

/** \brief Set of graph nodes */
class NodeBundle : public McSmallArray<GlobalNodeIdx, 2> {

public:

    /** \brief Get the index of that node of the bundle that is on triangle 'tri'
        \return -1 if none of the nodes is on triangle 'tri'
    */
    NodeIdx triToIdx(TriangleIdx tri) const {
        for (int i=0; i<size(); i++)
            if ((*this)[i].tri == tri)
                return (*this)[i].idx;

        return -1;
    }

    /** \brief Print the content for debugging */
    void print() const {
        for (int i=0; i<size(); i++)
            printf("triangle: %d,   index: %d\n", (*this)[i].tri, (*this)[i].idx);
    }

    /** \brief Two NodeBundles are equal if they consist of the same nodes in the same order */
    int operator==(const NodeBundle& other) const {
        if (size()!=other.size())
            return false;

        for (int i=0; i<size(); i++)
            if ((*this)[i].tri!=other[i].tri || (*this)[i].idx!=other[i].idx )
                return false;

        return true;
    }

    /** \brief Inequality */
    int operator!=(const NodeBundle& other) const {
        return !( (*this)==other);
    }

};

#endif

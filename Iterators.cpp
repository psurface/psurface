#include "config.h"

#include "PlaneParam.h"

using namespace psurface;

template <class ctype>
typename PlaneParam<ctype>::DirectedEdgeIterator& PlaneParam<ctype>::DirectedEdgeIterator::operator++()
{
    if (neighborIdx < (*nodes)[from()].degree()-1)
        neighborIdx++;
    else {
        do{
            fromNode++;
            if (!isValid())
                return *this;
        }while (!(*nodes)[fromNode].degree());

        neighborIdx = 0;
    }

    return *this;
}

template <class ctype>
void PlaneParam<ctype>::DirectedEdgeIterator::invert()
{
    int other = to();
    for (int i=0; i<(*nodes)[other].degree(); i++)
        if ((*nodes)[other].neighbors(i)==fromNode)
            neighborIdx = i;

    fromNode = other;
}


template <class ctype>
typename PlaneParam<ctype>::UndirectedEdgeIterator& PlaneParam<ctype>::UndirectedEdgeIterator::operator++()
{

    do{
        //nextPseudoEdge(edge);
        if (neighborIdx < (*nodes)[from()].degree()-1)
            neighborIdx++;
        else {
            do{
                fromNode++;
                if (!isValid())
                    return *this;
            }while (!(*nodes)[fromNode].degree());

            neighborIdx = 0;
        }

        if (!isValid())
            return *this;
    }while (!isCorrectlyOriented());

            return *this;
}

template <class ctype>
PlaneParam<ctype>::TriangleIterator::TriangleIterator(const DirectedEdgeIterator& firstEdge)
{
    cE = firstEdge;

    // regular handling
    while (cE.isValid() && !isCorrectlyOriented())
        ++cE;
}

template <class ctype>
typename PlaneParam<ctype>::TriangleIterator& PlaneParam<ctype>::TriangleIterator::operator++()
{
    do {
        ++cE;
    } while (cE.isValid() && !isCorrectlyOriented());

    return *this;
}

template <class ctype>
bool PlaneParam<ctype>::TriangleIterator::isCorrectlyOriented() const
{
    if (cE.getONext().to()!=cE.getDPrev().from() ||
             vertices(2) >= vertices(0) || vertices(2) >= vertices(1))
        return false;

    // if the plane graph contains isolated triangles these triangles
    // appear twice each.  The following code filters out the doubles.
    // This can lead to inconsistent orientation
    DirectedEdgeIterator cEInv = cE;
    cEInv.invert();

    return (PlaneParam::orientation((*cE.nodes)[vertices(0)].domainPos(),
                                 (*cE.nodes)[vertices(1)].domainPos(),
                                 (*cE.nodes)[vertices(2)].domainPos()) == 1);

    // This is a topological test and not well thought out
//     int i, j;

//     for (i=0; i<3; i++) {

//         const Node& p = (*cE.nodes)[vertices(i)];
//         const Node& q = (*cE.nodes)[vertices((i+1)%3)];

//         if (p.isINTERIOR_NODE() || q.isINTERIOR_NODE())
//             continue;

//         for (j=0; j<3; j++)
//             if (p.isOnEdge(j) && q.isOnEdge(j))
//                 return q.getDomainEdgePosition(j) > p.getDomainEdgePosition(j);

//     }

//     if (cEInv.getONext().to()!=cE.getDPrev().from() || cEInv.getDPrev().from()!=cE.getONext().to())
//         return true;

//     return vertices(0) < vertices(1);
}

#ifndef __APPLE__
// ////////////////////////////////////////////////////////
//   Explicit template instantiations.
//   If you need more, you can add them here.
// ////////////////////////////////////////////////////////
namespace psurface {
  template class PSURFACE_EXPORT PlaneParam<float>::DirectedEdgeIterator;
  template class PSURFACE_EXPORT PlaneParam<double>::DirectedEdgeIterator;

  template class PSURFACE_EXPORT PlaneParam<float>::UndirectedEdgeIterator;
  template class PSURFACE_EXPORT PlaneParam<double>::UndirectedEdgeIterator;

  template class PSURFACE_EXPORT PlaneParam<float>::TriangleIterator;
  template class PSURFACE_EXPORT PlaneParam<double>::TriangleIterator;
}

#endif

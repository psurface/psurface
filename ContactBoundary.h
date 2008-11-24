#ifndef CONTACT_BOUNDARY_H
#define CONTACT_BOUNDARY_H

#include <mclib/McBox3f.h>
#include <hxsurface/Surface.h>

class ContactBoundary {
public:

    ContactBoundary() {}

    ContactBoundary(const Surface* surface) {
        surf = surface;
    }

    void init(const Surface* surface) {
        surf = surface;
    }

    McBox3f getBoundingBox() const {
        if (vertices.size()==0) return McBox3f(1, 1, 1, -1, -1, -1);
        
        McBox3f box(surf->points[vertices[0]], surf->points[vertices[0]]);
        
        for (int i=1; i<vertices.size(); i++)
            box.extendBy(surf->points[vertices[i]]);

        return box; 
    }        

    const Surface::Triangle& triangles(int n) const {
        return surf->triangles[triIdx[n]];
    }

    McDArray<int> getVertexOffsets() const {
        McDArray<int> result(surf->points.size());
        result.fill(-1);
        for (int i=0; i<vertices.size(); i++)
            result[vertices[i]] = i;

        return result;
    }

    /// \todo Should be faster
    int containsEdge(int from, int to) const {
        int counter=0;
        for (int i=0; i<triIdx.size(); i++)
            for (int j=0; j<3; j++)
                if ((triangles(i).points[j]==from && triangles(i).points[(j+1)%3]==to) ||
                    (triangles(i).points[j]==to && triangles(i).points[(j+1)%3]==from))
                    counter++;

        return counter;
    }

    McDArray<int> vertices;

    McDArray<int> triIdx;

    const Surface* surf;
};

#endif

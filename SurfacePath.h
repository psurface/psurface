#ifndef SURFACE_PATH
#define SURFACE_PATH

    class SurfacePath {

    public:
        /** Indices to Vertex coordinates.  For closed lines the first
            index should be replicated as last one. */
        McDArray<VertexIdx> points;
        McDArray<char> isFix;
        int selected;

        ///
        bool isLoop() const {
            return points.size()>1 && points[0]==points.last();
        }

        ///
        void print() const {
            for (int i=0; i<points.size(); i++)
                printf("PointIdx: %d  isFix: %d\n", points[i], isFix[i]);
        }

        ///
        void removePoint(VertexIdx p) {

            // special handling for loops
            if (points.size()>1 && points[0]==p && points.last()==p) {
                points[0] = points[points.size()-2];
                isFix[0]  = isFix[isFix.size()-2];
                points.removeLast();
                isFix.removeLast();
            }
            
            for (int i=points.size()-1; i>=0; i--)
                if (points[i] == p) {
                    points.remove(i);
                    isFix.remove(i);
                }
            
        }

    };

#endif

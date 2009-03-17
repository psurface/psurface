#ifndef SURFACE_PATH
#define SURFACE_PATH

    class SurfacePath {

    public:
        /** Indices to Vertex coordinates.  For closed lines the first
            index should be replicated as last one. */
        std::vector<int> points;
        std::vector<char> isFix;
        int selected;

        ///
        bool isLoop() const {
            return points.size()>1 && points[0]==points.back();
        }

        ///
        void print() const {
            for (size_t i=0; i<points.size(); i++)
                printf("PointIdx: %d  isFix: %d\n", points[i], isFix[i]);
        }

        ///
        void removePoint(int p) {

            // special handling for loops
            if (points.size()>1 && points[0]==p && points.back()==p) {
                points[0] = points[points.size()-2];
                isFix[0]  = isFix[isFix.size()-2];
                points.pop_back();
                isFix.pop_back();
            }
            
            for (int i=points.size()-1; i>=0; i--)
                if (points[i] == p) {
                    points.erase(points.begin() + i);
                    isFix.erase(isFix.begin() + i);
                }
            
        }

    };

#endif

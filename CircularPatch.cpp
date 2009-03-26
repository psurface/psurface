#include <limits>

#include <psurface/CircularPatch.h>
#include <psurface/Parametrization.h>


bool CircularPatch::inducesTopologyChange() const
{
    int i;

    for (i=0; i<size()-1; i++){

        if (par->findEdge(par->triangles(triangles[i]).vertices[0], 
                          par->triangles(triangles[i]).vertices[2]) != -1){
            //printf("POSSIBLE TOPOLOGY CHANGE FOUND!\n");
            return true;
        }
    }

    return false;
}

bool CircularPatch::hasSmallDihedralAngles(float threshold, const Parametrization* par, 
                                           const DomainVertex* centerVertex) const
{
    printf("hasSmallDihedralAngles has been commented out!\n");
#if 0
    int i, j, k;

    for (i=0; i<triangles.size(); i++) {
        
        const DomainTriangle* cT = triangles[i];
        
        for (j=0; j<3; j++) {
            
            // The triangles are not edge-connected..., i.e. cT->edges is empty
            const DomainEdge* cE = par->findEdge(cT->points[j], cT->points[(j+1)%3]);

            if (cE){ 
                //////////////////////////////////////////////////////
                // this edge is a boundary between the CircularPatch and #par#
         
                for (k=0; k<cE->triangles.size(); k++) {
                    
                    if (cE->triangles[k] == cT || cE->triangles[k]->isConnectedTo(centerVertex))
                        continue;
                    
                    if (cT->dihedralAngle(cE->triangles[k]) < threshold)
                        return true;
                    
                }
            } else {
                ///////////////////////////////////////////////////////
                // this is an edge between two triangles of the CircularPatch
                for (k=i+1; k<triangles.size(); k++) {
                    
                    if (triangles[k]->isConnectedTo(cT->points[j]) &&
                        triangles[k]->isConnectedTo(cT->points[(j+1)%3])) {

                        if (cT->dihedralAngle(triangles[k]) < threshold)
                        return true;
                        
                    }
                }
            }
            
        }
        
    }

#endif
    return false;
}


//////////////////////////////////////////////////////////////////
// this routine returns the bounding box of the patch
// it is not well programmed.  Each vertex is checked three times
void CircularPatch::getBoundingBox(Box<float,3> &bbox) const
{
    assert(size());

    bbox.set(par->vertices(par->triangles(triangles[0]).vertices[0]), 
             par->vertices(par->triangles(triangles[0]).vertices[1]));
    bbox.extendBy( par->vertices(par->triangles(triangles[0]).vertices[2]));

    for (int i=1; i<size(); i++)
        for (int k=0; k<3; k++)
            bbox.extendBy(par->vertices(par->triangles(triangles[i]).vertices[k]));

}



//////////////////////////////////////////////////////////////////
// gives the distance of a point to the patch
float CircularPatch::distanceTo(const StaticVector<float,3> &p) const
{
    int i, j;
    float bestDist = std::numeric_limits<float>::max();
    
    // check point against triangles
    for (j=0; j<size(); j++){

        const DomainTriangle& cT = par->triangles(triangles[j]);
        
        StaticVector<float,3> triPoints[3];
        triPoints[0] = par->vertices(cT.vertices[0]);
        triPoints[1] = par->vertices(cT.vertices[1]);
        triPoints[2] = par->vertices(cT.vertices[2]);
        
        // local base
        StaticVector<float,3> a = triPoints[1] - triPoints[0];
        StaticVector<float,3> b = triPoints[2] - triPoints[0];
        StaticVector<float,3> c = a.cross(b);
        c.normalize();
        
        StaticVector<float,3> x = p - triPoints[0];
        
        // write x in the new base  (Cramer's rule)
        StaticMatrix<float,3> numerator(a, b, c);
        StaticMatrix<float,3> alphaMat(x, b, c);
        StaticMatrix<float,3> betaMat(a, x, c);
        StaticMatrix<float,3> gammaMat(a, b, x);
        
        float alpha = alphaMat.det()/numerator.det();
        float beta  = betaMat.det()/numerator.det();
        float gamma = gammaMat.det()/numerator.det();
        
        // check whether orthogonal projection onto the ab plane is in triangle
        bool isIn = alpha>=0 && beta>=0 && (1-alpha-beta)>=0;
        
        if (isIn && fabs(gamma)<bestDist){
            
            //              printf("a(%1.2f %1.2f %1.2f) b(%1.2f %1.2f %1.2f)  c(%1.2f %1.2f %1.2f)  x(%1.2f %1.2f %1.2f)\n",
            //                     a.x, a.y, a.z, b.x, b.y, b.z, c.x, c.y, c.z, x.x, x.y, x.z);
            
            //              printf("tri: %d, alpha = %f, beta = %f, gamma = %f\n", j, alpha, beta, gamma);
            
            bestDist = fabs(gamma);
        }
    }    
    
    // check point against edges
    for (i=0; i<size(); i++){
        for (j=0; j<3; j++){

            const DomainTriangle& cT = par->triangles(triangles[i]);

            StaticVector<float,3> from = par->vertices(cT.vertices[j]);
            StaticVector<float,3> to   = par->vertices(cT.vertices[(j+1)%3]);
            
            StaticVector<float,3> edge = to - from;
            
            float projectLength = edge.dot(p - from)/edge.length();
            StaticVector<float,3> projection = edge/edge.length() * projectLength;
            
            float orthoDist = ((p-from) - projection).length();
            
            if (projectLength>=0 && projectLength<=edge.length() && orthoDist<bestDist)
                bestDist = orthoDist;       
        }
    }
    
    // check point against vertices
    for (i=0; i<size(); i++){
        for (j=0; j<3; j++){
            float dist = (p - par->vertices(par->triangles(triangles[i]).vertices[j])).length();
            if (dist < bestDist){
                bestDist = dist;
            }
        }
    }
    
    return bestDist;
}

bool CircularPatch::intersectsParametrization(const std::vector<int> &closeEdges) const
{
    for (size_t i=0; i<closeEdges.size(); i++){
        
//         printf("closeEdge:  %d    %d --> %d\n", closeEdges[i],
//                par->edges(closeEdges[i]).from, par->edges(closeEdges[i]).to);

        int from = par->edges(closeEdges[i]).from;
        int to   = par->edges(closeEdges[i]).to;
        
        for (int j=0; j<size(); j++){
            
            // printf("tri: %d  von %d (%d %d %d)\n", j, size(),
//                    par->triangles(triangles[j]).vertices[0],
//                    par->triangles(triangles[j]).vertices[1],
//                    par->triangles(triangles[j]).vertices[2]);

            // check whether triangle and edge have one common point
            if (par->triangles(triangles[j]).isConnectedTo(from) || 
                par->triangles(triangles[j]).isConnectedTo(to) ) 
                continue;
            //printf("a\n");
            
            //if (triangles[j]->intersects(closeEdges[i], 0.00001)){
            if (par->intersectionTriangleEdge(triangles[j], &par->edges(closeEdges[i]), 0.00001)){
                return true;
            }
            //printf("a\n");

        }
    }
    
    return false; 
}


bool CircularPatch::hasSelfintersections() const
{
    DomainEdge tmpEdge;
    
    for (size_t i=0; i<innerEdges.size(); i++){
        
        tmpEdge.from = innerEdges[i][0];
        tmpEdge.to   = innerEdges[i][1];
        
        for (int j=0; j<size(); j++){
            
            // check whether triangle and edge have one common point
            if (par->triangles(triangles[j]).isConnectedTo(tmpEdge.from) 
                || par->triangles(triangles[j]).isConnectedTo(tmpEdge.to) ) 
                continue;
            
            //if (triangles[j]->intersects(&tmpEdge, 0.00001)){
            if (par->intersectionTriangleEdge(triangles[j], &tmpEdge, 0.00001)){
                return true;
            }
        }
    }
    
    return false;
}







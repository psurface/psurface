#ifndef QUALITY_REQUEST
#define QUALITY_REQUEST

namespace psurface {

/// A class encapsulating the exact definition of a 'good' point removal step
struct QualityRequest {

    /// Are surface intersections allowed?
    bool intersections;

    /// An edge length that mustn't be exceeded.  Negative values mean 'no restriction'.
    float maxEdgeLength;

    /// Should small dihedral angles be allowed?
    bool smallDihedralAngles;

    /// Minimum allowed dihedral angle.
    float dihedralAngleThreshold;

    /// Importance of the triangle aspect ratio in the new triangulation
    float aspectRatio;

    /// Importance of the Hausdorff distance in the new triangulation
    float hausdorffDistance;


    // //////////////////////
    // Member functions

    /** Default constructor.  Deactivates all restrictions except for paths and intents
     * equal rights for aspect ratio and Hausdorff distance.
     */
    QualityRequest() {
        intersections = false;
        maxEdgeLength = -1;
        smallDihedralAngles = false;
        dihedralAngleThreshold = 0;
        aspectRatio = 0.5;
        hausdorffDistance = 0.5;
    }

    void normalize() {
        float sum = aspectRatio + hausdorffDistance;
        aspectRatio       /= sum;
        hausdorffDistance /= sum;
    }

};

} // namespace

#endif

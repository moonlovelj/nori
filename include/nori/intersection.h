//
// Created by 郭彬 on 2022/3/18.
//

#ifndef NORI_INTERSECTION_H
#define NORI_INTERSECTION_H

#include <nori/common.h>
#include <nori/frame.h>
#include <nori/bbox.h>
#include <nori/dpdf.h>
#include <nori/phasefunction.h>
#include <nori/medium.h>
#include <nori/mediuminterface.h>

NORI_NAMESPACE_BEGIN

/**
 * \brief Intersection data structure
 *
 * This data structure records local information about a ray-triangle intersection.
 * This includes the position, traveled ray distance, uv coordinates, as well
 * as well as two local coordinate frames (one that corresponds to the true
 * geometry, and one that is used for shading computations).
 */
struct Intersection {
    /// Position of the surface intersection
    Point3f p;
    /// Unoccluded distance along the ray
    float t;
    /// UV coordinates, if any
    Point2f uv;
    /// Shading frame (based on the shading normal)
    Frame shFrame;
    /// Geometric frame (based on the true geometry)
    Frame geoFrame;
    /// Pointer to the associated mesh
    const Mesh *mesh;
    /// medium
    MediumInterface mediumInterface;
//    /// this if only for an interaction at a point in a scattering medium
//    std::shared_ptr<const PhaseFunction> phase;
    bool insideMedium;

    /// Create an uninitialized intersection record
    Intersection() : mesh(nullptr), mediumInterface(MediumInterface(nullptr)) { }

    /// Transform a direction vector into the local shading frame
    Vector3f toLocal(const Vector3f &d) const {
        return shFrame.toLocal(d);
    }

    /// Transform a direction vector from local to world coordinates
    Vector3f toWorld(const Vector3f &d) const {
        return shFrame.toWorld(d);
    }

    const Medium* getMedium(const Vector3f &w) const {
        return w.dot(shFrame.n) > 0 ? mediumInterface.m_outside :
               mediumInterface.m_inside;
    }

    /// For interactions that are known to be inside participating media
    const Medium* getMedium() const {
        assert(mediumInterface.m_inside == mediumInterface.m_outside);
        return mediumInterface.m_inside;
    }

    bool IsMedium() const {
        return insideMedium;
    }

    /// Return a human-readable summary of the intersection record
    std::string toString() const;
};

NORI_NAMESPACE_END

#endif //NORI_INTERSECTION_H

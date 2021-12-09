/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
*/

#pragma once

#include <nori/object.h>
#include <nori/mesh.h>

NORI_NAMESPACE_BEGIN

struct EmitterQueryRecord {
    /// Point in world space
    Point3f point;
    /// Normal in world space
    Vector3f normal;

    EmitterQueryRecord() {}

    EmitterQueryRecord(const Point3f &p, const Vector3f &n) : point(p), normal(n) {}
};

/**
 * \brief Superclass of all emitters
 */
class Emitter : public NoriObject {
public:
    virtual void build() = 0;

    virtual Color3f eval(const EmitterQueryRecord &eRec) const = 0;

    virtual Color3f sample(const Point3f &surfacePoint, EmitterQueryRecord &eRec, const Point2f &sample) const = 0;

    virtual float pdf(const EmitterQueryRecord &bRec) const = 0;

    /**
     * \brief Return the type of object (i.e. Mesh/Emitter/etc.) 
     * provided by this instance
     * */
    EClassType getClassType() const { return EEmitter; }

    virtual Color3f emission() const = 0;
};

NORI_NAMESPACE_END

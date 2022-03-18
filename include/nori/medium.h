//
// Created by 郭彬 on 2022/1/19.
//
#ifndef NORI_MEDIUM_H
#define NORI_MEDIUM_H

#include <nori/common.h>
#include <nori/Ray.h>
#include <nori/color.h>
#include <nori/intersection.h>
#include <nori/phasefunction.h>
#include <nori/object.h>

NORI_NAMESPACE_BEGIN

struct Intersection;

class Medium : public NoriObject{
public:
    EClassType getClassType() const override { return EMedium; }
    virtual Color3f tr(const Ray3f &ray, Sampler *sampler) const = 0;
    virtual Color3f sample(const Ray3f &ray, Sampler *sampler, Intersection &its) const = 0;

protected:
    std::shared_ptr<PhaseFunction> m_phase;
};


NORI_NAMESPACE_END

#endif //NORI_MEDIUM_H

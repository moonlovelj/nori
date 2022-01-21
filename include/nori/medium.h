//
// Created by 郭彬 on 2022/1/19.
//

#ifndef NORI_MEDIUM_H
#define NORI_MEDIUM_H

#include <nori/common.h>
#include <nori/Ray.h>
#include <nori/color.h>
#include <nori/mesh.h>

NORI_NAMESPACE_BEGIN

struct Intersection;

class Medium {
public:
    virtual Color3f tr(const Ray3f &ray, Sampler *sampler) const = 0;
    virtual Color3f sample(const Ray3f &ray, Sampler *sampler, Intersection &its) const = 0;
};

struct MediumInterface {
    MediumInterface(std::shared_ptr<const Medium> medium)
            : m_inside(medium), m_outside(medium) {}

    MediumInterface(std::shared_ptr<const Medium> inside, std::shared_ptr<const Medium> outside)
            : m_inside(inside), m_outside(outside) {}

    bool isMediumTransition() const { return m_inside != m_outside; }

    std::shared_ptr<const Medium> m_inside;
    std::shared_ptr<const Medium> m_outside;
};

NORI_NAMESPACE_END

#endif //NORI_MEDIUM_H

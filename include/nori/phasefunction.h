//
// Created by 郭彬 on 2022/1/18.
//

#ifndef NORI_PHASEFUNCTION_H
#define NORI_PHASEFUNCTION_H

#include <nori/common.h>
#include <nori/vector.h>

NORI_NAMESPACE_BEGIN

class PhaseFunction {
public:
    virtual ~PhaseFunction() {}

    virtual float p(const Vector3f &wi, const Vector3f &wo) const = 0;

    virtual float sample(const Vector3f &wi, Vector3f &wo, const Point2f &sample) const = 0;
};

class HenyeyGreenstein : public PhaseFunction {
public:
    HenyeyGreenstein(float g) : m_g(g) {}

    ~HenyeyGreenstein() {}

    float p(const Vector3f &wi, const Vector3f &wo) const override;

    float sample(const Vector3f &wi, Vector3f &wo, const Point2f &sample) const override;

private:
    float m_g;
};

NORI_NAMESPACE_END
#endif //NORI_PHASEFUNCTION_H

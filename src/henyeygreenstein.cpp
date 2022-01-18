//
// Created by 郭彬 on 2022/1/18.
//
#include <nori/phasefunction.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class HenyeyGreenstein : PhaseFunction {
public:
    HenyeyGreenstein(float g) : m_g(g) {}

    float p(const Vector3f &wi, const Vector3f &wo) const override {
        return Warp::phaseHG(wi.dot(wo), m_g);
    }

private:
    float m_g;
};


NORI_NAMESPACE_END
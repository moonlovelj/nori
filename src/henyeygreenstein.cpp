//
// Created by 郭彬 on 2022/1/18.
//
#include <nori/phasefunction.h>
#include <nori/warp.h>
#include <nori/vector.h>

NORI_NAMESPACE_BEGIN

float HenyeyGreenstein::p(const Vector3f &wi, const Vector3f &wo) const {
    return Warp::phaseHG(wi.dot(wo), m_g);
}

float HenyeyGreenstein::sample(const Vector3f &wi, Vector3f &wo, const Point2f &sample) const {
    float cosTheta;
    if (std::abs(m_g) < 1e-3)
        cosTheta = 1 - 2 * sample.x();
    else {
        float sqrTerm = (1 - m_g * m_g) /
                        (1 - m_g + 2 * m_g * sample.x());
        cosTheta = (1 + m_g * m_g - sqrTerm * sqrTerm) / (2 * m_g);
    }

    float phi = 2 * M_PI * sample.y();
    Vector3f v1, v2;
    coordinateSystem(wi, v1, v2);
    Vector3f sphericalCoord = sphericalDirection(std::acos(cosTheta), phi);
    wo = v1 * sphericalCoord.x() + v2 * sphericalCoord.y() + wi * sphericalCoord.z();

    return Warp::phaseHG(cosTheta, m_g);
}

NORI_NAMESPACE_END
//
// Created by 郭彬 on 2022/1/18.
//
#include <nori/phasefunction.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

float HenyeyGreenstein::p(const Vector3f &wi, const Vector3f &wo) const {
return Warp::phaseHG(wi.dot(wo), m_g);
}

float HenyeyGreenstein::sample(const Vector3f &wi, Vector3f &wo, const Point2f &sample) const {
return 1.f;
}

NORI_NAMESPACE_END
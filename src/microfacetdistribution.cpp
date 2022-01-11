//
// Created by 郭彬 on 2022/1/11.
//
#include <nori/microfacetdistribution.h>
#include <nori/frame.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

float BeckmannDistribution::D(const Vector3f &wh, const Vector3f &wi) const {
    float pdf = Warp::squareToBeckmannPdf(wh, scaleAlpha(wi));
    if (pdf == 0) return 0;
    return pdf / Frame::cosTheta(wh);
}

float BeckmannDistribution::G1(const Vector3f &w, const Vector3f &wh, const Vector3f &wi) const {
    float cosThetaVN = Frame::cosTheta(w);
    if (cosThetaVN == 0.f) return 0;
    float sinThetaVN = std::sqrtf(1 - cosThetaVN * cosThetaVN);
    float tanThetaVN = std::fabsf(sinThetaVN / cosThetaVN);
    if (tanThetaVN == 0.f) return 0;

    float b = 1.f / (scaleAlpha(wi) * tanThetaVN);
    if (w.dot(wh) * cosThetaVN < 0) return 0;
    return (b < 1.6f ? (3.535f * b + 2.181f * b * b) / (1 + 2.276f * b + 2.577 * b * b) : 1);
}

Vector3f BeckmannDistribution::sample_wh(const Vector3f &wi, const Point2f &sample) const {
    return Warp::squareToBeckmann(sample, scaleAlpha(wi));
}

std::string BeckmannDistribution::toString() const {
    return "";
}

NORI_NAMESPACE_END

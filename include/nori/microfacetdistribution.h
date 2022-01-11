//
// Created by 郭彬 on 2022/1/11.
//

#ifndef NORI_MICROFACETDISTRIBUTION_H
#define NORI_MICROFACETDISTRIBUTION_H

#include <nori/common.h>
#include <nori/frame.h>

NORI_NAMESPACE_BEGIN
class MicrofacetDistribution {
public:
    virtual ~MicrofacetDistribution() {}

    virtual float D(const Vector3f &wh, const Vector3f &wi) const = 0;

    virtual float G1(const Vector3f &w, const Vector3f &wh, const Vector3f &wi) const = 0;

    virtual Vector3f sample_wh(const Vector3f &wi, const Point2f &sample) const = 0;

    virtual std::string toString() const = 0;

    float G(const Vector3f &wi, const Vector3f &wo, const Vector3f &wh) const {
        return G1(wi, wh, wi) * G1(wo, wh, wi);
    }

    float pdf(const Vector3f &wo, const Vector3f &wh, const Vector3f &wi) const {
        return D(wh, wi) * Frame::cosTheta(wh);
    }
};

class BeckmannDistribution : public MicrofacetDistribution {
public:
    BeckmannDistribution(float alpha) : m_alpha(alpha) {}

    float D(const Vector3f &wh, const Vector3f &wi) const override;

    float G1(const Vector3f &w, const Vector3f &wh, const Vector3f &wi) const override;

    Vector3f sample_wh(const Vector3f &wi, const Point2f &sample) const override;

    std::string toString() const override;

private:
    float scaleAlpha(const Vector3f &wi) const {
        return (1.2f - .2f * sqrtf(Frame::absCosTheta(wi))) * m_alpha;
    }

    float m_alpha;
};

NORI_NAMESPACE_END
#endif //NORI_MICROFACETDISTRIBUTION_H

/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
*/

#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/microfacetdistribution.h>

NORI_NAMESPACE_BEGIN

class MicrofacetTransmission : public BSDF {
public:
    MicrofacetTransmission(const PropertyList &propList) {
        /* RMS surface roughness */
        m_alpha = propList.getFloat("alpha", 0.1f);

        /* Interior IOR (default: BK7 borosilicate optical glass) */
        m_intIOR = propList.getFloat("intIOR", 1.5046f);

        /* Exterior IOR (default: air) */
        m_extIOR = propList.getFloat("extIOR", 1.000277f);

        m_microfacetDistribution = std::make_shared<BeckmannDistribution>(m_alpha);
    }

    /// Evaluate the BRDF for the given pair of directions
    Color3f eval(const BSDFQueryRecord &bRec) const {
        if (bRec.measure != ESolidAngle) {
            return Color3f(0);
        }

        Vector3f wi = bRec.wi;
        Vector3f wo = bRec.wo;
        if (sameHemisphere(wi, wo)) return 0;
        if ((wo.array() == 0).all()) return 0;

        float absCosThetaI = std::fabsf(Frame::cosTheta(wi));
        float absCosThetaO = std::fabsf(Frame::cosTheta(wo));
        if (absCosThetaI == 0.f || absCosThetaO == 0.f) return 0;

        float eta = Frame::cosTheta(wi) > 0 ? m_intIOR / m_extIOR : m_extIOR / m_intIOR;
        Vector3f wh = -(eta * wo + wi).normalized();
        if (wh.z() < 0) wh = -wh;
        if (wo.dot(wh) * wi.dot(wh) > 0) return 0;

//        float absCosThetaH = std::fabsf(Frame::cosTheta(wh));
//        if (absCosThetaH == 0.f) return 0;

        float D = m_microfacetDistribution->D(wh, wi);
        float F = fresnel(wh.dot(wi), m_extIOR, m_intIOR);
        float sqrtDenom = eta * wo.dot(wh) + wi.dot(wh);
        float dwh_dwo_dwi = std::fabsf(wo.dot(wh)) * std::fabsf(wi.dot(wh)) / (sqrtDenom * sqrtDenom * absCosThetaI * absCosThetaO);
        float G = m_microfacetDistribution->G(wi, wo, wh);
        return (1.f - F) * D * G * dwh_dwo_dwi;
    }

    /// Evaluate the sampling density of \ref sample() wrt. solid angles
    float pdf(const BSDFQueryRecord &bRec) const {
        if (bRec.measure != ESolidAngle) {
            return 0.f;
        }

        Vector3f wi = bRec.wi;
        Vector3f wo = bRec.wo;
        if (sameHemisphere(wo, wi)) return 0;

        float eta = Frame::cosTheta(wi) > 0 ? m_intIOR / m_extIOR : m_extIOR / m_intIOR;
        Vector3f wh = -(eta * wo + wi).normalized();
        if (wo.dot(wh) * wi.dot(wh) > 0) return 0;

        float sqrtDenom = eta * wo.dot(wh) + wi.dot(wh);
        float dwh_dwo = eta * eta * std::fabsf(wo.dot(wh)) / (sqrtDenom * sqrtDenom);
        return m_microfacetDistribution->pdf(wo, wh, wi) * dwh_dwo;
    }

    /// Sample the BRDF
    Color3f sample(BSDFQueryRecord &bRec, const Point2f &_sample) const {
        // Note: Once you have implemented the part that computes the scattered
        // direction, the last part of this function should simply return the
        // BRDF value divided by the solid angle density and multiplied by the
        // cosine factor from the reflection equation, i.e.
        // return eval(bRec) * Frame::cosTheta(bRec.wo) / pdf(bRec);

        if (Frame::cosTheta(bRec.wi) == 0.f) return Color3f(0);
        Vector3f m = m_microfacetDistribution->sample_wh(bRec.wi, _sample);
        if ((Frame::cosTheta(bRec.wi) > 0.f ? m : -m).dot(bRec.wi) < 0){
            return Color3f(0);
        }

        //refract
        bRec.wo = refract(m, bRec.wi, m_extIOR, m_intIOR);
        if ((bRec.wo.array() == 0).all()) return Color3f(0);

        bRec.measure = ESolidAngle;
        float p = pdf(bRec);
        if (p == 0.f) return Color3f(0);

        return eval(bRec) * Frame::absCosTheta(bRec.wo) / p;
    }

    bool isDiffuse() const {
        /* While microfacet BRDFs are not perfectly diffuse, they can be
           handled by sampling techniques for diffuse/non-specular materials,
           hence we return true here */
        return true;
    }

    std::string toString() const {
        return tfm::format(
                "MicrofacetTransmission[\n"
                "  alpha = %f,\n"
                "  intIOR = %f,\n"
                "  extIOR = %f,\n"
                "]",
                m_alpha,
                m_intIOR,
                m_extIOR
        );
    }
private:

//    float scaleAlpha(const Vector3f &wi) const {
//        return (1.2f-.2f*sqrtf(Frame::absCosTheta(wi))) * m_alpha;
//    }

    float m_alpha;
    float m_intIOR, m_extIOR;
    std::shared_ptr<MicrofacetDistribution> m_microfacetDistribution;
};

NORI_REGISTER_CLASS(MicrofacetTransmission, "microfacettransmission");
NORI_NAMESPACE_END

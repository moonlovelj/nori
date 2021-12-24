/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
*/

#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class MicrofacetDielectric : public BSDF {
public:
    MicrofacetDielectric(const PropertyList &propList) {
        /* RMS surface roughness */
        m_alpha = propList.getFloat("alpha", 0.1f);

        /* Interior IOR (default: BK7 borosilicate optical glass) */
        m_intIOR = propList.getFloat("intIOR", 1.5046f);

        /* Exterior IOR (default: air) */
        m_extIOR = propList.getFloat("extIOR", 1.000277f);
    }

    /// Evaluate the BRDF for the given pair of directions
    Color3f eval(const BSDFQueryRecord &bRec) const {
        if (bRec.measure != ESolidAngle) {
            return Color3f(0);
        }

        return Color3f(evalR(bRec.wi, bRec.wo) + evalT(bRec.wi, bRec.wo));
    }

    /// Evaluate the sampling density of \ref sample() wrt. solid angles
    float pdf(const BSDFQueryRecord &bRec) const {
        if (bRec.measure != ESolidAngle) {
            return 0.f;
        }
        if (Frame::cosTheta(bRec.wi) == 0.f) return 0;

        return pdfR(bRec.wi, bRec.wo) + pdfT(bRec.wi, bRec.wo);
    }

    /// Sample the BRDF
    Color3f sample(BSDFQueryRecord &bRec, const Point2f &_sample) const {
        // Note: Once you have implemented the part that computes the scattered
        // direction, the last part of this function should simply return the
        // BRDF value divided by the solid angle density and multiplied by the
        // cosine factor from the reflection equation, i.e.
        // return eval(bRec) * Frame::cosTheta(bRec.wo) / pdf(bRec);


        if (Frame::cosTheta(bRec.wi) == 0.f) return Color3f(0);

        Vector3f m = Warp::squareToBeckmann(_sample, scaleAlpha(bRec.wi));
        if ((Frame::cosTheta(bRec.wi) > 0.f ? m : -m).dot(bRec.wi) < 0){
            return Color3f(0);
        }
        if (Warp::squareToBeckmannPdf(m, scaleAlpha(bRec.wi)) == 0.f) return Color3f(0);

        float F = fresnel(m.dot(bRec.wi), m_extIOR, m_intIOR);
        std::srand(std::time(nullptr)); // use current time as seed for random generator
        float randomVariable = std::rand() / (float)RAND_MAX;
        if (randomVariable <= F) {
            // reflect
            bRec.wo = reflect(bRec.wi, m);
            if (!sameHemisphere(bRec.wi, bRec.wo)) return Color3f(0);
        } else {
            //refract
            bRec.wo = refract(m, bRec.wi, m_extIOR, m_intIOR);
            if (sameHemisphere(bRec.wi, bRec.wo)) return Color3f(0);
            if ((bRec.wo.array() == 0).all()) return Color3f(0);
        }

        bRec.measure = ESolidAngle;
        float p = pdf(bRec);
        if (p == 0) return Color3f(0);

        return eval(bRec) * std::fabsf(Frame::cosTheta(bRec.wo)) / p;
    }

    bool isDiffuse() const {
        /* While microfacet BRDFs are not perfectly diffuse, they can be
           handled by sampling techniques for diffuse/non-specular materials,
           hence we return true here */
        return true;
    }

    std::string toString() const {
        return tfm::format(
                "MicrofacetDielectric[\n"
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
    float G1(const Vector3f &wv, const Vector3f &wh, float alpha) const {
        float cosThetaVN = Frame::cosTheta(wv);
        if (cosThetaVN == 0.f) return 0;
        float sinThetaVN = std::sqrtf(1 - cosThetaVN * cosThetaVN);
        float tanThetaVN = std::fabsf(sinThetaVN / cosThetaVN);
        if (tanThetaVN == 0.f) return 0;

        float b = 1.f / (alpha * tanThetaVN);
        if (wv.dot(wh) * cosThetaVN < 0) return 0;
        return (b < 1.6f ? (3.535f*b+2.181f*b*b)/(1+2.276f*b+2.577*b*b):1);
    }

    float G(const Vector3f &wi, const Vector3f &wo, const Vector3f &wh) const {
        float alpha = scaleAlpha(wi);
        return G1(wi, wh, alpha) * G1(wo, wh, alpha);
    }

    float scaleAlpha(const Vector3f &wi) const {
        return (1.2f-.2f*sqrtf(Frame::absCosTheta(wi))) * m_alpha;
    }

    float evalR(const Vector3f &wi, const Vector3f &wo) const {
        if (!sameHemisphere(wi, wo)) return 0;

        float absCosThetaI = std::fabsf(Frame::cosTheta(wi));
        float absCosThetaO = std::fabsf(Frame::cosTheta(wo));
        if (absCosThetaI == 0.f || absCosThetaO == 0.f) return 0;

        Vector3f wh = (wi + wo).normalized();
        if (wh.z() < 0) wh = -wh;
        float absCosThetaH = std::fabsf(Frame::cosTheta(wh));
        if (absCosThetaH == 0.f) return 0;

        float D = Warp::squareToBeckmannPdf(wh, scaleAlpha(wi));
        float F = fresnel(wh.dot(wi), m_extIOR, m_intIOR);
        return D * F * G(wi, wo, wh) / (4 * absCosThetaI * absCosThetaO * absCosThetaH);
    }


    float pdfR(const Vector3f &wi, const Vector3f &wo) const {
        if (!sameHemisphere(wo, wi)) return 0;
        Vector3f wh = (wo + wi).normalized();
        if (wh.z() < 0) wh = -wh;

        float F = fresnel(wh.dot(wi), m_extIOR, m_intIOR);
        float D = Warp::squareToBeckmannPdf(wh, scaleAlpha(wi));
        return F * D / (4 * fabs(wo.dot(wh)));
    }

    float evalT(const Vector3f &wi, const Vector3f &wo) const {
        if (sameHemisphere(wi, wo)) return 0;
        if ((wo.array() == 0).all()) return 0;

        float absCosThetaI = std::fabsf(Frame::cosTheta(wi));
        float absCosThetaO = std::fabsf(Frame::cosTheta(wo));
        if (absCosThetaI == 0.f || absCosThetaO == 0.f) return 0;

        float eta = Frame::cosTheta(wi) > 0 ? m_intIOR / m_extIOR : m_extIOR / m_intIOR;
        Vector3f wh = -(eta * wo + wi).normalized();
        if (wh.z() < 0) wh = -wh;
        if (wo.dot(wh) * wi.dot(wh) > 0) return 0;

        float absCosThetaH = std::fabsf(Frame::cosTheta(wh));
        if (absCosThetaH == 0.f) return 0;

        float D = Warp::squareToBeckmannPdf(wh, scaleAlpha(wi)) / absCosThetaH;
        float F = fresnel(wh.dot(wi), m_extIOR, m_intIOR);
        float sqrtDenom = eta * wo.dot(wh) + wi.dot(wh);
        float dwh_dwo_dwi = std::fabsf(wo.dot(wh)) * std::fabsf(wi.dot(wh)) / (sqrtDenom * sqrtDenom * absCosThetaI * absCosThetaO);
        return (1.f - F) * D * G(wi, wo, wh) * dwh_dwo_dwi;
    }


    float pdfT(const Vector3f &wi, const Vector3f &wo) const {
        if (sameHemisphere(wo, wi)) return 0;

        float eta = Frame::cosTheta(wi) > 0 ? m_intIOR / m_extIOR : m_extIOR / m_intIOR;
        Vector3f wh = -(eta * wo + wi).normalized();
        if (wo.dot(wh) * wi.dot(wh) > 0) return 0;

        float sqrtDenom = eta * wo.dot(wh) + wi.dot(wh);
        float dwh_dwo = eta * eta * std::fabsf(wo.dot(wh)) / (sqrtDenom * sqrtDenom);
        float F = fresnel(wh.dot(wi), m_extIOR, m_intIOR);
        float D = Warp::squareToBeckmannPdf(wh.z() > 0.f ? wh : -wh, scaleAlpha(wi));
        return (1-F) * D * dwh_dwo;
    }

    float m_alpha;
    float m_intIOR, m_extIOR;
};

NORI_REGISTER_CLASS(MicrofacetDielectric, "microfacetdielectric");
NORI_NAMESPACE_END

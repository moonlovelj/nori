/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
*/

#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

/// Ideal dielectric BSDF
class Dielectric : public BSDF {
public:
    Dielectric(const PropertyList &propList) {
        /* Interior IOR (default: BK7 borosilicate optical glass) */
        m_intIOR = propList.getFloat("intIOR", 1.5046f);

        /* Exterior IOR (default: air) */
        m_extIOR = propList.getFloat("extIOR", 1.000277f);
    }

    Color3f eval(const BSDFQueryRecord &) const {
        /* Discrete BRDFs always evaluate to zero in Nori */
        return Color3f(0.0f);
    }

    float pdf(const BSDFQueryRecord &) const {
        /* Discrete BRDFs always evaluate to zero in Nori */
        return 0.0f;
    }

    Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample) const {
	  bRec.measure = EDiscrete;
	  auto cosThetaI = Frame::cosTheta(bRec.wi);
	  if (cosThetaI < 0) {
		bRec.eta = m_intIOR / m_extIOR;
	  } else {
		bRec.eta = m_extIOR / m_intIOR;
	  }
	  auto reflectance = fresnel(cosThetaI, m_extIOR, m_intIOR);
	  if (sample.x() <= reflectance) {
//		auto reflect = Vector3f(-bRec.wi.x(),-bRec.wi.y(),bRec.wi.z());
//		auto wo = Warp::squareToCosineHemisphere(sample);
//		auto cosThetaOAndR = Frame::cosTheta(wo);
//		auto pdf = Warp::squareToCosineHemispherePdf(wo);
//		Frame shReflect(reflect);
//		Frame shNormal(Vector3f(0,0,1));
//		wo = shReflect.toLocal(shReflect.toWorld(wo));
//		bRec.wo = wo;
//		auto cosThetaO = cosThetaI < 0 ? -Frame::cosTheta(bRec.wo) : Frame::cosTheta(bRec.wo);
//		return 1.5 * M_PI * cosThetaOAndR * std::max(0.f, cosThetaO) / pdf / reflectance;
		bRec.wo = Vector3f(-bRec.wi.x(),-bRec.wi.y(),bRec.wi.z());
		return Color3f(1);
	  } else {
		bRec.wo = refract(Vector3f(0,0,1), bRec.wi, m_extIOR, m_intIOR);
		return Color3f(1);
	  }
    }

    std::string toString() const {
        return tfm::format(
            "Dielectric[\n"
            "  intIOR = %f,\n"
            "  extIOR = %f\n"
            "]",
            m_intIOR, m_extIOR);
    }
private:
    float m_intIOR, m_extIOR;
};

NORI_REGISTER_CLASS(Dielectric, "dielectric");
NORI_NAMESPACE_END

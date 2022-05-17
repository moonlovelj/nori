/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
*/

#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class Microfacet : public BSDF {
public:
    Microfacet(const PropertyList &propList) {
        /* RMS surface roughness */
        m_alpha = propList.getFloat("alpha", 0.1f);

        /* Interior IOR (default: BK7 borosilicate optical glass) */
        m_intIOR = propList.getFloat("intIOR", 1.5046f);

        /* Exterior IOR (default: air) */
        m_extIOR = propList.getFloat("extIOR", 1.000277f);

        /* Albedo of the diffuse base material (a.k.a "kd") */
        m_kd = propList.getColor("kd", Color3f(0.5f));

        /* To ensure energy conservation, we must scale the 
           specular component by 1-kd. 

           While that is not a particularly realistic model of what 
           happens in reality, this will greatly simplify the 
           implementation. Please see the course staff if you're 
           interested in implementing a more realistic version 
           of this BRDF. */
        m_ks = 1 - m_kd.maxCoeff();
    }

    /// Evaluate the BRDF for the given pair of directions
    Color3f eval(const BSDFQueryRecord &bRec) const {
    	if (Frame::cosTheta(bRec.wi) < 0) return Color3f(0);

    	Vector3f wh = (bRec.wi + bRec.wo).normalized();
    	float D = Warp::squareToBeckmannPdf(wh, m_alpha);
    	float F = fresnel(wh.dot(bRec.wi), m_extIOR, m_intIOR);
		float G = GTerm(bRec.wi, wh) * GTerm(bRec.wo, wh);
		return m_kd * M_1_PI + Color3f(
			m_ks * D * F * G / (4 * Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) * Frame::cosTheta(wh)));
    }

    /// Evaluate the sampling density of \ref sample() wrt. solid angles
	float pdf(const BSDFQueryRecord &bRec) const {
	  if (bRec.measure != ESolidAngle) {
		return 0.f;
	  }
//
	  Vector3f wh = (bRec.wi + bRec.wo).normalized();
//	  if (bRec.wi.dot(wh) / Frame::cosTheta(bRec.wi) <= 0 ||
//		  bRec.wo.dot(wh) / Frame::cosTheta(bRec.wo) <= 0) {
//		return 0.f;
//	  }

	  if (Frame::cosTheta(bRec.wi) <= 0 ||
		  Frame::cosTheta(bRec.wo) <= 0) {
		return 0.f;
	  }

	  float D = Warp::squareToBeckmannPdf(wh, m_alpha);
	  return m_ks * D / (4 * wh.dot(bRec.wo)) + (1 - m_ks) * Frame::cosTheta(bRec.wo) * M_1_PI;
	}

    /// Sample the BRDF
    Color3f sample(BSDFQueryRecord &bRec, const Point2f &_sample) const {
	  	if (Frame::cosTheta(bRec.wi) <= 0)
			return Color3f(0.0f);
        // Note: Once you have implemented the part that computes the scattered
        // direction, the last part of this function should simply return the
        // BRDF value divided by the solid angle density and multiplied by the
        // cosine factor from the reflection equation, i.e.
        // return eval(bRec) * Frame::cosTheta(bRec.wo) / pdf(bRec);
	    bRec.measure = ESolidAngle;
        if (_sample.x() < m_ks) {
		  Point2f reuseSample(_sample.x() / m_ks, _sample.y());
		  Vector3f h = Warp::squareToBeckmann(reuseSample, m_alpha);
		  bRec.wo = reflect(bRec.wi, h);
		  if (Frame::cosTheta(bRec.wo) <= 0) return Color3f(0);
        } else {
          //diffuse
		  Point2f reuseSample((_sample.x()-m_ks)/(1-m_ks), _sample.y());
		  bRec.wo = Warp::squareToCosineHemisphere(reuseSample);
        }
	  	return eval(bRec) * Frame::cosTheta(bRec.wo) / pdf(bRec);
    }

    bool isDiffuse() const {
        /* While microfacet BRDFs are not perfectly diffuse, they can be
           handled by sampling techniques for diffuse/non-specular materials,
           hence we return true here */
        return true;
    }

    std::string toString() const {
        return tfm::format(
            "Microfacet[\n"
            "  alpha = %f,\n"
            "  intIOR = %f,\n"
            "  extIOR = %f,\n"
            "  kd = %s,\n"
            "  ks = %f\n"
            "]",
            m_alpha,
            m_intIOR,
            m_extIOR,
            m_kd.toString(),
            m_ks
        );
    }
private:
  	float GTerm(const Vector3f &wv, const Vector3f &wh) const {
		float cosThetaVN = Frame::cosTheta(wv);
		float sinThetaVN = std::sqrtf(1 - cosThetaVN * cosThetaVN);
		float tanThetaVN = sinThetaVN / cosThetaVN;
		float b = 1.f / (m_alpha * tanThetaVN);
		float c = wv.dot(wh) / cosThetaVN;
		float Xc = c > 0 ? 1 : 0;
		return Xc * (b < 1.6f ? (3.535f*b+2.181f*b*b)/(1+2.276f*b+2.577*b*b):1);
    }

    float m_alpha;
    float m_intIOR, m_extIOR;
    float m_ks;
    Color3f m_kd;
};

NORI_REGISTER_CLASS(Microfacet, "roughplastic");
NORI_NAMESPACE_END

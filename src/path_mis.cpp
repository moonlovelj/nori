//
// Created by 郭彬 on 2021/12/1.
//
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/bsdf.h>
#include <nori/sampler.h>

NORI_NAMESPACE_BEGIN

class PathMultipleImportanceSamplingIntegrator : public Integrator {
public:
  PathMultipleImportanceSamplingIntegrator(const PropertyList &props) {
	/* No parameters this time */
  }

  Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
    /// 先选择一种采样策略
	Intersection its;
	if (!scene->rayIntersect(ray, its)) {
	  return Color3f(0);
	}

	float nLight, nBRDF;
	if (its.mesh->getBSDF()->isDiffuse()) {
	  nLight = 0.5f;
	  nBRDF = 0.5f;
	} else {
	  nLight = 0;
	  nBRDF = 1;
	}

	return Color3f(0);
  }

  std::string toString() const {
	return "PathMultipleImportanceSamplingIntegrator[]";
  }

private:
  float BalanceHeuristic(int nf, float fPdf, int ng, float gPdf) const {
	return (nf * fPdf) / (nf * fPdf + ng * gPdf);
  }
};

NORI_REGISTER_CLASS(PathMultipleImportanceSamplingIntegrator, "path_mis");
NORI_NAMESPACE_END
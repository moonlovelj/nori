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
    return Color3f(0);
  }

  std::string toString() const {
	return "PathMultipleImportanceSamplingIntegrator[]";
  }
};

NORI_REGISTER_CLASS(PathMultipleImportanceSamplingIntegrator, "path_mis");
NORI_NAMESPACE_END
//
// Created by 郭彬 on 2021/11/22.
//


#include <nori/integrator.h>
#include <nori/scene.h>

NORI_NAMESPACE_BEGIN

class WhittedIntegrator : public Integrator {
 public:
  WhittedIntegrator(const PropertyList &props) {
	  /* No parameters this time */
  }

  Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
	  /* Find the surface that is visible in the requested direction */
	  Intersection its;
	  if (!scene->rayIntersect(ray, its))
		  return Color3f(0.0f);

	  return Color3f(0);
  }

  std::string toString() const {
	  return "WhittedIntegrator[]";
  }
};

NORI_REGISTER_CLASS(WhittedIntegrator, "whitted");
NORI_NAMESPACE_END
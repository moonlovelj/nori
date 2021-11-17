
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class AO : public Integrator {
 public:
  AO(const PropertyList &props) {
	  /* No parameters this time */
  }

  Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
	  /* Find the surface that is visible in the requested direction */
	  Intersection its;
	  if (!scene->rayIntersect(ray, its))
		  return Color3f(0.0f);

	  auto sample_dir = Warp::squareToCosineHemisphere(sampler->next2D());
	  auto dir = its.shFrame.toWorld(sample_dir).normalized();
	  auto origin = its.p + Epsilon * its.shFrame.n;
	  auto visibility_ray = Ray3f(origin, dir);
	  if (scene->rayIntersect(visibility_ray))
		  return Color3f(0.0f);

	  auto cos_theta = sample_dir.z();
	  auto pdf = Warp::squareToCosineHemispherePdf(sample_dir);
	  return Color3f(cos_theta * M_1_PI / pdf);
  }

  std::string toString() const {
	  return "AO[]";
  }
};

NORI_REGISTER_CLASS(AO, "ao");
NORI_NAMESPACE_END
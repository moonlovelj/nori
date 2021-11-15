//
// Created by 郭彬 on 2021/11/15.
//

#include <nori/integrator.h>
#include <nori/scene.h>

NORI_NAMESPACE_BEGIN

class SimpleIntegrator : public Integrator {
 public:
  SimpleIntegrator(const PropertyList &props){
	  this->position_ = props.getPoint("position");
	  this->intensity_ = props.getColor("energy");
  }

  Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
	  /* Find the surface that is visible in the requested direction */
	  Intersection its;
	  if (!scene->rayIntersect(ray, its))
		  return Color3f(0.0f);

	  auto origin = its.p + Epsilon * its.shFrame.n;
	  Ray3f shadow_ray(origin, (position_-origin).normalized());
	  if (scene->rayIntersect(shadow_ray))
	  	  return Color3f(0.0f);

	  auto dir = position_-its.p;
	  auto cos_theta = std::max(0.f, its.shFrame.n.dot(dir.normalized()));
	  auto color = intensity_ * cos_theta * M_1_PI * M_1_PI * 0.25f / dir.squaredNorm();
	  return Color3f(color);
  }

  std::string toString() const {
	  return "SimpleIntegrator[]";
  }

 private:
  	Color3f intensity_;
  	Point3f position_;
};

NORI_REGISTER_CLASS(SimpleIntegrator, "simple");
NORI_NAMESPACE_END
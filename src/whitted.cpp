//
// Created by 郭彬 on 2021/11/22.
//


#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/bsdf.h>
#include <nori/sampler.h>

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

	Color3f l_e = its.mesh->getEmission(its, -ray.d);
	Color3f l_dir(0);
	auto lights = scene->getEmitters();
	Point3f point_on_light;
	Vector3f normal_from_light;
	float pdf;
	auto l_i = lights[std::rand() % lights.size()]->SampleLight(point_on_light, normal_from_light, pdf, sampler->next2D());
	auto dis = point_on_light - its.p;
	auto wi = dis.normalized();
	auto shadow_ray = Ray3f(its.p, wi);
	shadow_ray.mint = Epsilon;
	shadow_ray.maxt = dis.norm() - Epsilon;
	if (scene->rayIntersect(shadow_ray)) {
	  return l_e;
	}

	BSDFQueryRecord bsdf_record(its.shFrame.toLocal(wi), its.shFrame.toLocal(-ray.d), ESolidAngle);
	l_dir = l_i * its.mesh->getBSDF()->eval(bsdf_record) * std::max(0.f, its.shFrame.n.dot(wi))
		* std::max(0.f, normal_from_light.dot(-wi)) / dis.squaredNorm() / pdf * lights.size();

	return l_e + l_dir;
  }

  std::string toString() const {
	return "WhittedIntegrator[]";
  }
};

NORI_REGISTER_CLASS(WhittedIntegrator, "whitted");
NORI_NAMESPACE_END
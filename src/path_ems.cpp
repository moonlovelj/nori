//
// Created by 郭彬 on 2021/11/30.
//

#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/bsdf.h>
#include <nori/sampler.h>

NORI_NAMESPACE_BEGIN

class PathEmitterSamplingIntegrator : public Integrator {
public:
  PathEmitterSamplingIntegrator(const PropertyList &props) {
	/* No parameters this time */
  }

  Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
	return Li(scene, sampler, ray, true);
  }

  std::string toString() const {
	return "PathEmitterSamplingIntegrator[]";
  }
private:
  Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray, bool includeEmitter) const {
	Intersection its;
	if (!scene->rayIntersect(ray, its)) {
	  return Color3f(0);
	}

	Color3f L_e(0.f);
	if (includeEmitter) {
	  L_e = its.mesh->getEmission(its, -ray.d);
	}

	if (sampler->next1D() > 0.95f) {
	  return L_e;
	}

	Color3f L_dir(0);
	if(its.mesh->getBSDF()->isDiffuse()) {
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
	  if (!scene->rayIntersect(shadow_ray)) {
		BSDFQueryRecord bsdf_record(its.shFrame.toLocal(wi), its.shFrame.toLocal(-ray.d), ESolidAngle);
		L_dir = l_i * its.mesh->getBSDF()->eval(bsdf_record) * std::max(0.f, its.shFrame.n.dot(wi))
			* std::max(0.f, normal_from_light.dot(-wi)) / dis.squaredNorm() / pdf * lights.size() / 0.95f;
	  }
	}

	BSDFQueryRecord bsdfRec(its.shFrame.toLocal(-ray.d));
	Color3f L = its.mesh->getBSDF()->sample(bsdfRec, sampler->next2D());
	Ray3f new_ray(its.p, its.shFrame.toWorld(bsdfRec.wo));
	new_ray.mint = Epsilon;
	Color3f l_ind = Li(scene, sampler, new_ray, its.mesh->getBSDF()->isDiffuse() ? false : true) * L / 0.95f;
	return L_e + L_dir + l_ind;
  }
};

NORI_REGISTER_CLASS(PathEmitterSamplingIntegrator, "path_ems");
NORI_NAMESPACE_END
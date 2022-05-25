//
// Created by 郭彬 on 2021/11/30.
//

#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/bsdf.h>
#include <nori/sampler.h>
#include <nori/emitter.h>

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
	  Emitter *pLight = lights[std::rand() % lights.size()]->getEmitter();
	  EmitterQueryRecord eRec;
	  Color3f l_i = pLight->sample(its.p, eRec, sampler->next2D());
	  Vector3f wi = (eRec.point - its.p).normalized();
	  if (scene->illuminatedEachOther(its.p, eRec.point)) {
		BSDFQueryRecord sampleLightRecord(its.shFrame.toLocal(-ray.d), its.shFrame.toLocal(wi),  ESolidAngle, sampler);
		L_dir = l_i * its.mesh->getBSDF()->eval(sampleLightRecord) * std::max(0.f, its.shFrame.n.dot(wi)) * lights.size() / 0.95f;
	  }
	}

	BSDFQueryRecord sampleBRDFRecord(its.shFrame.toLocal(-ray.d), sampler);
	Color3f L = its.mesh->getBSDF()->sample(sampleBRDFRecord);
	Ray3f new_ray(its.p, its.shFrame.toWorld(sampleBRDFRecord.wo));
	new_ray.mint = Epsilon;
	Color3f l_ind = Li(scene, sampler, new_ray, its.mesh->getBSDF()->isDiffuse() ? false : true) * L / 0.95f;
	return L_e + L_dir + l_ind;
  }
};

NORI_REGISTER_CLASS(PathEmitterSamplingIntegrator, "path_ems");
NORI_NAMESPACE_END
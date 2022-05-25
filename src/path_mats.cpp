//
// Created by 郭彬 on 2021/11/30.
//

#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/bsdf.h>
#include <nori/sampler.h>

NORI_NAMESPACE_BEGIN

class PathMaterialSamplingIntegrator : public Integrator {
public:
    PathMaterialSamplingIntegrator(const PropertyList &props) {
        /* No parameters this time */
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        /* Find the surface that is visible in the requested direction */
        Intersection its;
        if (!scene->rayIntersect(ray, its)) {
            return Color3f(0);
        }

        Color3f sampleColor = its.mesh->getEmission(its, -ray.d);
        if (sampler->next1D() > 0.95f) {
            return sampleColor;
        }

        BSDFQueryRecord bsdfRec(its.shFrame.toLocal(-ray.d), sampler);
        Color3f L = its.mesh->getBSDF()->sample(bsdfRec);
        Ray3f new_ray(its.p, its.shFrame.toWorld(bsdfRec.wo));
        new_ray.mint = Epsilon;
        return sampleColor + Li(scene, sampler, new_ray) * L / 0.95f;

//	Color3f L(0.f);
//	Ray3f sampleRay(ray);
//	int bounces = 0;
//	float continuation_probability = 1.f;
//	while (true) {
//	  if (bounces > 2) {
//	    //Russian Roulette
//		continuation_probability = 0.95f;
//		if (sampler->next1D() > continuation_probability){
//		  break;
//		}
//	  }
//	  ++bounces;
//	  if (!scene->rayIntersect(sampleRay, its)) {
//	    break;
//	  }
//
//	  Color3f sampleColor = its.mesh->getEmission(its, -sampleRay.d);
//	}
//	return L;
    }

    std::string toString() const {
        return "PathMaterialSamplingIntegrator[]";
    }
};

NORI_REGISTER_CLASS(PathMaterialSamplingIntegrator, "path_mats");
NORI_NAMESPACE_END
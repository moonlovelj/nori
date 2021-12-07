//
// Created by 郭彬 on 2021/11/22.
//


#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/bsdf.h>
#include <nori/sampler.h>
#include <nori/emitter.h>

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
        auto bsdf = its.mesh->getBSDF();
        if (!bsdf->isDiffuse()) {
            if (sampler->next1D() > 0.95f) return l_e;

            BSDFQueryRecord bsdf_query_record(its.shFrame.toLocal(-ray.d));
            Color3f sampleBSDF = bsdf->sample(bsdf_query_record, sampler->next2D()) / 0.95f;
            auto newRay = Ray3f(its.p, its.shFrame.toWorld(bsdf_query_record.wo).normalized());
            newRay.mint = Epsilon;
            auto next_sample = Li(scene, sampler, newRay);
            l_dir = sampleBSDF * next_sample;
            return l_e + l_dir;
        }

        auto lights = scene->getEmitters();
        EmitterQueryRecord sampleLightRecord;
        Emitter *pLight = lights[std::rand() % lights.size()]->getEmitter();
        auto l_i = pLight->sample(its.p, sampleLightRecord, sampler->next2D());
        auto wi = (sampleLightRecord.point - its.p).normalized();
        if (!scene->illuminatedEachOther(its.p, sampleLightRecord.point)) {
            return l_e;
        }
        BSDFQueryRecord bsdf_record(its.shFrame.toLocal(wi), its.shFrame.toLocal(-ray.d), ESolidAngle);
        l_dir = l_i * its.mesh->getBSDF()->eval(bsdf_record) * std::max(0.f, its.shFrame.n.dot(wi)) * lights.size();

        return l_e + l_dir;
    }

    std::string toString() const {
        return "WhittedIntegrator[]";
    }
};

NORI_REGISTER_CLASS(WhittedIntegrator, "whitted");
NORI_NAMESPACE_END
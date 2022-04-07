//
// Created by 郭彬 on 2021/11/30.
//

#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/bsdf.h>
#include <nori/sampler.h>
#include <nori/emitter.h>

NORI_NAMESPACE_BEGIN

class PathVolumeIntegrator : public Integrator {
public:
    PathVolumeIntegrator(const PropertyList &props) {
        /* No parameters this time */
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        return Li(scene, sampler, ray, true);
    }

    std::string toString() const {
        return "PathVolumeIntegrator[]";
    }

private:
    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray, bool includeEmitter) const {
        Intersection its;
        if (!scene->rayIntersect(ray, its)) {
            return Color3f(0);
        }

        Intersection newIts;
        Ray3f newRay(ray.o, ray.d, Epsilon, its.t);
        Color3f tr_pdf = scene->getMedium()->sample(newRay, sampler, newIts);
        if (newIts.isMedium()) {
            newIts.p = newRay(newIts.t);
            const Medium *medium = scene->getMedium();
            Point3f sampleMediumPoint = newIts.t * ray.d + ray.o;
            Color3f le = medium->getEmittance(sampleMediumPoint, -ray.d) * medium->sigmaA();
            if (sampler->next1D() > 0.95f) {
                tr_pdf * le;
            }

            Color3f L_dir = estimateDirect(newIts, -newRay.d, scene, sampler);
            Vector3f wo;
            Color3f phaseSample = medium->getPhase()->sample(-ray.d, wo, sampler->next2D());
            Ray3f nextRay(sampleMediumPoint, wo);
            return tr_pdf *
                   (le + medium->sigmaS() * (L_dir + phaseSample * Li(scene, sampler, nextRay, false)) / 0.95f);
        } else {
            Color3f le = includeEmitter ? its.mesh->getEmission(its, -ray.d) : Color3f(0);
            if (sampler->next1D() > 0.95f) {
                return tr_pdf * le;
            }

            BSDFQueryRecord sampleBRDFRecord(its.shFrame.toLocal(-ray.d));
            Color3f bsdf = its.mesh->getBSDF()->sample(sampleBRDFRecord, sampler->next2D());
            Ray3f nextRay(its.p, its.shFrame.toWorld(sampleBRDFRecord.wo));

            if (its.mesh->getBSDF()->isDiffuse()) {
                Color3f L_dir = estimateDirect(its, -ray.d, scene, sampler);
                return tr_pdf * (le + (L_dir + bsdf * Li(scene, sampler, nextRay, false)) / 0.95f);
            } else {
                return tr_pdf * (le + bsdf * Li(scene, sampler, nextRay, true) / 0.95f);
            }
        }
    }
};

NORI_REGISTER_CLASS(PathVolumeIntegrator, "path_volume");
NORI_NAMESPACE_END
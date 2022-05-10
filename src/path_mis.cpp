//
// Created by 郭彬 on 2021/12/1.
//
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/bsdf.h>
#include <nori/sampler.h>
#include <nori/emitter.h>

NORI_NAMESPACE_BEGIN

class PathMultipleImportanceSamplingIntegrator : public Integrator {
public:
    PathMultipleImportanceSamplingIntegrator(const PropertyList &props) {
        /* No parameters this time */
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        return Li(scene, sampler, ray, true);
    }

    std::string toString() const {
        return "PathMultipleImportanceSamplingIntegrator[]";
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
        Color3f l_ind(0);
        if (its.mesh->getBSDF()->isDiffuse()) {
            float sampleLightProbability = 0.5f;
            if (sampler->next1D() < sampleLightProbability) {
                // sample light
                auto lights = scene->getEmitters();
                Emitter *pLight = lights[std::rand() % lights.size()]->getEmitter();
                EmitterQueryRecord eRec;
                pLight->sample(its.p, eRec, sampler->next2D());
                if (scene->illuminatedEachOther(its.p, eRec.point)) {
                    Vector3f wi = (eRec.point - its.p).normalized();
                    float pdfLight =
                            pLight->pdf(eRec) * (eRec.point - its.p).squaredNorm() / std::fabsf(eRec.normal.dot(-wi));
                    BSDFQueryRecord sampleLightRecord(its.shFrame.toLocal(-ray.d), its.shFrame.toLocal(wi),
                                                      ESolidAngle);
                    sampleLightRecord.sampler = sampler;
                    float pdfBSDF = its.mesh->getBSDF()->pdf(sampleLightRecord);
                    L_dir = pLight->eval(eRec) *
                            its.mesh->getBSDF()->eval(sampleLightRecord) * std::max(0.f, its.shFrame.n.dot(wi))
                            * lights.size() / 0.95f /
                            (sampleLightProbability * pdfLight + (1 - sampleLightProbability) * pdfBSDF);
                }
            } else {
                // sample bsdf
                BSDFQueryRecord sampleBRDFRecord(its.shFrame.toLocal(-ray.d));
                sampleBRDFRecord.sampler = sampler;
                its.mesh->getBSDF()->sample(sampleBRDFRecord, sampler);
                float pdfBSDF = its.mesh->getBSDF()->pdf(sampleBRDFRecord);
                float pdfLight = 0.f;
                Ray3f nextRay(its.p, its.shFrame.toWorld(sampleBRDFRecord.wo), Epsilon,
                              std::numeric_limits<float>::infinity());
                Intersection itsNext;
                if (scene->rayIntersect(nextRay, itsNext) && itsNext.mesh->isEmitter()) {
                    Vector3f wi = (itsNext.p - its.p).normalized();
                    pdfLight = itsNext.mesh->getEmitter()->pdf(EmitterQueryRecord(itsNext.p, itsNext.shFrame.n)) *
                               (itsNext.p - its.p).squaredNorm() / std::fabsf(itsNext.shFrame.n.dot(-wi));
                    L_dir = itsNext.mesh->getEmitter()->eval(EmitterQueryRecord(itsNext.p, itsNext.shFrame.n)) *
                            std::max(0.f, its.shFrame.n.dot(nextRay.d)) *
                            its.mesh->getBSDF()->eval(sampleBRDFRecord) / 0.95f /
                            (sampleLightProbability * pdfLight + (1 - sampleLightProbability) * pdfBSDF);

                }
            }
            BSDFQueryRecord sampleIndirectBSDF(its.shFrame.toLocal(-ray.d));
            sampleIndirectBSDF.sampler = sampler;
            l_ind = its.mesh->getBSDF()->sample(sampleIndirectBSDF, sampler) *
                    Li(scene, sampler, Ray3f(its.p, its.shFrame.toWorld(sampleIndirectBSDF.wo), Epsilon,
                                             std::numeric_limits<float>::infinity()), false) / 0.95f;

        } else {
            BSDFQueryRecord sampleBRDFRecord(its.shFrame.toLocal(-ray.d));
            sampleBRDFRecord.sampler = sampler;
            Color3f L = its.mesh->getBSDF()->sample(sampleBRDFRecord, sampler);
            Ray3f new_ray(its.p, its.shFrame.toWorld(sampleBRDFRecord.wo));
            new_ray.mint = Epsilon;
            l_ind = Li(scene, sampler, new_ray, true) * L / 0.95f;
        }

        return L_e + L_dir + l_ind;
    }
};

NORI_REGISTER_CLASS(PathMultipleImportanceSamplingIntegrator, "path_mis");
NORI_NAMESPACE_END
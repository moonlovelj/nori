//
// Created by 郭彬 on 2022/4/7.
//

#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>
#include <nori/sampler.h>

NORI_NAMESPACE_BEGIN

Color3f Integrator::estimateDirect(const Intersection &its,
                                   const Vector3f &w, const Scene *scene, Sampler *sampler) const {
    Color3f L_dir(0);
    auto lights = scene->getEmitters();
    Emitter *pLight = lights[std::rand() % lights.size()]->getEmitter();
    EmitterQueryRecord eRec;
    Color3f l_i = pLight->sample(its.p, eRec, sampler->next2D()) * lights.size();
    Vector3f wi = (eRec.point - its.p).normalized();

    if (!scene->illuminatedEachOther(its.p, eRec.point)) {
        return L_dir;
    }

    if (its.isMedium()) {
        L_dir = l_i * its.getMedium()->tr(Ray3f(its.p,wi, Epsilon, (eRec.point - its.p).norm()), sampler) *
                its.getMedium()->getPhase()->p(w, wi);
    } else {
        BSDFQueryRecord sampleLightRecord(its.shFrame.toLocal(w), its.shFrame.toLocal(wi), ESolidAngle);
        L_dir = l_i * its.mesh->getBSDF()->eval(sampleLightRecord) * std::max(0.f, its.shFrame.n.dot(wi));
    }
    return L_dir;
}

NORI_NAMESPACE_END
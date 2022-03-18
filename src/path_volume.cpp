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
        return Color3f(0);
    }

    std::string toString() const {
        return "PathVolumeIntegrator[]";
    }
};

NORI_REGISTER_CLASS(PathVolumeIntegrator, "path_volume");
NORI_NAMESPACE_END
//
// Created by 郭彬 on 2022/3/31.
//

#include <nori/medium.h>
#include <nori/sampler.h>
#include <nori/phasefunction.h>
#include <nori/color.h>

NORI_NAMESPACE_BEGIN

class GridDensityMedium : public Medium {
public:

    GridDensityMedium(const PropertyList &propList) {
        sigma_a = propList.getColor("sigma_a", Color3f(0.1f));
        sigma_s = propList.getColor("sigma_s", Color3f(0.2f));
        g = propList.getFloat("g", 0.1f);
        sigma_t = sigma_s + sigma_a;
        m_phase = std::make_shared<HenyeyGreenstein>(g);
    }

    Color3f tr(const Ray3f &ray, Sampler *sampler) const override {
        return Color3f(1);
    }

    Color3f sample(const Ray3f &ray, Sampler *sampler, Intersection &its) const override {
//        int channel = std::min((int) (sampler->next1D() * 3), 2);
//        float dist = -std::log(1 - sampler->next1D()) / sigma_t[channel];
//        float t = std::min(dist * ray.d.norm(), ray.maxt);
//        its.t = t;
//        bool sampledMedium = t < ray.maxt;
//        if (sampledMedium) {
//            its.mediumInterface = this;
//            its.insideMedium = true;
//        }
//
//        Color3f e = -sigma_t * t;
//        Color3f transmission(std::exp(e.x()), std::exp(e.y()), std::exp(e.z()));
//
//        Color3f density = sampledMedium ? (sigma_t * transmission) : transmission;
//        float pdf = density.x() + density.y() + density.z();
//        pdf /= 3;
//        return (Color3f) (transmission / pdf);
        return Color3f(1);
    }

    Color3f sigmaA() const override {
        return sigma_a;
    }

    Color3f sigmaS() const override {
        return sigma_s;
    }

    std::string toString() const override {
        return tfm::format("GridDensityMedium, sigma_a:%s, sigma_s:%s, g:%s \n", sigma_a, sigma_s, g);
    }

private:
    Color3f sigma_a, sigma_s, sigma_t;
    float g;
};

NORI_REGISTER_CLASS(GridDensityMedium, "griddensitymedium");
NORI_NAMESPACE_END

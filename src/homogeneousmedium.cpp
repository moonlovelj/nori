//
// Created by 郭彬 on 2022/1/19.
//

#include <nori/medium.h>
#include <nori/sampler.h>

NORI_NAMESPACE_BEGIN

class HenyeyGreenstein;

class HomogeneousMedium : public Medium {
public:
    HomogeneousMedium(const Color3f &sigma_a, const Color3f &sigma_s,
                      float g)
            : sigma_a(sigma_a), sigma_s(sigma_s), sigma_t(sigma_s + sigma_a),
              g(g) {}

    Color3f tr(const Ray3f &ray, Sampler *sampler) const override {
        Color3f exponent = -sigma_t * std::min(ray.maxt * ray.d.norm(), std::numeric_limits<float>::max());
        return Color3f(std::exp(exponent.x()), std::exp(exponent.y()), std::exp(exponent.z()));
    }

    Color3f sample(const Ray3f &ray, Sampler *sampler, Intersection &its) const override {
        int channel = std::min((int)(sampler->next1D() * 3), 2);
        float dist = -std::log(1 - sampler->next1D()) / sigma_t[channel];
        float t = std::min(dist * ray.d.norm(), ray.maxt);
        bool sampledMedium = t < ray.maxt;
        if (sampledMedium) {
            its = Intersection();
            its.t = t;
            its.mediumInterface = MediumInterface(this);
            its.phase = std::make_shared<HenyeyGreenstein>(g);
        }

        Color3f e = -sigma_t * std::min(t, std::numeric_limits<float>::max()) * ray.d.norm();
        Color3f transmission(std::exp(e.x()), std::exp(e.y()), std::exp(e.z()));

        Color3f density = sampledMedium ? (sigma_t * transmission) : transmission;
        float pdf = density.x() + density.y() + density.z();
        pdf /= 3;
        return sampledMedium ? ((Color3f)(transmission * sigma_s) / pdf) : (transmission / pdf);
    }

private:
    const Color3f sigma_a, sigma_s, sigma_t;
    const float g;
};

NORI_NAMESPACE_END



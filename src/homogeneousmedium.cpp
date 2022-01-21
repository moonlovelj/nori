//
// Created by 郭彬 on 2022/1/19.
//

#include <nori/medium.h>

NORI_NAMESPACE_BEGIN

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

    }

private:
    const Color3f sigma_a, sigma_s, sigma_t;
    const float g;
};

NORI_NAMESPACE_END



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
        if (sigma_t.x() != sigma_t.y() || sigma_t.x() != sigma_t.z()) {
            std::cout << "GridDensityMedium requires a spectrally uniform attenuation \n";
        }

        m_phase = std::make_shared<HenyeyGreenstein>(g);
        worldToMedium = propList.getTransform("toWorld", Transform()).inverse();
        invMaxDensity = 1;
    }

    Color3f tr(const Ray3f &ray, Sampler *sampler) const override {
        return Color3f(1);
    }

    Color3f sample(const Ray3f &ray, Sampler *sampler, Intersection &its) const override {
        Ray3f localRay = worldToMedium * ray;
        const BoundingBox3f b(Point3f(0, 0, 0), Point3f(1, 1, 1));
        float tMin, tMax;
        if (!b.rayIntersect(localRay, tMin, tMax)) return Color3f(1);

        float t = tMin;
        while (true) {
            t -= std::log(1.f - sampler->next1D()) * invMaxDensity / sigma_t[0];
            if (t >= tMax) break;
            if (density(localRay(t)) * invMaxDensity > sampler->next1D()) {
                its.mediumInterface = this;
                its.insideMedium = true;
                its.t = t;
                return Color3f(1.f / sigma_t[0]);
            }
        }

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
    float invMaxDensity;
    int nx, ny, nz;
    Transform worldToMedium;

    float density(const Point3f &p) const {
        // Compute voxel coordinates and offsets for _p_
//        Point3f pSamples(p.x * nx - .5f, p.y * ny - .5f, p.z * nz - .5f);
//        Point3i pi = (Point3i)Floor(pSamples);
//        Vector3f d = pSamples - (Point3f)pi;
//
//        // Trilinearly interpolate density values to compute local density
//        Float d00 = Lerp(d.x, D(pi), D(pi + Vector3i(1, 0, 0)));
//        Float d10 = Lerp(d.x, D(pi + Vector3i(0, 1, 0)), D(pi + Vector3i(1, 1, 0)));
//        Float d01 = Lerp(d.x, D(pi + Vector3i(0, 0, 1)), D(pi + Vector3i(1, 0, 1)));
//        Float d11 = Lerp(d.x, D(pi + Vector3i(0, 1, 1)), D(pi + Vector3i(1, 1, 1)));
//        Float d0 = Lerp(d.y, d00, d10);
//        Float d1 = Lerp(d.y, d01, d11);
//        return Lerp(d.z, d0, d1);
          return 1.f;
    }
};

NORI_REGISTER_CLASS(GridDensityMedium, "griddensitymedium");
NORI_NAMESPACE_END

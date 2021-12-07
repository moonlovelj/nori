//
// Created by 郭彬 on 2021/11/22.
//

#include <nori/emitter.h>
#include <nori/mesh.h>
#include <nori/dpdf.h>
#include <Eigen/Dense>

NORI_NAMESPACE_BEGIN

class AreaLight : public Emitter {
public:
    AreaLight(const PropertyList &props) {
        m_radiance = props.getColor("radiance", Color3f(1.f));
    }

    void setParent(NoriObject *parent) override {
        m_parent = parent;
    }

    void build() override {
        Mesh *pMesh = dynamic_cast<Mesh *>(m_parent);
        auto face_count = pMesh->getTriangleCount();
        m_DPDF = std::make_shared<DiscretePDF>(face_count);
        for (uint32_t face_index = 0; face_index < face_count; ++face_index) {
            m_DPDF->append(pMesh->surfaceArea(face_index));
        }
        m_DPDF->normalize();
    }

    Color3f eval(const EmitterQueryRecord &eRec) const override {
        return m_radiance;
    }

    Color3f sample(const Point3f &surfacePoint, EmitterQueryRecord &eRec, const Point2f &sample) const override {
        Mesh *pMesh = dynamic_cast<Mesh *>(m_parent);
        float eps1 = sample.x();
        float esp2 = sample.y();
        auto face_index = m_DPDF->sampleReuse(eps1);
        auto alpha = 1.f - std::sqrtf(1.f - eps1);
        auto beta = esp2 * std::sqrtf(1.f - eps1);
        auto gamma = 1.f - alpha - beta;

        auto i0 = pMesh->getIndices()(0, face_index);
        auto i1 = pMesh->getIndices()(1, face_index);
        auto i2 = pMesh->getIndices()(2, face_index);
        Point3f p0 = pMesh->getVertexPositions().col(i0);
        Point3f p1 = pMesh->getVertexPositions().col(i1);
        Point3f p2 = pMesh->getVertexPositions().col(i2);
        eRec.point = alpha * p0 + beta * p1 + gamma * p2;
        if (pMesh->getVertexNormals().size() > 0) {
            eRec.normal = (alpha * pMesh->getVertexNormals().col(i0) +
                           beta * pMesh->getVertexNormals().col(i1) +
                           gamma * pMesh->getVertexNormals().col(i2)).normalized();
        } else {
            eRec.normal = (p1 - p0).cross(p2 - p0).normalized();
        }

        Vector3f wo = surfacePoint - eRec.point;
        float squaredDis = wo.squaredNorm();
        wo = wo.normalized();
        return eval(eRec) * std::max(0.f, eRec.normal.dot(wo)) / squaredDis / pdf(eRec);
    }

    float pdf(const EmitterQueryRecord &bRec) const override {
        return m_DPDF->getNormalization();
    }

    Color3f emission() const override {
        return m_radiance;
    }

    std::string toString() const override {
        return "AreaLight[]";
    }

private:
    NoriObject *m_parent{nullptr};
    Color3f m_radiance;
    std::shared_ptr<DiscretePDF> m_DPDF = nullptr;/// for sampling point from mesh
};

NORI_REGISTER_CLASS(AreaLight, "area");
NORI_NAMESPACE_END
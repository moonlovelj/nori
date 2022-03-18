/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
*/

#include <nori/scene.h>
#include <nori/bitmap.h>
#include <nori/integrator.h>
#include <nori/sampler.h>
#include <nori/camera.h>
#include <nori/emitter.h>

NORI_NAMESPACE_BEGIN

Scene::Scene(const PropertyList &) {
    m_accel = new Accel();
}

Scene::~Scene() {
    delete m_accel;
    delete m_sampler;
    delete m_camera;
    delete m_integrator;
}

void Scene::activate() {
    m_accel->build();

    if (!m_integrator)
        throw NoriException("No integrator was specified!");
    if (!m_camera)
        throw NoriException("No camera was specified!");
    
    if (!m_sampler) {
        /* Create a default (independent) sampler */
        m_sampler = static_cast<Sampler*>(
            NoriObjectFactory::createInstance("independent", PropertyList()));
    }

    cout << endl;
    cout << "Configuration: " << toString() << endl;
    cout << endl;
}

void Scene::addChild(NoriObject *obj) {
    switch (obj->getClassType()) {
        case EMesh: {
                Mesh *mesh = static_cast<Mesh *>(obj);
                m_accel->addMesh(mesh);
                m_meshes.push_back(mesh);
                if (mesh->isEmitter()) {
                  m_emitters.push_back(mesh);
                }
            }
            break;
        
        case EEmitter: {
                //Emitter *emitter = static_cast<Emitter *>(obj);
                /* TBD */
                throw NoriException("Scene::addChild(): You need to implement this for emitters");
            }
            break;

        case ESampler:
            if (m_sampler)
                throw NoriException("There can only be one sampler per scene!");
            m_sampler = static_cast<Sampler *>(obj);
            break;

        case ECamera:
            if (m_camera)
                throw NoriException("There can only be one camera per scene!");
            m_camera = static_cast<Camera *>(obj);
            break;
        case EMedium:
            if (m_medium)
                throw NoriException("There can only be one medium per scene!");
            m_medium = static_cast<Medium *>(obj);
            break;
        
        case EIntegrator:
            if (m_integrator)
                throw NoriException("There can only be one integrator per scene!");
            m_integrator = static_cast<Integrator *>(obj);
            break;

        default:
            throw NoriException("Scene::addChild(<%s>) is not supported!",
                classTypeName(obj->getClassType()));
    }
}

std::string Scene::toString() const {
    std::string meshes;
    for (size_t i=0; i<m_meshes.size(); ++i) {
        meshes += std::string("  ") + indent(m_meshes[i]->toString(), 2);
        if (i + 1 < m_meshes.size())
            meshes += ",";
        meshes += "\n";
    }

    return tfm::format(
        "Scene[\n"
        "  integrator = %s,\n"
        "  sampler = %s\n"
        "  camera = %s,\n"
        "  meshes = {\n"
        "  %s  }\n"
        "]",
        indent(m_integrator->toString()),
        indent(m_sampler->toString()),
        indent(m_camera->toString()),
        indent(meshes, 2)
    );
}

bool Scene::illuminatedEachOther(const Point3f &p0, const Point3f &p1) const {
    Vector3f dir = p1 - p0;
    Ray3f ray(p0, dir.normalized());
    ray.mint = Epsilon;
    ray.maxt = dir.norm() - Epsilon;
    return !this->rayIntersect(ray);
}

bool Scene::rayIntersectTr(const Ray3f &ray, Sampler *sampler,
                           Intersection &its, Color3f &tr) const {
    tr = Color3f(1.f);
    Ray3f r = ray;
    while (true) {
        bool hit = rayIntersect(r, its);
        if (m_medium)
            tr *= m_medium->tr(r, sampler);
        if (!hit)
            return false;
        if (its.mesh->getBSDF() != nullptr)
            return true;
        r = Ray3f(its.p + Epsilon*r.d, r.d);
    }
}

NORI_REGISTER_CLASS(Scene, "scene");
NORI_NAMESPACE_END

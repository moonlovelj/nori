/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
*/
#include <chrono>

#include <nori/accel.h>
#include <Eigen/Geometry>
#include <algorithm>
#include <cmath>

NORI_NAMESPACE_BEGIN

void Accel::addMesh(Mesh *mesh) {
  	meshes_.push_back(mesh);
    m_bbox.expandBy(mesh->getBoundingBox());
}

void Accel::build() {
  auto start_time = std::chrono::system_clock::now();
  accel_struct_ =
	  std::make_shared<Octree>(this->meshes_);

  auto end_time = std::chrono::system_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
	  end_time - start_time)
	  .count();
  std::cout << "\nAccel::build cost total time " << duration << "ms\n";
  std::cout << accel_struct_->ToString() << std::endl;
}

bool Accel::rayIntersect(const Ray3f &ray_, Intersection &its, bool shadowRay) const {
    bool foundIntersection = false;  // Was an intersection found so far?
    uint32_t f = (uint32_t) -1;      // Triangle index of the closest intersection
    Ray3f ray(ray_); /// Make a copy of the ray (we will need to update its '.maxt' value)

    /* Brute force search through all triangles */
//    for (auto mesh : meshes_) {
//	  for (uint32_t idx = 0; idx < mesh->getTriangleCount(); ++idx) {
//		float u, v, t;
//		if (mesh->rayIntersect(idx, ray, u, v, t)) {
//		  /* An intersection was found! Can terminate
//			 immediately if this is a shadow ray query */
//		  if (shadowRay)
//			return true;
//		  ray.maxt = its.t = t;
//		  its.uv = Point2f(u, v);
//		  its.mesh = mesh;
//		  f = idx;
//		  foundIntersection = true;
//		}
//	  }
//    }


    // Octree accel.
    foundIntersection = accel_struct_->RayIntersect(ray, its, shadowRay, f);
    if (shadowRay)
    	return foundIntersection;

    if (foundIntersection) {
        /* At this point, we now know that there is an intersection,
           and we know the triangle index of the closest such intersection.

           The following computes a number of additional properties which
           characterize the intersection (normals, texture coordinates, etc..)
        */

        /* Find the barycentric coordinates */
        Vector3f bary;
        bary << 1-its.uv.sum(), its.uv;

        /* References to all relevant mesh buffers */
        const Mesh *mesh   = its.mesh;
        const MatrixXf &V  = mesh->getVertexPositions();
        const MatrixXf &N  = mesh->getVertexNormals();
        const MatrixXf &UV = mesh->getVertexTexCoords();
        const MatrixXu &F  = mesh->getIndices();

        /* Vertex indices of the triangle */
        uint32_t idx0 = F(0, f), idx1 = F(1, f), idx2 = F(2, f);

        Point3f p0 = V.col(idx0), p1 = V.col(idx1), p2 = V.col(idx2);

        /* Compute the intersection positon accurately
           using barycentric coordinates */
        its.p = bary.x() * p0 + bary.y() * p1 + bary.z() * p2;

        /* Compute proper texture coordinates if provided by the mesh */
        if (UV.size() > 0)
            its.uv = bary.x() * UV.col(idx0) +
                bary.y() * UV.col(idx1) +
                bary.z() * UV.col(idx2);

        /* Compute the geometry frame */
        its.geoFrame = Frame((p1-p0).cross(p2-p0).normalized());

        if (N.size() > 0) {
            /* Compute the shading frame. Note that for simplicity,
               the current implementation doesn't attempt to provide
               tangents that are continuous across the surface. That
               means that this code will need to be modified to be able
               use anisotropic BRDFs, which need tangent continuity */

            its.shFrame = Frame(
                (bary.x() * N.col(idx0) +
                 bary.y() * N.col(idx1) +
                 bary.z() * N.col(idx2)).normalized());
        } else {
            its.shFrame = its.geoFrame;
        }
    }

    return foundIntersection;
}

NORI_NAMESPACE_END


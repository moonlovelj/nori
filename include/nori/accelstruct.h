//
// Created by 郭彬 on 2021/11/23.
//

#pragma once

#include <vector>
#include <nori/mesh.h>

NORI_NAMESPACE_BEGIN

class AccelStruct {
 public:
  AccelStruct(const std::vector<Mesh *> &meshes);
  virtual ~AccelStruct() {}
  virtual bool RayIntersect(Ray3f &ray,
							Intersection &its,
							bool shadowRay,
							uint32_t &face) const = 0;
  virtual std::string ToString() const { return ""; }
 protected:
  virtual void Build() = 0;
  std::vector<Mesh *> meshes_;
};

NORI_NAMESPACE_END

#pragma once

#include <iostream>
#include <memory>
#include <queue>
#include <vector>
#include <tuple>

#include <nori/bbox.h>
#include <nori/accelstruct.h>

NORI_NAMESPACE_BEGIN

const uint32_t kOctreeNodePrimitivesLimit = 20;
const uint32_t kOctreeMaxDepth = 20;

class Octree : public AccelStruct {
public:
  Octree(const std::vector<Mesh *> &meshes);

  Octree(const std::vector<Mesh *> &meshes,
		 const BoundingBox3f &bbox,
		 const std::vector<uint64_t> &face_indices,
		 uint32_t depth);

  bool RayIntersect(Ray3f &ray,
					Intersection &its,
					bool shadowRay,
					uint32_t &face) const override;

  std::string ToString() const override;

protected:
  void Build() override;

private:
  std::tuple<uint32_t, uint32_t> ParseFaceIndex(uint64_t value) const {
	return std::make_tuple(value >> 32, value & 0x00000000FFFFFFFF);
  }

  uint64_t EncodeFaceIndex(uint32_t mesh_index, uint32_t face_index) const {
	return static_cast<uint64_t>(mesh_index) << 32 | face_index;
  }

  void Build(const BoundingBox3f &bbox, const std::vector<uint64_t> &facesIndices, uint32_t depth);

  std::vector<std::shared_ptr<Octree>> children_;
  std::vector<uint64_t> facesIndices_; /// mesh index is stored in the first 32 bits and face index in the last 32 bits.
  BoundingBox3f bbox_;
};

NORI_NAMESPACE_END
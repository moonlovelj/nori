#pragma once

#include <iostream>
#include <memory>
#include <nori/bbox.h>
#include <queue>
#include <vector>

NORI_NAMESPACE_BEGIN

const uint32_t kOctreeNodePrimitivesLimit = 20;
const uint32_t kOctreeMaxDepth = 20;

class Octree {
public:
  Octree(const BoundingBox3f &bbox, const MatrixXf &vertices,
         const MatrixXu &faces, const std::vector<uint32_t> &facesIndices,
         uint32_t depth = 0);

  void OutputDebugInfo();

  std::vector<std::shared_ptr<Octree>> children_;
  std::vector<uint32_t> facesIndices_;
  BoundingBox3f bbox_;

private:
  void Build(const BoundingBox3f &bbox, const MatrixXf &vertices,
             const MatrixXu &faces, const std::vector<uint32_t> &facesIndices,
             uint32_t depth);
};

NORI_NAMESPACE_END
#include <nori/octree.h>
#include <queue>

NORI_NAMESPACE_BEGIN

Octree::Octree(const BoundingBox3f &bbox, const MatrixXf &vertices,
               const MatrixXu &faces, const std::vector<uint32_t> &facesIndices,
               uint32_t depth) {
  Build(bbox, vertices, faces, facesIndices, depth);
}

void Octree::Build(const BoundingBox3f &bbox, const MatrixXf &vertices,
                   const MatrixXu &faces,
                   const std::vector<uint32_t> &facesIndices, uint32_t depth) {
  if (facesIndices.size() == 0)
    return;

  bbox_ = bbox;
  if (facesIndices.size() < kOctreeNodePrimitivesLimit ||
      depth > kOctreeMaxDepth) {
    facesIndices_ = facesIndices;
    return;
  }

  auto bbox_center = bbox.getCenter();
  std::vector<std::vector<uint32_t>> face_list(8);
  std::vector<BoundingBox3f> bboxs(8);
  bboxs[0] = BoundingBox3f(bbox.min, bbox_center);
  bboxs[1] =
      BoundingBox3f(Point3f(bbox.min.x(), bbox.min.y(), bbox_center.z()),
                    Point3f(bbox_center.x(), bbox_center.y(), bbox.max.z()));
  bboxs[2] =
      BoundingBox3f(Point3f(bbox_center.x(), bbox.min.y(), bbox_center.z()),
                    Point3f(bbox.max.x(), bbox_center.y(), bbox.max.z()));
  bboxs[3] =
      BoundingBox3f(Point3f(bbox_center.x(), bbox.min.y(), bbox.min.z()),
                    Point3f(bbox.max.x(), bbox_center.y(), bbox_center.z()));
  bboxs[4] =
      BoundingBox3f(Point3f(bbox.min.x(), bbox_center.y(), bbox.min.z()),
                    Point3f(bbox_center.x(), bbox.max.y(), bbox_center.z()));
  bboxs[5] =
      BoundingBox3f(Point3f(bbox.min.x(), bbox_center.y(), bbox_center.z()),
                    Point3f(bbox_center.x(), bbox.max.y(), bbox.max.z()));
  bboxs[6] =
      BoundingBox3f(Point3f(bbox_center.x(), bbox_center.y(), bbox_center.z()),
                    Point3f(bbox.max.x(), bbox.max.y(), bbox.max.z()));
  bboxs[7] =
      BoundingBox3f(Point3f(bbox_center.x(), bbox_center.y(), bbox.min.z()),
                    Point3f(bbox.max.x(), bbox.max.y(), bbox_center.z()));

  for (size_t i = 0; i < facesIndices.size(); i++) {
    for (size_t j = 0; j < bboxs.size(); j++) {
      // for (size_t index = 0; index < 3; index++) {
      //   if (bboxs[j].contains(vertices.col(faces(index, facesIndices[i])))) {
      //     face_list[j].emplace_back(facesIndices[i]);
      //     break;
      //   }
      // }
      BoundingBox3f tb(vertices.col(faces(0, facesIndices[i])));
      tb.expandBy(vertices.col(faces(1, facesIndices[i])));
      tb.expandBy(vertices.col(faces(2, facesIndices[i])));
      if (bboxs[j].overlaps(tb)) {
        face_list[j].emplace_back(facesIndices[i]);
      }
    }
  }

  children_.resize(8);
  for (size_t i = 0; i < 8; i++) {
    children_[i] = std::make_shared<Octree>(bboxs[i], vertices, faces,
                                            face_list[i], depth + 1);
  }
}

void Octree::OutputDebugInfo() {
  int interior_node_num = 0;
  int leaf_node_num = 0;
  int total_triangle_num = 0;
  std::cout << std::endl;
  std::queue<const Octree *> record;
  record.push(this);
  while (!record.empty()) {
    auto node = record.front();
    record.pop();
    if (node->children_.size() == 0) {
      ++leaf_node_num;
      total_triangle_num += node->facesIndices_.size();
    } else {
      ++interior_node_num;
      for (auto child : node->children_) {
        record.push(child.get());
      }
    }
  }

  std::cout << "Interior node num is : " << interior_node_num << std::endl;
  std::cout << "Leaf node num is : " << leaf_node_num << std::endl;
  // std::cout << "total_triangle_num : " << total_triangle_num << std::endl;
  // std::cout << "Average number of triangles per leaf node is : "
  //           << this->facesIndices_.size() / (double)leaf_node_num <<
  //           std::endl;
  std::cout << std::endl;
}

NORI_NAMESPACE_END
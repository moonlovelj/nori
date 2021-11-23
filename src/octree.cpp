#include <nori/octree.h>
#include <queue>

NORI_NAMESPACE_BEGIN

Octree::Octree(const std::vector<Mesh *> &meshes) : AccelStruct(meshes) {
  Build();
}

void Octree::Build() {
  BoundingBox3f bbox;
  std::vector<uint64_t> facesIndices;
  for (uint32_t i = 0; i < meshes_.size(); ++i) {
	bbox.expandBy(meshes_[i]->getBoundingBox());
	auto face_count = meshes_[i]->getTriangleCount();
	for (uint32_t face_index = 0; face_index < face_count; ++face_index) {
	  facesIndices.push_back(EncodeFaceIndex(i, face_index));
	}
  }
  Build(bbox, facesIndices, 0);
}

Octree::Octree(const std::vector<Mesh *> &meshes,
			   const BoundingBox3f &bbox,
			   const std::vector<uint64_t> &face_indices,
			   uint32_t depth) : AccelStruct(meshes) {
  Build(bbox, face_indices, depth);
}

void Octree::Build(const BoundingBox3f &bbox, const std::vector<uint64_t> &facesIndices, uint32_t depth) {
  if (facesIndices.size() == 0)
	return;

  bbox_ = bbox;
  if (facesIndices.size() < kOctreeNodePrimitivesLimit ||
	  depth > kOctreeMaxDepth) {
	facesIndices_ = facesIndices;
	return;
  }

  auto bbox_center = bbox.getCenter();
  std::vector<std::vector<uint64_t>> face_list(8);
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
	auto[mesh_index, face_index] = ParseFaceIndex(facesIndices[i]);
	auto box_face = this->meshes_[mesh_index]->getBoundingBox(face_index);
	for (size_t j = 0; j < bboxs.size(); j++) {
	  if (bboxs[j].overlaps(box_face)) {
		face_list[j].emplace_back(EncodeFaceIndex(mesh_index, face_index));
	  }
	}
  }

  children_.resize(8);
  for (size_t i = 0; i < 8; i++) {
	children_[i] = std::make_shared<Octree>(this->meshes_, bboxs[i], face_list[i], depth + 1);
  }
}

bool Octree::RayIntersect(Ray3f &ray,
						  Intersection &its,
						  bool shadowRay,
						  uint32_t &face) const {
  if (this->facesIndices_.size() > 0) {
	// leaf node
	bool foundIntersection = false;
	for (const auto &faces : this->facesIndices_) {
	  auto[mesh_index, face_index] = ParseFaceIndex(faces);
	  float u, v, t;
	  if (meshes_[mesh_index]->rayIntersect(face_index, ray, u, v, t)) {
		/* An intersection was found! Can terminate
		   immediately if this is a shadow ray query */
		if (shadowRay)
		  return true;
		ray.maxt = its.t = t;
		its.uv = Point2f(u, v);
		its.mesh = meshes_[mesh_index];
		face = face_index;
		foundIntersection = true;
	  }
	}
	return foundIntersection;
  }

  bool intersectChild = false;
  // Base version
  for (auto child : this->children_) {
	if (child->bbox_.rayIntersect(ray) && child->RayIntersect(ray, its, shadowRay, face)) {
	  if (shadowRay)
		return true;
	  intersectChild = true;
	}
  }

  //TODO this have problem of multiple thread.

  // std::vector<std::shared_ptr<Octree>> clone_children(octree->children_);
  // std::sort(clone_children.begin(), clone_children.end(),
  //           [&ray](auto a, auto b) {
  //             float ta, tb;
  //             bool intersect_a = a->bbox_.rayIntersect(ray, ta);
  //             bool intersect_b = b->bbox_.rayIntersect(ray, tb);
  //             return std::isless(ta, tb);
  //           });
  // for (auto child : clone_children) {
  //   if (!child->bbox_.rayIntersect(ray))
  //     break;

  //   if (rayIntersectOctree(child, ray, its, shadowRay, f)) {
  //     intersectChild = true;
  //     break;
  //   }
  // }

  return intersectChild;
}

std::string Octree::ToString() const {
  std::string str;
  int interior_node_num = 0;
  int leaf_node_num = 0;
  int total_triangle_num = 0;
  std::queue<const Octree *> record;
  record.push(this);
  while (!record.empty()) {
	auto node = record.front();
	record.pop();
	if (node->facesIndices_.size() > 0) {
	  ++leaf_node_num;
	  total_triangle_num += node->facesIndices_.size();
	} else {
	  ++interior_node_num;
	  for (auto child : node->children_) {
		record.push(child.get());
	  }
	}
  }

  str += "Name : Octree\n";

  str += "Interior node num is : ";
  str += std::to_string(interior_node_num);
  str += "\n";

  str += "Leaf node num is : ";
  str += std::to_string(leaf_node_num);
  str += "\n";

  str += "Total_triangle_num is : ";
  str += std::to_string(total_triangle_num);
  str += "\n";

  str += "Average number of triangles per leaf node is : ";
  str += std::to_string(total_triangle_num / (double)leaf_node_num);
  str += "\n";

  return str;
}

NORI_NAMESPACE_END
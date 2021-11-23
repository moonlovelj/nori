//
// Created by 郭彬 on 2021/11/22.
//


#include <nori/emitter.h>

NORI_NAMESPACE_BEGIN

class AreaLight : public Emitter {
 public:
  AreaLight(const PropertyList &props) {
	  radiance_ = props.getColor("radiance", Color3f(1.f));
  }

  Color3f Emission() const override {
  	  return radiance_;
  }

  std::string toString() const override {
	  return "AreaLight[]";
  }
 private:
  	Color3f radiance_;
};

NORI_REGISTER_CLASS(AreaLight, "area");
NORI_NAMESPACE_END
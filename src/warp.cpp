/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
*/

#include <nori/warp.h>
#include <nori/vector.h>
#include <nori/frame.h>

NORI_NAMESPACE_BEGIN

static double Pi = 3.14159265358979;

Point2f Warp::squareToUniformSquare(const Point2f &sample) {
    return sample;
}

float Warp::squareToUniformSquarePdf(const Point2f &sample) {
    return ((sample.array() >= 0).all() && (sample.array() <= 1).all()) ? 1.0f : 0.0f;
}

Point2f Warp::squareToTent(const Point2f &sample) {
	Point2f p;
	p(0) = sample.x() >= 0.5f ? 1-std::sqrt(2*(1.f-sample.x())) : -1 + std::sqrt(2*sample.x());
	p(1) = sample.y() >= 0.5f ? 1-std::sqrt(2*(1.f-sample.y())) : -1 + std::sqrt(2*sample.y());
	return p;
}

float Warp::squareToTentPdf(const Point2f &p) {
    return ((p.array() >= -1).all() && (p.array() <= 1).all()) ? (1.0f-std::fabs(p.x()))*(1.0f-std::fabs(p.y())) : 0.0f;
}

Point2f Warp::squareToUniformDisk(const Point2f &sample) {
	auto r = std::sqrt(sample.x());
	auto theta = sample.y() * 2 * Pi;
	return Point2f(r * std::cosf(theta), r * std::sinf(theta));
}

float Warp::squareToUniformDiskPdf(const Point2f &p) {
    return p.x() * p.x() + p.y() * p.y() <= 1.f ? 1.f / Pi : 0;
}

Vector3f Warp::squareToUniformSphere(const Point2f &sample) {
    auto theta = std::acosf(1.f-2*sample.x());
    auto phi = 2 * Pi * sample.y();
    return Vector3f(std::sinf(theta)*std::cosf(phi),std::cosf(theta),std::sinf(theta)*std::sinf(phi));
}

float Warp::squareToUniformSpherePdf(const Vector3f &v) {
	return 1.f / (4 * Pi);
}

Vector3f Warp::squareToUniformHemisphere(const Point2f &sample) {
	auto theta = std::acosf(1.f-2*sample.x());
	auto phi = Pi * sample.y();
	return Vector3f(std::sinf(theta)*std::cosf(phi),std::cosf(theta),std::sinf(theta)*std::sinf(phi));
}

float Warp::squareToUniformHemispherePdf(const Vector3f &v) {
	return v.z() >= 0 ? 1.f / (2 * Pi) : 0;
}

Vector3f Warp::squareToCosineHemisphere(const Point2f &sample) {
	auto r = std::sqrt(sample.x());
	auto theta = sample.y() * 2 * Pi;
	auto x = r * std::cosf(theta);
	auto y = r * std::sinf(theta);
	return Vector3f(x, y, std::sqrtf(1.f - x*x - y*y));
}

float Warp::squareToCosineHemispherePdf(const Vector3f &v) {
    auto cos_theta = Vector3f(0,0,1).dot(v.normalized());
    return cos_theta >= 0.f ? cos_theta/Pi : 0;
}

Vector3f Warp::squareToBeckmann(const Point2f &sample, float alpha) {
    throw NoriException("Warp::squareToBeckmann() is not yet implemented!");
}

float Warp::squareToBeckmannPdf(const Vector3f &m, float alpha) {
    if (m.z() < 0) return 0;

    auto theta = std::acosf(m.z());
	return 1 / (2 * Pi) * 2 * std::expf(-std::tanf(theta)*std::tanf(theta)/(alpha*alpha)) / (alpha*alpha*m.z()*m.z()*m.z());
}

NORI_NAMESPACE_END

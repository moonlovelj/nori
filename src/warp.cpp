/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
*/

#include <nori/warp.h>
#include <nori/vector.h>
#include <nori/frame.h>

NORI_NAMESPACE_BEGIN

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
	auto theta = sample.y() * 2 * M_PI;
	return Point2f(r * std::cosf(theta), r * std::sinf(theta));
}

float Warp::squareToUniformDiskPdf(const Point2f &p) {
    return p.x() * p.x() + p.y() * p.y() <= 1.f ? M_1_PI : 0;
}

Vector3f Warp::squareToUniformSphere(const Point2f &sample) {
    auto theta = std::acosf(1.f-2*sample.x());
    auto phi = 2 * M_PI * sample.y();
    return Vector3f(std::sinf(theta)*std::cosf(phi),std::cosf(theta),std::sinf(theta)*std::sinf(phi));
}

float Warp::squareToUniformSpherePdf(const Vector3f &v) {
	return 1.f / (4 * M_PI);
}

Vector3f Warp::squareToUniformHemisphere(const Point2f &sample) {
	auto theta = std::acosf(1.f-2*sample.x());
	auto phi = M_PI * sample.y();
	return Vector3f(std::sinf(theta)*std::cosf(phi),std::cosf(theta),std::sinf(theta)*std::sinf(phi));
}

float Warp::squareToUniformHemispherePdf(const Vector3f &v) {
	return v.z() >= 0 ? 1.f / (2 * M_PI) : 0;
}

Vector3f Warp::squareToCosineHemisphere(const Point2f &sample) {
	auto r = std::sqrt(sample.x());
	auto theta = sample.y() * 2 * M_PI;
	auto x = r * std::cosf(theta);
	auto y = r * std::sinf(theta);
	return Vector3f(x, y, std::sqrtf(1.f - x*x - y*y));
}

float Warp::squareToCosineHemispherePdf(const Vector3f &v) {
    auto cos_theta = Vector3f(0,0,1).dot(v.normalized());
    return cos_theta >= 0.f ? cos_theta*M_1_PI : 0;
}

Vector3f Warp::squareToBeckmann(const Point2f &sample, float alpha) {
	auto ln = logf(1.f - sample.x());
	auto tan_theta2 = -alpha*alpha*ln;
	auto cos_theta2 = 1/(1+tan_theta2);
    auto sin_theta = std::sqrtf(std::max(0.f,1-cos_theta2));
    auto phi = 2 * M_PI * sample.y();
    return Vector3f(sin_theta*std::cosf(phi), sin_theta*std::sinf(phi), std::sqrtf(cos_theta2));
}

float Warp::squareToBeckmannPdf(const Vector3f &m, float alpha) {
    if (m.z() <= 0) return 0;
    auto cos_theta = m.z();
    auto cos_theta2 = cos_theta * cos_theta;
    auto tan_theta2 = (1-cos_theta2)/cos_theta2;
    auto alpha2 = alpha * alpha;
    float cos_theta3 = cos_theta2*cos_theta;
    float pdf = std::expf(-tan_theta2/alpha2) / (M_PI * alpha2*cos_theta3);
    return pdf > 1e-20 ? pdf : 0;
}

Vector3f Warp::squareToGGX(const Point2f &sample, float alpha) {
    float tan_theta = alpha * sqrt(sample.x()) / sqrt(1.f - sample.x());
    float theta = std::atan(tan_theta);
    float cos_theta = std::cos(theta);
    float sin_theta = std::sqrtf(std::max(0.f, 1 - cos_theta * cos_theta));
    float phi = 2 * M_PI * sample.y();
    return Vector3f(sin_theta * std::cosf(phi), sin_theta * std::sinf(phi), cos_theta);
}

float Warp::squareToGGXPdf(const Vector3f &m, float alpha) {
    if (m.z() <= 0) return 0;
    float cos_theta = m.z();
    float cos_theta2 = cos_theta * cos_theta;
    float tan_theta2 = std::max(0.f, (1 - cos_theta2) / cos_theta2);
    float alpha2 = alpha * alpha;
    float denominator = M_PI * cos_theta * cos_theta2 * (alpha2 + tan_theta2) * (alpha2 + tan_theta2);
    float pdf = alpha2 / denominator;
    return pdf > 1e-20 ? pdf : 0;
}

float Warp::phaseHG(float cosTheta, float g) {
    float denom = 1 + g * g + 2 * g * cosTheta;
    return M_1_PI * 0.25f * (1 - g * g) / (denom * std::sqrt(denom));
}

NORI_NAMESPACE_END

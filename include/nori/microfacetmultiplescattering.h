//
// Created by 郭彬 on 2022/4/20.
//
// https://eheitzresearch.wordpress.com/240-2/

#ifndef NORI_MICROFACETMULTIPLESCATTERING_H
#define NORI_MICROFACETMULTIPLESCATTERING_H

#include <memory>

#include <nori/bsdf.h>

NORI_NAMESPACE_BEGIN

/************* MICROSURFACE HEIGHT DISTRIBUTION *************/

/* API */
class MicrosurfaceHeight
{
public:
    virtual ~MicrosurfaceHeight(){}
    // height PDF
    virtual float P1(const float h) const=0;
    // height CDF
    virtual float C1(const float h) const=0;
    // inverse of the height CDF
    virtual float invC1(const float U) const=0;
};

/* Uniform height distribution in [-1, 1] */
class MicrosurfaceHeightUniform : public MicrosurfaceHeight
{
public:
    // height PDF
    virtual float P1(const float h) const;
    // height CDF
    virtual float C1(const float h) const;
    // inverse of the height CDF
    virtual float invC1(const float U) const;
};

/* Gaussian height distribution N(0,1) */
class MicrosurfaceHeightGaussian : public MicrosurfaceHeight
{
public:
    // height PDF
    virtual float P1(const float h) const;
    // height CDF
    virtual float C1(const float h) const;
    // inverse of the height CDF
    virtual float invC1(const float U) const;
};

/************* MICROSURFACE SLOPE DISTRIBUTION *************/

/* API */
class MicrosurfaceSlope
{
public:
    MicrosurfaceSlope(const float alpha_x=1.0f, const float alpha_y=1.0f)
            : m_alpha_x(alpha_x), m_alpha_y(alpha_y)
    {}
    virtual ~MicrosurfaceSlope(){}
public:
    // roughness
    const float m_alpha_x, m_alpha_y;
    // projected roughness in wi
    float alpha_i(const Vector3f& wi) const;

public:
    // distribution of normals (NDF)
    float D(const Vector3f& wm) const;
    // distribution of visible normals (VNDF)
    float D_wi(const Vector3f& wi, const Vector3f& wm) const;
    // sample the VNDF
    Vector3f sampleD_wi(const Vector3f& wi, const float U1, const float U2) const;

public:
    // distribution of slopes
    virtual float P22(const float slope_x, const float slope_y) const=0;
    // Smith's Lambda function
    virtual float Lambda(const Vector3f& wi) const=0;
    // projected area towards incident direction
    virtual float projectedArea(const Vector3f& wi) const=0;
    // sample the distribution of visible slopes with alpha=1.0
    virtual Vector2f sampleP22_11(const float theta_i, const float U1, const float U2) const=0;
};

/* Beckmann slope distribution */
class MicrosurfaceSlopeBeckmann : public MicrosurfaceSlope
{
public:
    MicrosurfaceSlopeBeckmann(const float alpha_x=1.0f, const float alpha_y=1.0f)
            : MicrosurfaceSlope(alpha_x, alpha_y)
    {}

    // distribution of slopes
    virtual float P22(const float slope_x, const float slope_y) const;
    // Smith's Lambda function
    virtual float Lambda(const Vector3f& wi) const;
    // projected area towards incident direction
    virtual float projectedArea(const Vector3f& wi) const;
    // sample the distribution of visible slopes with alpha=1.0
    virtual Vector2f sampleP22_11(const float theta_i, const float U1, const float U2) const;
};

/* GGX slope distribution */
class MicrosurfaceSlopeGGX : public MicrosurfaceSlope
{
public:
    MicrosurfaceSlopeGGX(const float alpha_x=1.0f, const float alpha_y=1.0f)
            : MicrosurfaceSlope(alpha_x, alpha_y)
    {}

    // distribution of slopes
    virtual float P22(const float slope_x, const float slope_y) const;
    // Smith's Lambda function
    virtual float Lambda(const Vector3f& wi) const;
    // projected area towards incident direction
    virtual float projectedArea(const Vector3f& wi) const;
    // sample the distribution of visible slopes with alpha=1.0
    virtual Vector2f sampleP22_11(const float theta_i, const float U1, const float U2) const;
};


/************* MICROSURFACE *************/

/* API */
class Microsurface : public BSDF
{
public:
    Microsurface(const PropertyList &);

    virtual ~Microsurface() {}

    // evaluate BSDF with a random walk (stochastic but unbiased)
    // scatteringOrder=0 --> contribution from all scattering events
    // scatteringOrder=1 --> contribution from 1st bounce only
    // scatteringOrder=2 --> contribution from 2nd bounce only, etc..
    virtual float eval(const Vector3f& wi, const Vector3f& wo, Sampler *sampler, const int scatteringOrder=0) const;

    // sample BSDF with a random walk
    // scatteringOrder is set to the number of bounces computed for this sample
    virtual Vector3f sample(const Vector3f& wi, int& scatteringOrder, Sampler *sampler) const;
    Vector3f sample(const Vector3f& wi, Sampler *sampler) const {int scatteringOrder; return sample(wi, scatteringOrder, sampler);}

    virtual Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample) const = 0;
    virtual Color3f eval(const BSDFQueryRecord &bRec) const = 0;
    virtual float pdf(const BSDFQueryRecord &bRec) const = 0;
public:
    // masking function
    float G_1(const Vector3f& wi) const;
    // masking function at height h0
    float G_1(const Vector3f& wi, const float h0) const;
    // sample height in outgoing direction
    float sampleHeight(const Vector3f& wo, const float h0, const float U) const;

public:
    // evaluate local phase function 
    virtual float evalPhaseFunction(const Vector3f& wi, const Vector3f& wo) const=0;
    // sample local phase function
    virtual Vector3f samplePhaseFunction(const Vector3f& wi) const=0;

    // evaluate BSDF limited to single scattering 
    // this is in average equivalent to eval(wi, wo, 1);
    virtual float evalSingleScattering(const Vector3f& wi, const Vector3f& wo) const=0;
public:
    // height distribution
    std::unique_ptr<const MicrosurfaceHeight> m_microsurfaceheight;
    // slope distribution
    std::unique_ptr<const MicrosurfaceSlope> m_microsurfaceslope;
};

NORI_NAMESPACE_END

#endif //NORI_MICROFACETMULTIPLESCATTERING_H

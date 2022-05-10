//
// Created by 郭彬 on 2022/4/20.
//

#include <nori/microfacetmultiplescattering.h>
#include <nori/sampler.h>

NORI_NAMESPACE_BEGIN

#define SQRT_2			M_SQRT2 /* sqrt(2) */
#define INV_M_PI		0.31830988618379067153f /* 1/pi */
#define SQRT_M_PI		1.77245385090551602729f /* sqrt(pi) */
#define INV_SQRT_M_PI	0.56418958354775628694f /* 1/sqrt(pi) */
#define INV_2_SQRT_M_PI	0.28209479177387814347f /* 0.5/sqrt(pi) */
#define INV_SQRT_2_M_PI 0.3989422804014326779f /* 1/sqrt(2*pi) */
#define INV_SQRT_2		0.7071067811865475244f /* 1/sqrt(2) */

static float sign(float a)
{
    return a > 0 ? 1 : -1;
}

static bool IsFiniteNumber(float x)
{
    return (x <= std::numeric_limits<float>::max() && x >= -std::numeric_limits<float>::max());
}

static double erf(double x)
{
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;

    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x);

    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

    return sign*y;
}

static float erfinv(float x)
{
    float w, p;
    w = - logf((1.0f-x)*(1.0f+x));
    if ( w < 5.000000f ) {
        w = w - 2.500000f;
        p = 2.81022636e-08f;
        p = 3.43273939e-07f + p*w;
        p = -3.5233877e-06f + p*w;
        p = -4.39150654e-06f + p*w;
        p = 0.00021858087f + p*w;
        p = -0.00125372503f + p*w;
        p = -0.00417768164f + p*w;
        p = 0.246640727f + p*w;
        p = 1.50140941f + p*w;
    }
    else {
        w = sqrtf(w) - 3.000000f;
        p = -0.000200214257f;
        p = 0.000100950558f + p*w;
        p = 0.00134934322f + p*w;
        p = -0.00367342844f + p*w;
        p = 0.00573950773f + p*w;
        p = -0.0076224613f + p*w;
        p = 0.00943887047f + p*w;
        p = 1.00167406f + p*w;
        p = 2.83297682f + p*w;
    }
    return p*x;
}

/*
 * A method to compute the gamma() function.
 *
 */

static double  abgam (double x)
{
    double  gam[10],
            temp;

    gam[0] = 1./ 12.;
    gam[1] = 1./ 30.;
    gam[2] = 53./ 210.;
    gam[3] = 195./ 371.;
    gam[4] = 22999./ 22737.;
    gam[5] = 29944523./ 19733142.;
    gam[6] = 109535241009./ 48264275462.;
    temp = 0.5*log (2*M_PI) - x + (x - 0.5)*log (x)
           + gam[0]/(x + gam[1]/(x + gam[2]/(x + gam[3]/(x + gam[4] /
                                                             (x + gam[5]/(x + gam[6]/x))))));

    return temp;
}

static double  gamma (double x)
{
    double  result;
    result = exp (abgam (x + 5))/(x*(x + 1)*(x + 2)*(x + 3)*(x + 4));
    return result;
}

static double  beta (double m, double n)
{
    return (gamma (m)*gamma (n)/gamma (m + n));
}


/************* MICROSURFACE HEIGHT DISTRIBUTION *************/

float MicrosurfaceHeightUniform::P1(const float h) const
{
    const float value = (h >= -1.0f && h <= 1.0f) ? 0.5f : 0.0f;
    return value;
}

float MicrosurfaceHeightUniform::C1(const float h) const
{
    const float value = std::min(1.0f, std::max(0.0f, 0.5f*(h+1.0f)));
    return value;
}

float MicrosurfaceHeightUniform::invC1(const float U) const
{
    const float h = std::max(-1.0f, std::min(1.0f, 2.0f*U-1.0f));
    return h;
}

float MicrosurfaceHeightGaussian::P1(const float h) const
{
    const float value = INV_SQRT_2_M_PI * expf(-0.5f * h*h);
    return value;
}

float MicrosurfaceHeightGaussian::C1(const float h) const
{
    const float value = 0.5f + 0.5f * (float)erf(INV_SQRT_2*h);
    return value;
}

float MicrosurfaceHeightGaussian::invC1(const float U) const
{
    const float h = SQRT_2 * erfinv(2.0f*U - 1.0f);
    return h;
}


/************* MICROSURFACE SLOPE DISTRIBUTION *************/

float MicrosurfaceSlope::D(const Vector3f& wm) const {
    if( wm.z() <= 0.0f)
        return 0.0f;

    // slope of wm
    const float slope_x = -wm.x()/wm.z();
    const float slope_y = -wm.y()/wm.z();

    // value
    const float value = P22(slope_x, slope_y) / (wm.z()*wm.z()*wm.z()*wm.z());
    return value;
}

float MicrosurfaceSlope::D_wi(const Vector3f& wi, const Vector3f& wm) const {
    if( wm.z() <= 0.0f)
        return 0.0f;

    // normalization coefficient
    const float projectedarea = projectedArea(wi);
    if(projectedarea == 0)
        return 0;
    const float c = 1.0f / projectedarea;

    // value
    const float value = c * std::max(0.0f, wi.dot(wm)) * D(wm);
    return value;
}

Vector3f MicrosurfaceSlope::sampleD_wi(const Vector3f& wi, const float U1, const float U2) const {

    // stretch to match configuration with alpha=1.0
    const Vector3f wi_11 = Vector3f(m_alpha_x * wi.x(), m_alpha_y * wi.y(), wi.z()).normalized();

    // sample visible slope with alpha=1.0
    Vector2f slope_11 = sampleP22_11(acosf(wi_11.z()), U1, U2);

    // align with view direction
    const float phi = atan2(wi_11.y(), wi_11.x());
    Vector2f slope(cosf(phi)*slope_11.x() - sinf(phi)*slope_11.y(), sinf(phi)*slope_11.x() + cos(phi)*slope_11.y());

    // stretch back
    slope.x() *= m_alpha_x;
    slope.y() *= m_alpha_y;

    // if numerical instability
    if( (slope.x() != slope.x()) || !IsFiniteNumber(slope.x()) )
    {
        if(wi.z() > 0) return Vector3f(0.0f,0.0f,1.0f);
        else return Vector3f(wi.x(), wi.y(), 0.0f).normalized();
    }

    // compute normal
    const Vector3f wm = Vector3f(-slope.x(), -slope.y(), 1.0f).normalized();
    return wm;
}

float MicrosurfaceSlope::alpha_i(const Vector3f& wi) const
{
    const float invSinTheta2 = 1.0f / (1.0f - wi.z()*wi.z());
    const float cosPhi2 = wi.x()*wi.x()*invSinTheta2;
    const float sinPhi2 = wi.y()*wi.y()*invSinTheta2;
    const float alpha_i = sqrtf( cosPhi2*m_alpha_x*m_alpha_x + sinPhi2*m_alpha_y*m_alpha_y );
    return alpha_i;
}

float MicrosurfaceSlopeBeckmann::P22(const float slope_x, const float slope_y) const
{
    const float value = 1.0f / (M_PI * m_alpha_x * m_alpha_y) * expf(-slope_x*slope_x/(m_alpha_x*m_alpha_x) - slope_y*slope_y/(m_alpha_y*m_alpha_y) );
    return value;
}

float MicrosurfaceSlopeBeckmann::Lambda(const Vector3f& wi) const
{
    if(wi.z() > 0.9999f)
        return 0.0f;
    if(wi.z() < -0.9999f)
        return -1.0f;

    // a
    const float theta_i = acosf(wi.z());
    const float a = 1.0f/tanf(theta_i)/alpha_i(wi);

    // value
    const float value = 0.5f*((float)erf(a) - 1.0f) + INV_2_SQRT_M_PI / a * expf(-a*a);

    return value;
}

float MicrosurfaceSlopeBeckmann::projectedArea(const Vector3f& wi) const
{
    if(wi.z() > 0.9999f)
        return 1.0f;
    if(wi.z() < -0.9999f)
        return 0.0f;

    // a
    const float alphai = alpha_i(wi);
    const float theta_i = acosf(wi.z());
    const float a = 1.0f/tanf(theta_i)/alphai;

    // value
    const float value = 0.5f*((float)erf(a) + 1.0f)*wi.z() + INV_2_SQRT_M_PI * alphai * sinf(theta_i) * expf(-a*a);

    return value;
}

Vector2f MicrosurfaceSlopeBeckmann::sampleP22_11(const float theta_i, const float U, const float U_2) const
{
    Vector2f slope;

    if(theta_i < 0.0001f)
    {
        const float r = sqrtf(-logf(U));
        const float phi = 6.28318530718f * U_2;
        slope.x() = r * cosf(phi);
        slope.y() = r * sinf(phi);
        return slope;
    }

    // constant
    const float sin_theta_i = sinf(theta_i);
    const float cos_theta_i = cosf(theta_i);

    // slope associated to theta_i
    const float slope_i = cos_theta_i/sin_theta_i;

    // projected area
    const float a = cos_theta_i/sin_theta_i;
    const float projectedarea = 0.5f*((float)erf(a) + 1.0f)*cos_theta_i + INV_2_SQRT_M_PI * sin_theta_i * expf(-a*a);
    if(projectedarea < 0.0001f || projectedarea!=projectedarea)
        return Vector2f(0,0);
    // VNDF normalization factor
    const float c = 1.0f / projectedarea;

    // search
    float erf_min = -0.9999f;
    float erf_max = std::max(erf_min, (float)erf(slope_i));
    float erf_current = 0.5f * (erf_min+erf_max);

    while(erf_max-erf_min > 0.00001f)
    {
        if (!(erf_current >= erf_min && erf_current <= erf_max))
            erf_current = 0.5f * (erf_min + erf_max);

        // evaluate slope
        const float slope = erfinv(erf_current);

        // CDF
        const float CDF = (slope>=slope_i) ? 1.0f : c * (INV_2_SQRT_M_PI*sin_theta_i*expf(-slope*slope) + cos_theta_i*(0.5f+0.5f*(float)erf(slope)));
        const float diff = CDF - U;

        // test estimate
        if( abs(diff) < 0.00001f )
            break;

        // update bounds
        if(diff > 0.0f)
        {
            if(erf_max == erf_current)
                break;
            erf_max = erf_current;
        }
        else
        {
            if(erf_min == erf_current)
                break;
            erf_min = erf_current;
        }

        // update estimate
        const float derivative = 0.5f*c*cos_theta_i - 0.5f*c*sin_theta_i * slope;
        erf_current -= diff/derivative;
    }

    slope.x() = erfinv(std::min(erf_max, std::max(erf_min, erf_current)));
    slope.y() = erfinv(2.0f*U_2-1.0f);
    return slope;
}

float MicrosurfaceSlopeGGX::P22(const float slope_x, const float slope_y) const
{
    const float tmp = 1.0f + slope_x*slope_x/(m_alpha_x*m_alpha_x) + slope_y*slope_y/(m_alpha_y*m_alpha_y);
    const float value = 1.0f / (M_PI * m_alpha_x * m_alpha_y) / (tmp * tmp);
    return value;
}

float MicrosurfaceSlopeGGX::Lambda(const Vector3f& wi) const
{
    if(wi.z() > 0.9999f)
        return 0.0f;
    if(wi.z() < -0.9999f)
        return -1.0f;

    // a
    const float theta_i = acosf(wi.z());
    const float a = 1.0f/tanf(theta_i)/alpha_i(wi);

    // value
    const float value = 0.5f*(-1.0f + sign(a) * sqrtf(1 + 1/(a*a)));

    return value;
}

float MicrosurfaceSlopeGGX::projectedArea(const Vector3f& wi) const
{
    if(wi.z() > 0.9999f)
        return 1.0f;
    if( wi.z() < -0.9999f)
        return 0.0f;

    // a
    const float theta_i = acosf(wi.z());
    const float sin_theta_i = sinf(theta_i);

    const float alphai = alpha_i(wi);

    // value
    const float value = 0.5f * (wi.z() + sqrtf(wi.z()*wi.z() + sin_theta_i*sin_theta_i*alphai*alphai));

    return value;
}

Vector2f MicrosurfaceSlopeGGX::sampleP22_11(const float theta_i, const float U, const float U_2) const
{
    Vector2f slope;

    if(theta_i < 0.0001f)
    {
        const float r = sqrtf(U/(1.0f-U));
        const float phi = 6.28318530718f * U_2;
        slope.x() = r * cosf(phi);
        slope.y() = r * sinf(phi);
        return slope;
    }

    // constant
    const float sin_theta_i = sinf(theta_i);
    const float cos_theta_i = cosf(theta_i);
    const float tan_theta_i = sin_theta_i/cos_theta_i;

    // slope associated to theta_i
    const float slope_i = cos_theta_i/sin_theta_i;

    // projected area
    const float projectedarea = 0.5f * (cos_theta_i + 1.0f);
    if(projectedarea < 0.0001f || projectedarea!=projectedarea)
        return Vector2f(0,0);
    // normalization coefficient
    const float c = 1.0f / projectedarea;

    const float A = 2.0f*U/cos_theta_i/c - 1.0f;
    const float B = tan_theta_i;
    const float tmp = 1.0f / (A*A-1.0f);

    const float D = sqrtf(std::max(0.0f, B*B*tmp*tmp - (A*A-B*B)*tmp));
    const float slope_x_1 = B*tmp - D;
    const float slope_x_2 = B*tmp + D;
    slope.x() = (A < 0.0f || slope_x_2 > 1.0f/tan_theta_i) ? slope_x_1 : slope_x_2;

    float U2;
    float S;
    if(U_2 > 0.5f)
    {
        S = 1.0f;
        U2 = 2.0f*(U_2-0.5f);
    }
    else
    {
        S = -1.0f;
        U2 = 2.0f*(0.5f-U_2);
    }
    const float z = (U2*(U2*(U2*0.27385f-0.73369f)+0.46341f)) / (U2*(U2*(U2*0.093073f+0.309420f)-1.000000f)+0.597999f);
    slope.y() = S * z * sqrtf(1.0f+slope.x()*slope.x());

    return slope;
}

/************* MICROSURFACE *************/

Microsurface::Microsurface(const PropertyList &props) {
    bool height_uniform = props.getBoolean("height_uniform", true);
    bool slope_beckmann = props.getBoolean("slope_beckmann", false);
    float alpha_x = props.getFloat("alpha");
    float alpha_y = props.getFloat("alpha");
    if (height_uniform) {
        m_microsurfaceheight = std::make_unique<MicrosurfaceHeightUniform>();
    } else {
        m_microsurfaceheight = std::make_unique<MicrosurfaceHeightGaussian>();
    }

    if (slope_beckmann) {
        m_microsurfaceslope = std::make_unique<MicrosurfaceSlopeBeckmann>(alpha_x, alpha_y);
    } else {
        m_microsurfaceslope = std::make_unique<MicrosurfaceSlopeGGX>(alpha_x, alpha_y);
    }
}

float Microsurface::G_1(const Vector3f& wi) const
{
    if(wi.z() > 0.9999f)
        return 1.0f;
    if(wi.z() <= 0.0f)
        return 0.0f;

    // Lambda
    const float Lambda = m_microsurfaceslope->Lambda(wi);
    // value
    const float value = 1.0f / (1.0f + Lambda);
    return value;
}

float Microsurface::G_1(const Vector3f& wi, const float h0) const
{
    if(wi.z() > 0.9999f)
        return 1.0f;
    if(wi.z() <= 0.0f)
        return 0.0f;

    // height CDF
    const float C1_h0 = m_microsurfaceheight->C1(h0);
    // Lambda
    const float Lambda = m_microsurfaceslope->Lambda(wi);
    // value
    const float value = powf(C1_h0, Lambda);
    return value;
}

float Microsurface::sampleHeight(const Vector3f& wr, const float hr, const float U) const
{
    if(wr.z() > 0.9999f)
        return std::numeric_limits<float>::max();
    if(wr.z() < -0.9999f)
    {
        const float value = m_microsurfaceheight->invC1(U*m_microsurfaceheight->C1(hr));
        return value;
    }
    if(fabsf(wr.z()) < 0.0001f)
        return hr;

    // probability of intersection
    const float G_1_ = G_1(wr, hr);

    if (U > 1.0f - G_1_) // leave the microsurface
        return std::numeric_limits<float>::max();

    const float h = m_microsurfaceheight->invC1(
            m_microsurfaceheight->C1(hr) / powf((1.0f-U),1.0f/m_microsurfaceslope->Lambda(wr))
    );
    return h;
}

Vector3f Microsurface::sample(const Vector3f& wi, int& scatteringOrder, Sampler *sampler) const
{
    // init
    Vector3f wr = -wi;
    float hr = 1.0f + m_microsurfaceheight->invC1(0.999f);

    // random walk
    scatteringOrder = 0;
    while(true)
    {
        // next height
        float U = sampler->next1D();
        hr = sampleHeight(wr, hr, U);

        // leave the microsurface?
        if( hr == std::numeric_limits<float>::max() )
            break;
        else
            scatteringOrder++;

        // next direction
        wr = samplePhaseFunction(-wr, sampler);

        // if NaN (should not happen, just in case)
        if( (hr != hr) || (wr.z() != wr.z()) )
            return Vector3f(0,0,1);
    }

    return wr;
}

float Microsurface::eval(const Vector3f& wi, const Vector3f& wo, Sampler *sampler, const int scatteringOrder) const
{
    if(wo.z() < 0)
        return 0;
    // init
    Vector3f wr = -wi;
    float hr = 1.0f + m_microsurfaceheight->invC1(0.999f);

    float sum = 0;

    // random walk
    int current_scatteringOrder = 0;
    while(scatteringOrder==0 || current_scatteringOrder <= scatteringOrder)
    {
        // next height
        float U = sampler->next1D();
        hr = sampleHeight(wr, hr, U);

        // leave the microsurface?
        if( hr == std::numeric_limits<float>::max())
            break;
        else
            current_scatteringOrder++;

        // next event estimation
        float phasefunction = evalPhaseFunction(-wr, wo, sampler);
        float shadowing = G_1(wo, hr);
        float I = phasefunction * shadowing;

        if ( IsFiniteNumber(I) && (scatteringOrder==0 || current_scatteringOrder==scatteringOrder) )
            sum += I;

        // next direction
        wr = samplePhaseFunction(-wr, sampler);

        // if NaN (should not happen, just in case)
        if( (hr != hr) || (wr.z() != wr.z()) )
            return 0.0f;
    }

    return sum;
}

float MicrosurfaceConductor::evalPhaseFunction(const Vector3f& wi, const Vector3f& wo, Sampler *sampler) const
{
    // half vector
    const Vector3f wh = (wi+wo).normalized();
    if(wh.z() < 0.0f)
        return 0.0f;

    // value
    const float value = 0.25f * m_microsurfaceslope->D_wi(wi, wh) / wi.dot(wh);
    return value;
}

Vector3f MicrosurfaceConductor::samplePhaseFunction(const Vector3f& wi, Sampler *sampler) const
{
    const float U1 = sampler->next1D();
    const float U2 = sampler->next1D();

    Vector3f wm = m_microsurfaceslope->sampleD_wi(wi, U1, U2);

    // reflect
    const Vector3f wo = -wi + 2.0f * wm * wi.dot(wm);

    return wo;
}

float MicrosurfaceConductor::evalSingleScattering(const Vector3f& wi, const Vector3f& wo, Sampler *sampler) const
{
    // half-vector
    const Vector3f wh = (wi+wo).normalized();
    const float D = m_microsurfaceslope->D(wh);

    // masking-shadowing
    const float G2 = 1.0f / (1.0f + m_microsurfaceslope->Lambda(wi) + m_microsurfaceslope->Lambda(wo));

    // BRDF * cos
    const float value = D * G2 / (4.0f * wi.z());

    return value;
}

Vector3f MicrosurfaceDielectric::refract(const Vector3f &wi, const Vector3f &wm, const float eta) const
{
    const float cos_theta_i = wi.dot(wm);
    const float cos_theta_t2 = 1.0f - (1.0f-cos_theta_i*cos_theta_i) / (eta*eta);
    const float cos_theta_t = -sqrtf(std::max(0.0f,cos_theta_t2));

    return wm * (wi.dot(wm) / eta + cos_theta_t) - wi / eta;
}


float MicrosurfaceDielectric::Fresnel(const Vector3f& wi, const Vector3f& wm, const float eta) const
{
    const float cos_theta_i = wi.dot(wm);
    const float cos_theta_t2 = 1.0f - (1.0f-cos_theta_i*cos_theta_i) / (eta*eta);

    // total internal reflection
    if (cos_theta_t2 <= 0.0f) return 1.0f;

    const float cos_theta_t = sqrtf(cos_theta_t2);

    const float Rs = (cos_theta_i - eta * cos_theta_t) / (cos_theta_i + eta * cos_theta_t);
    const float Rp = (eta * cos_theta_i - cos_theta_t) / (eta * cos_theta_i + cos_theta_t);

    const float F = 0.5f * (Rs * Rs + Rp * Rp);
    return F;
}

// wrapper (only for the API and testing)
float MicrosurfaceDielectric::evalPhaseFunction(const Vector3f& wi, const Vector3f& wo, Sampler *sampler) const
{
    return evalPhaseFunction(wi, wo, true, true, sampler) + evalPhaseFunction(wi, wo, true, false, sampler);
}

float MicrosurfaceDielectric::evalPhaseFunction(const Vector3f& wi, const Vector3f& wo, const bool wi_outside, const bool wo_outside, Sampler *sampler) const
{
    const float eta = wi_outside ? m_eta : 1.0f / m_eta;

    if( wi_outside == wo_outside ) // reflection
    {
        // half vector
        const Vector3f wh = (wi+wo).normalized();
        // value
        const float value = (wi_outside) ?
                            (0.25f * m_microsurfaceslope->D_wi(wi, wh) / wi.dot(wh) * Fresnel(wi, wh, eta)) :
                            (0.25f * m_microsurfaceslope->D_wi(-wi, -wh) / wi.dot(wh) * Fresnel(-wi, -wh, eta)) ;
        return value;
    }
    else // transmission
    {
        Vector3f wh = -(wi+wo*eta).normalized();
        wh *= (wi_outside) ? (sign(wh.z())) : (-sign(wh.z()));

        if(wh.dot(wi) < 0)
            return 0;

        float value;
        if(wi_outside){
            value = eta*eta * (1.0f-Fresnel(wi, wh, eta)) *
                    m_microsurfaceslope->D_wi(wi, wh) * std::max(0.0f, -wo.dot(wh)) *
                    1.0f / powf(wi.dot(wh)+eta*wo.dot(wh), 2.0f);
        }
        else
        {
            value = eta*eta * (1.0f-Fresnel(-wi, -wh, eta)) *
                    m_microsurfaceslope->D_wi(-wi, -wh) * std::max(0.0f, -wo.dot(wh)) *
                    1.0f / powf(wi.dot(wh)+eta*wo.dot(wh), 2.0f);
        }

        return value;
    }
}

Vector3f MicrosurfaceDielectric::samplePhaseFunction(const Vector3f& wi, Sampler *sampler) const
{
    bool wo_outside;
    return samplePhaseFunction(wi, wi.z() > 0, wo_outside, sampler);
}

Vector3f MicrosurfaceDielectric::samplePhaseFunction(const Vector3f& wi, const bool wi_outside, bool& wo_outside, Sampler *sampler) const
{
    const float U1 = sampler->next1D();
    const float U2 = sampler->next1D();

    const float eta = wi_outside ? m_eta : 1.0f / m_eta;

    Vector3f wm = wi_outside ? (m_microsurfaceslope->sampleD_wi(wi, U1, U2)) :
              (-m_microsurfaceslope->sampleD_wi(-wi, U1, U2)) ;

    const float F = Fresnel(wi, wm, eta);

    if( sampler->next1D() < F )
    {
        const Vector3f wo = -wi + 2.0f * wm * wi.dot(wm); // reflect
        return wo;
    }
    else
    {
        wo_outside = !wi_outside;
        const Vector3f wo = refract(wi, wm, eta);
        return (wo).normalized();
    }
}

float MicrosurfaceDielectric::evalSingleScattering(const Vector3f& wi, const Vector3f& wo, Sampler *sampler) const
{
    bool wi_outside = true;
    bool wo_outside = wo.z() > 0;

    const float eta = m_eta;

    if(wo_outside) // reflection
    {
        // D
        const Vector3f wh = (wi+wo).normalized();
        const float D = m_microsurfaceslope->D(wh);

        // masking shadowing
        const float Lambda_i = m_microsurfaceslope->Lambda(wi);
        const float Lambda_o = m_microsurfaceslope->Lambda(wo);
        const float G2 = 1.0f / (1.0f + Lambda_i + Lambda_o);

        // BRDF
        const float value = Fresnel(wi, wh, eta) * D * G2 / (4.0f * wi.z());
        return value;
    }
    else // refraction
    {
        // D
        Vector3f wh = -(wi+wo*eta).normalized();
        if(eta<1.0f)
            wh = -wh;
        const float D = m_microsurfaceslope->D(wh);

        // G2
        const float Lambda_i = m_microsurfaceslope->Lambda(wi);
        const float Lambda_o = m_microsurfaceslope->Lambda(-wo);
        const float G2 = (float) beta(1.0f+Lambda_i, 1.0f+Lambda_o);

        // BSDF
        const float value = std::max(0.0f, wi.dot(wh)) * std::max(0.0f, -wo.dot(wh)) *
                            1.0f / wi.z() * eta*eta * (1.0f-Fresnel(wi, wh, eta)) *
                            G2 * D / powf(wi.dot(wh)+eta*wo.dot(wh), 2.0f);
        return value;
    }
}

float MicrosurfaceDielectric::eval(const Vector3f& wi, const Vector3f& wo, Sampler *sampler, const int scatteringOrder) const
{
    // init
    Vector3f wr = -wi;
    float hr = 1.0f + m_microsurfaceheight->invC1(0.999f);
    bool outside = wi.z() > 0;

    float sum = 0.0f;

    // random walk
    int current_scatteringOrder = 0;
    while(scatteringOrder==0 || current_scatteringOrder <= scatteringOrder)
    {
        // next height
        float U = sampler->next1D();
        hr = (outside) ? sampleHeight(wr, hr, U) : -sampleHeight(-wr, -hr, U);

        // leave the microsurface?
        if( hr == std::numeric_limits<float>::max() || hr == -std::numeric_limits<float>::max())
            break;
        else
            current_scatteringOrder++;

        // next event estimation
        float phasefunction = evalPhaseFunction(-wr, wo, outside, (wo.z()>0), sampler);
        float shadowing = (wo.z()>0) ? G_1(wo, hr) : G_1(-wo, -hr);
        float I = phasefunction * shadowing;

        if ( IsFiniteNumber(I) && (scatteringOrder==0 || current_scatteringOrder==scatteringOrder) )
            sum += I;

        // next direction
        wr = samplePhaseFunction(-wr, outside, outside, sampler);

        // if NaN (should not happen, just in case)
        if( (hr != hr) || (wr.z() != wr.z()) )
            return 0.0f;
    }

    return sum;
}

Vector3f MicrosurfaceDielectric::sample(const Vector3f& wi, int& scatteringOrder, Sampler *sampler) const
{
    // init
    Vector3f wr = -wi;
    float hr = 1.0f + m_microsurfaceheight->invC1(0.999f);
    bool outside = wi.z() > 0;

    // random walk
    scatteringOrder = 0;
    while(true)
    {
        // next height
        float U = sampler->next1D();
        hr = (outside) ? sampleHeight(wr, hr, U) : -sampleHeight(-wr, -hr, U);

        // leave the microsurface?
        if( hr == std::numeric_limits<float>::max() || hr == -std::numeric_limits<float>::max())
            break;
        else
            scatteringOrder++;

        // next direction
        wr = samplePhaseFunction(-wr, outside, outside, sampler);

        // if NaN (should not happen, just in case)
        if( (hr != hr) || (wr.z() != wr.z()) )
            return Vector3f(0,0,1);
    }

    return wr;
}

float MicrosurfaceDiffuse::evalPhaseFunction(const Vector3f& wi, const Vector3f& wo, Sampler *sampler) const
{
    const float U1 = sampler->next1D();
    const float U2 = sampler->next1D();
    Vector3f wm = m_microsurfaceslope->sampleD_wi(wi, U1, U2);

    // value
    const float value = 1.0f/M_PI * std::max(0.0f, wo.dot(wm));
    return value;
}

// build orthonormal basis (Building an Orthonormal Basis from a 3D Unit Vector Without Normalization, [Frisvad2012])
void buildOrthonormalBasis(Vector3f& omega_1, Vector3f& omega_2, const Vector3f& omega_3)
{
    if(omega_3.z() < -0.9999999f)
    {
        omega_1 = Vector3f ( 0.0f , -1.0f , 0.0f );
        omega_2 = Vector3f ( -1.0f , 0.0f , 0.0f );
    } else {
        const float a = 1.0f /(1.0f + omega_3.z() );
        const float b = -omega_3.x()*omega_3 .y()*a ;
        omega_1 = Vector3f (1.0f - omega_3.x()*omega_3. x()*a , b , -omega_3.x() );
        omega_2 = Vector3f (b , 1.0f - omega_3.y()*omega_3.y()*a , -omega_3.y() );
    }
}


Vector3f MicrosurfaceDiffuse::samplePhaseFunction(const Vector3f& wi, Sampler *sampler) const
{
    const float U1 = sampler->next1D();
    const float U2 = sampler->next1D();
    const float U3 = sampler->next1D();
    const float U4 = sampler->next1D();

    Vector3f wm = m_microsurfaceslope->sampleD_wi(wi, U1, U2);

    // sample diffuse reflection
    Vector3f w1, w2;
    buildOrthonormalBasis(w1, w2, wm);

    float r1 = 2.0f*U3 - 1.0f;
    float r2 = 2.0f*U4 - 1.0f;

    // concentric map code from
    // http://psgraphics.blogspot.ch/2011/01/improved-code-for-concentric-map.html
    float phi, r;
    if (r1 == 0 && r2 == 0) {
        r = phi = 0;
    } else if (r1*r1 > r2*r2) {
        r = r1;
        phi = (M_PI/4.0f) * (r2/r1);
    } else {
        r = r2;
        phi = (M_PI/2.0f) - (r1/r2) * (M_PI/4.0f);
    }
    float x = r*cosf(phi);
    float y = r*sinf(phi);
    float z = sqrtf(std::max(0.0f, 1.0f - x*x - y*y));
    Vector3f wo = x*w1 + y*w2 + z*wm;

    return wo;
}

// stochastic evaluation
// Heitz and Dupuy 2015
// Implementing a Simple Anisotropic Rough Diffuse Material with Stochastic Evaluation
float MicrosurfaceDiffuse::evalSingleScattering(const Vector3f& wi, const Vector3f& wo, Sampler *sampler) const
{
    // sample visible microfacet
    const float U1 = sampler->next1D();
    const float U2 = sampler->next1D();
    const Vector3f wm = m_microsurfaceslope->sampleD_wi(wi, U1, U2);

    // shadowing given masking
    const float Lambda_i = m_microsurfaceslope->Lambda(wi);
    const float Lambda_o = m_microsurfaceslope->Lambda(wo);
    float G2_given_G1 = (1.0f + Lambda_i) / (1.0f + Lambda_i + Lambda_o);

    // evaluate diffuse and shadowing given masking
    const float value = 1.0f / (float)M_PI * std::max(0.0f, wm.dot(wo)) * G2_given_G1;

    return value;
}

NORI_NAMESPACE_END


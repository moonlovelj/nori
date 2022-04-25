//
// Created by 郭彬 on 2022/4/20.
//
// https://eheitzresearch.wordpress.com/240-2/

#ifndef NORI_MICROFACETMULTIPLESCATTERING_H
#define NORI_MICROFACETMULTIPLESCATTERING_H

#include <memory>

#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/warp.h>
#include <nori/sampler.h>

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
    float pdf(const Vector3f &wi, const Vector3f &wh) const {
        return D_wi(wi, wh) * Frame::cosTheta(wh);
    }
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


    virtual Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample) const override{
        throw NoriException("This function should not be called!");
    }
    Color3f sample(BSDFQueryRecord &bRec, Sampler *sampler) const override {
        if (Frame::cosTheta(bRec.wi) <= 0)
            return Color3f(0.0f);

//        bRec.measure = ESolidAngle;
//
//        /* Warp a uniformly distributed sample on [0,1]^2
//           to a direction on a cosine-weighted hemisphere */
//        bRec.wo = Warp::squareToCosineHemisphere(sampler->next2D());
//
//        /* Relative index of refraction: no change */
//        bRec.eta = 1.0f;
//
//        /* eval() / pdf() * cos(theta) = albedo. There
//           is no need to call these functions. */
//        return Color3f(1);

        int scatteringOrder;
        bRec.wo = sample(bRec.wi, scatteringOrder, sampler);

        //test_masking(bRec.wo, sampler);
        //test_masking_shadowing(bRec.wi, bRec.wo, sampler);
        //test_sample_phasefunction(bRec.wi, sampler);

        if (bRec.wo == Vector3f(0,0,1)) {
            return Color3f(0);
        }
        bRec.measure = ESolidAngle;

        float f = pdf(bRec);
        if (f == 0) return Color3f(0);
        return eval(bRec) / f;
    }

    Color3f eval(const BSDFQueryRecord &bRec) const override {
        if (bRec.sampler == nullptr) throw NoriException("sampler can not be null");
        return bRec.wi.z() < bRec.wo.z() ? eval(bRec.wi, bRec.wo, bRec.sampler) :
               eval(bRec.wo, bRec.wi, bRec.sampler) / Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo);
    }

    virtual std::string toString() const override {
        return "";
    }

    void test_masking(const Vector3f &wo, Sampler *sampler) const
    {
        float G1 = 1.0f / (1.f + m_microsurfaceslope->Lambda(wo));

        const int N = 100000;
        double V = 0;
        for (int i = 0; i < N ; ++i) {
            const float U1 = sampler->next1D();
            const float h = m_microsurfaceheight->invC1(U1);

            const float U2 =  sampler->next1D();
            const float h_wo = sampleHeight(wo, h, U2);

            if (h_wo == std::numeric_limits<float>::max())
                V += 1.0 / N;
        }

        std::cout <<" analytic = " << G1 << endl;
        std::cout <<" stochastic = " << V << endl;
    }

    void test_masking_shadowing(const  Vector3f &wi, const Vector3f &wo, Sampler *sampler) const {
        float G2 = 1.0f / (1.0f + m_microsurfaceslope->Lambda(wi) + m_microsurfaceslope->Lambda(wo));

        const int N = 100000;
        double V = 0;
        for (int i = 0; i < N ; ++i) {
            const float U1 = sampler->next1D();
            const float h = m_microsurfaceheight->invC1(U1);
            const float h_wi = sampleHeight(wi, h, sampler->next1D());
            const float h_wo = sampleHeight(wo, h, sampler->next1D());
            if (h_wi == std::numeric_limits<float>::max() && h_wo == std::numeric_limits<float>::max())
                V += 1.0 / N;
        }

        std::cout <<" analytic = " << G2 << endl;
        std::cout <<" stochastic = " << V << endl;
    }

    void test_norimalization_D_wi(const Vector3f &wi, Sampler *sampler) const
    {

        double value_quadrature = 0;
        for(double theta_m=0 ; theta_m < 0.5*M_PI ; theta_m += 0.005){
            for(double phi_m=0 ; phi_m < 2.0*M_PI ; phi_m += 0.005)
            {
                const Vector3f wm(cos(phi_m)*sin(theta_m), sin(phi_m)*sin(theta_m), cos(theta_m));
                value_quadrature += 0.005*0.005*abs(sin(theta_m)) * (double)m_microsurfaceslope->D_wi(wi, wm);
            }
        }

        // display
        cout << "\\int␣D_wi(wm)␣dwm␣=␣\t\t" << value_quadrature << endl;
    }

    void test_sample_D_wi(const Vector3f &wi, Sampler *sampler) const {
        double quadrature_int_x = 0;
        double quadrature_int_y = 0;
        double quadrature_int_z = 0;
        double quadrature_int_x2 = 0;
        double quadrature_int_y2 = 0;
        double quadrature_int_z2 = 0;
        for(double theta_m=0 ; theta_m < 0.5*M_PI ; theta_m += 0.005) {
            for(double phi_m=0 ; phi_m < 2.0*M_PI ; phi_m += 0.005) {
                const Vector3f wm(cos(phi_m)*sin(theta_m), sin(phi_m)*sin(theta_m), cos(theta_m));
                const double d = 0.005*0.005*abs(sin(theta_m)) * (double)m_microsurfaceslope->D_wi(wi, wm);
                quadrature_int_x += d * wm.x();
                quadrature_int_y += d * wm.y();
                quadrature_int_z += d * wm.z();
                quadrature_int_x2 += d * wm.x() * wm.x();
                quadrature_int_y2 += d * wm.y() * wm.y();
                quadrature_int_z2 += d * wm.z() * wm.z();
            }
        }

        double stochastic_int_x = 0;
        double stochastic_int_y = 0;
        double stochastic_int_z = 0;
        double stochastic_int_x2 = 0;
        double stochastic_int_y2 = 0;
        double stochastic_int_z2 = 0;

        for(int n=0 ; n<100000 ; ++n) {
            const float U1 = sampler->next1D();
            const float U2 = sampler->next1D();
            const Vector3f wm = m_microsurfaceslope->sampleD_wi(wi, U1, U2);
            stochastic_int_x += wm.x() / 100000.0;
            stochastic_int_y += wm.y() / 100000.0;
            stochastic_int_z += wm.z() / 100000.0;
            stochastic_int_x2 += wm.x() * wm.x() / 100000.0;
            stochastic_int_y2 += wm.y() * wm.y() / 100000.0;
            stochastic_int_z2 += wm.z() * wm.z() / 100000.0;
        }

        std::cout << "quadrature␣\\int␣D_wi(wm)␣wm.x␣dwm␣=␣\t\t" << quadrature_int_x<< endl;
        std::cout << "quadrature␣\\int␣D_wi(wm)␣wm.y␣dwm␣=␣\t\t" << quadrature_int_y<< endl;
        std::cout << "quadrature␣\\int␣D_wi(wm)␣wm.z␣dwm␣=␣\t\t" << quadrature_int_z<< endl;
        std::cout << "quadrature␣\\int␣D_wi(wm)␣wm.x^2␣dwm␣=␣\t\t" << quadrature_int_x2 << endl;
        std::cout << "quadrature␣\\int␣D_wi(wm)␣wm.x^2␣dwm␣=␣\t\t" << quadrature_int_y2 << endl;
        std::cout << "quadrature␣\\int␣D_wi(wm)␣wm.x^2␣dwm␣=␣\t\t" << quadrature_int_z2 << endl;
        std::cout << endl;

        std::cout << "stochastic␣\\int␣D_wi(wm)␣wm.x␣dwm␣=␣\t\t" << stochastic_int_x << endl;
        std::cout << "stochastic␣\\int␣D_wi(wm)␣wm.y␣dwm␣=␣\t\t" << stochastic_int_y << endl;
        std::cout << "stochastic␣\\int␣D_wi(wm)␣wm.z␣dwm␣=␣\t\t" << stochastic_int_z << endl;
        std::cout << "stochastic␣\\int␣D_wi(wm)␣wm.x^2␣dwm␣=␣\t\t" << stochastic_int_x2 << endl;
        std::cout << "stochastic␣\\int␣D_wi(wm)␣wm.x^2␣dwm␣=␣\t\t" << stochastic_int_y2 << endl;
        std::cout << "stochastic␣\\int␣D_wi(wm)␣wm.x^2␣dwm␣=␣\t\t" << stochastic_int_z2 << endl;

        std::cout<<" ===**********====\n";
    }

    void test_phase_function(const Vector3f &wi, Sampler *sampler) const {
        // quadrature \int p(wi, wo) dwo
        double value_quadrature = 0;
        for(double theta_o=0 ; theta_o < M_PI ; theta_o += 0.005)
        {
            for(double phi_o=0 ; phi_o < 2.0*M_PI ; phi_o += 0.005)
            {
                const Vector3f wo(cos(phi_o)*sin(theta_o), sin(phi_o)*sin(theta_o), cos(theta_o));
                value_quadrature += 0.005*0.005*abs(sin(theta_o)) * (double)evalPhaseFunction(wi, wo, sampler);
            }
        }

// display
        cout << "\\int␣p(wi,␣wo)␣dwo␣=␣\t\t" << value_quadrature << endl;
    }

    void test_sample_phasefunction(const Vector3f &wi, Sampler *sampler) const {
        if(m_microsurfaceslope->projectedArea(wi) < 0.01f) {
            cout << "Warning:␣the␣projected␣area␣is␣too␣small" << endl;
            cout << "the␣ray␣cannot␣intersect␣the␣microsurface␣in␣this␣configuration" << endl;
            cout << "and␣the␣phase␣function␣should␣not␣be␣called." << endl << endl;
            throw NoriException("test_sample_phasefunction error \n");
        }
// quadrature with eval p(wi, wo)
        double quadrature_int_x = 0;
        double quadrature_int_y = 0;
        double quadrature_int_z = 0;
        double quadrature_int_x2 = 0;
        double quadrature_int_y2 = 0;
        double quadrature_int_z2 = 0;

        for(double theta_o=0 ; theta_o < M_PI ; theta_o += 0.005){
            for(double phi_o=0 ; phi_o < 2.0*M_PI ; phi_o += 0.005) {
                const Vector3f wo(cos(phi_o)*sin(theta_o), sin(phi_o)*sin(theta_o), cos(theta_o));
                const double d = 0.005*0.005*abs(sin(theta_o)) * (double)evalPhaseFunction(wi, wo, sampler);
                quadrature_int_x += d * wo.x();
                quadrature_int_y += d * wo.y();
                quadrature_int_z += d * wo.z();
                quadrature_int_x2 += d * wo.x() * wo.x();
                quadrature_int_y2 += d * wo.y() * wo.y();
                quadrature_int_z2 += d * wo.z() * wo.z();


            }
        }

// sampling \p(wi, wo)
        const int N = 100000;
        double stochastic_int_x = 0;
        double stochastic_int_y = 0;
        double stochastic_int_z = 0;
        double stochastic_int_x2 = 0;
        double stochastic_int_y2 = 0;
        double stochastic_int_z2 = 0;

        for(int n=0 ; n<N ; ++n) {
            const Vector3f wo = samplePhaseFunction(wi, sampler);
            stochastic_int_x += wo.x() / (double) N;
            stochastic_int_y += wo.y() / (double) N;
            stochastic_int_z += wo.z() / (double) N;
            stochastic_int_x2 += wo.x() * wo.x() / (double) N;
            stochastic_int_y2 += wo.y() * wo.y() / (double) N;
            stochastic_int_z2 += wo.z() * wo.z() / (double) N;
        }

// display
        cout << "quadrature␣\\int␣p(wi,␣wo)␣wo.x␣dwo␣=␣\t\t" << quadrature_int_x << endl;
        cout << "quadrature␣\\int␣p(wi,␣wo)␣wo.y␣dwo␣=␣\t\t" << quadrature_int_y << endl;
        cout << "quadrature␣\\int␣p(wi,␣wo)␣wo.z␣dwo␣=␣\t\t" << quadrature_int_z << endl;
        cout << "quadrature␣\\int␣p(wi,␣wo)␣wo.x^2␣dwo␣=␣\t" << quadrature_int_x2 << endl;
        cout << "quadrature␣\\int␣p(wi,␣wo)␣wo.x^2␣dwo␣=␣\t" << quadrature_int_y2 << endl;
        cout << "quadrature␣\\int␣p(wi,␣wo)␣wo.x^2␣dwo␣=␣\t" << quadrature_int_z2 << endl;
        cout << endl;
        cout << "stochastic␣\\int␣p(wi,␣wo)␣wo.x␣dwo␣=␣\t\t" << stochastic_int_x << endl;
        cout << "stochastic␣\\int␣p(wi,␣wo)␣wo.y␣dwo␣=␣\t\t" << stochastic_int_y << endl;
        cout << "stochastic␣\\int␣p(wi,␣wo)␣wo.z␣dwo␣=␣\t\t" << stochastic_int_z << endl;
        cout << "stochastic␣\\int␣p(wi,␣wo)␣wo.x^2␣dwo␣=␣\t" << stochastic_int_x2 << endl;
        cout << "stochastic␣\\int␣p(wi,␣wo)␣wo.x^2␣dwo␣=␣\t" << stochastic_int_y2 << endl;
        cout << "stochastic␣\\int␣p(wi,␣wo)␣wo.x^2␣dwo␣=␣\t" << stochastic_int_z2 << endl;
        cout << " *****=======********\n";
    }

protected:
    // masking function
    float G_1(const Vector3f& wi) const;
    // masking function at height h0
    float G_1(const Vector3f& wi, const float h0) const;
    // sample height in outgoing direction
    float sampleHeight(const Vector3f& wo, const float h0, const float U) const;

    // evaluate BSDF with a random walk (stochastic but unbiased)
    // scatteringOrder=0 --> contribution from all scattering events
    // scatteringOrder=1 --> contribution from 1st bounce only
    // scatteringOrder=2 --> contribution from 2nd bounce only, etc..
    virtual float eval(const Vector3f& wi, const Vector3f& wo, Sampler *sampler, const int scatteringOrder=0) const;

    // sample BSDF with a random walk
    // scatteringOrder is set to the number of bounces computed for this sample
    virtual Vector3f sample(const Vector3f& wi, int& scatteringOrder, Sampler *sampler) const;
    Vector3f sample(const Vector3f& wi, Sampler *sampler) const {int scatteringOrder; return sample(wi, scatteringOrder, sampler);}

    // evaluate local phase function 
    virtual float evalPhaseFunction(const Vector3f& wi, const Vector3f& wo, Sampler *sampler) const=0;
    // sample local phase function
    virtual Vector3f samplePhaseFunction(const Vector3f& wi, Sampler *sampler) const=0;

    // evaluate BSDF limited to single scattering 
    // this is in average equivalent to eval(wi, wo, 1);
    virtual float evalSingleScattering(const Vector3f& wi, const Vector3f& wo, Sampler *sampler) const=0;

    // height distribution
    std::unique_ptr<const MicrosurfaceHeight> m_microsurfaceheight;
    // slope distribution
    std::unique_ptr<const MicrosurfaceSlope> m_microsurfaceslope;
};

/* Microsurface made of conductor material */
class MicrosurfaceConductor : public Microsurface
{
public:
    MicrosurfaceConductor(const PropertyList & props) : Microsurface(props)
    {
    }

public:
    // evaluate local phase function
    virtual float evalPhaseFunction(const Vector3f& wi, const Vector3f& wo, Sampler *sampler) const override ;
    // sample local phase function
    virtual Vector3f samplePhaseFunction(const Vector3f& wi, Sampler *sampler) const override ;

    // evaluate BSDF limited to single scattering
    // this is in average equivalent to eval(wi, wo, 1);
    virtual float evalSingleScattering(const Vector3f& wi, const Vector3f& wo, Sampler *sampler) const override ;

    float pdf(const BSDFQueryRecord &bRec) const override {
        return 0;
    }
};
NORI_REGISTER_CLASS(MicrosurfaceConductor, "mulmicrofacetconductor");

/* Microsurface made of conductor material */
class MicrosurfaceDielectric : public Microsurface
{
public:
    float m_eta;
public:
    MicrosurfaceDielectric(const PropertyList & props) : Microsurface(props)
    {
        /* Interior IOR (default: BK7 borosilicate optical glass) */
        const float inIOR = props.getFloat("intIOR", 1.5046f);

        /* Exterior IOR (default: air) */
        const float exIOR = props.getFloat("extIOR", 1.000277f);

        m_eta = inIOR / exIOR;
    }
//    MicrosurfaceDielectric(const bool height_uniform, // uniform or Gaussian
//                           const bool slope_beckmann, // Beckmann or GGX
//                           const float alpha_x,
//                           const float alpha_y,
//                           const float eta = 1.5f)
//            : Microsurface(height_uniform, slope_beckmann, alpha_x, alpha_y),
//              m_eta(eta)
//    {}

    // evaluate BSDF with a random walk (stochastic but unbiased)
    // scatteringOrder=0 --> contribution from all scattering events
    // scatteringOrder=1 --> contribution from 1st bounce only
    // scatteringOrder=2 --> contribution from 2nd bounce only, etc..
    virtual float eval(const Vector3f& wi, const Vector3f& wo, Sampler *sampler, const int scatteringOrder=0) const override ;

    // sample final BSDF with a random walk
    // scatteringOrder is set to the number of bounces computed for this sample
    virtual Vector3f sample(const Vector3f& wi, int& scatteringOrder, Sampler *sampler) const override ;

    float pdf(const BSDFQueryRecord &bRec) const override {
        return 0;
    }

protected:
    // evaluate local phase function
    virtual float evalPhaseFunction(const Vector3f& wi, const Vector3f& wo, Sampler *sampler) const override ;
    float evalPhaseFunction(const Vector3f& wi, const Vector3f& wo, const bool wi_outside, const bool wo_outside, Sampler *sampler) const;
    // sample local phase function
    virtual Vector3f samplePhaseFunction(const Vector3f& wi, Sampler *sampler) const override ;
    Vector3f samplePhaseFunction(const Vector3f& wi, const bool wi_outside, bool& wo_outside, Sampler *sampler) const;

    // evaluate BSDF limited to single scattering
    // this is in average equivalent to eval(wi, wo, 1);
    virtual float evalSingleScattering(const Vector3f& wi, const Vector3f& wo, Sampler *sampler) const override ;

    float Fresnel(const Vector3f& wi, const Vector3f& wm, const float eta) const;
    Vector3f refract(const Vector3f &wi, const Vector3f &wm, const float eta) const;
};
NORI_REGISTER_CLASS(MicrosurfaceDielectric, "mulmicrofacetdielectric");

/* Microsurface made of conductor material */
class MicrosurfaceDiffuse : public Microsurface
{
public:
    MicrosurfaceDiffuse(const PropertyList & props) : Microsurface(props)
    {}

    float pdf(const BSDFQueryRecord &bRec) const override {
        if (bRec.measure != ESolidAngle) {
            return 0.f;
        }
        if (Frame::cosTheta(bRec.wi) <= 0 ||
            Frame::cosTheta(bRec.wo) <= 0) {
            return 0.f;
        }
        return Frame::cosTheta(bRec.wo) * M_1_PI;
    }

protected:
    // evaluate local phase function
    virtual float evalPhaseFunction(const Vector3f& wi, const Vector3f& wo, Sampler *sampler) const override ;
    // sample local phase function
    virtual Vector3f samplePhaseFunction(const Vector3f& wi, Sampler *sampler) const override ;

    // evaluate BSDF limited to single scattering
    // this is in average equivalent to eval(wi, wo, 1);
    virtual float evalSingleScattering(const Vector3f& wi, const Vector3f& wo, Sampler *sampler) const override ;
};

NORI_REGISTER_CLASS(MicrosurfaceDiffuse, "mulmicrofacetdiffuse");
NORI_NAMESPACE_END

#endif //NORI_MICROFACETMULTIPLESCATTERING_H

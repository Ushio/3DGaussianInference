#include "pr.hpp"
#include <iostream>
#include <memory>
#include <ppl.h>


#define SPLAT_BOUNDS 3.0f
#define MIN_THROUGHPUT ( 1.0f / 256.0f )

template <class T>
inline T ss_max(T x, T y)
{
    return (x < y) ? y : x;
}

template <class T>
inline T ss_min(T x, T y)
{
    return (y < x) ? y : x;
}

float sign_of(float v)
{
    return v < 0.0f ? -1.0f : 1.0f;
}

// ax^2 + bx + c == 0
int solve_quadratic(float xs[2], float a, float b, float c)
{
    float det = b * b - 4.0f * a * c;
    if (det < 0.0f)
    {
        return 0;
    }

    float k = (-b - sign_of(b) * std::sqrtf(det)) / 2.0f;
    float x0 = k / a;
    float x1 = c / k;
    xs[0] = ss_min(x0, x1);
    xs[1] = ss_max(x0, x1);
    return 2;
}
auto sqr = [](float x) { return x * x; };

float intersect_ray_ellipsoid(glm::vec3 U, glm::vec3 V, glm::vec3 W, glm::vec3 ro, glm::vec3 rd)
{
    glm::vec3 u = U / glm::dot(U, U);
    glm::vec3 v = V / glm::dot(V, V);
    glm::vec3 w = W / glm::dot(W, W);

    float k = 1.0f;
    float t_delta = -glm::dot(rd, ro) / glm::dot(rd, rd);
    glm::vec3 ro_prime = ro + rd * t_delta;

    float urd = glm::dot(u, rd);
    float vrd = glm::dot(v, rd);
    float wrd = glm::dot(w, rd);
    float uro = glm::dot(u, ro_prime);
    float vro = glm::dot(v, ro_prime);
    float wro = glm::dot(w, ro_prime);
    float A = sqr(urd) + sqr(vrd) + sqr(wrd);
    float B = 2.0f * (urd * uro + vrd * vro + wrd * wro);
    float C = sqr(uro) + sqr(vro) + sqr(wro) - k * k;

    float xs[2];
    if (solve_quadratic(xs, A, B, C))
    {
        return xs[0] + t_delta;
    }
    return -1.0f;
}

float lengthSquared(glm::vec3 v)
{
    return v.x * v.x + v.y * v.y + v.z * v.z;
}
// lambda0 is larger
void eignValues(float* lambda0, float* lambda1, float* determinant, const glm::mat2& mat)
{
    float mean = (mat[0][0] + mat[1][1]) * 0.5f;
    float det = mat[0][0] * mat[1][1] - mat[1][0] * mat[0][1];
    float d = std::sqrtf(ss_max(mean * mean - det, 0.0f));
    *lambda0 = mean + d;
    *lambda1 = mean - d;
    *determinant = det;
}

void eigenVectors_of_symmetric(glm::vec2* eigen0, glm::vec2* eigen1, const glm::mat2& m, float lambda)
{
    float s11 = m[0][0];
    float s22 = m[1][1];
    float s12 = m[1][0];

    // to workaround lambda0 == lambda1
    float eps = 1e-15f;
    glm::vec2 e0 = glm::normalize(s11 < s22 ? glm::vec2(s12 + eps, lambda - s11) : glm::vec2(lambda - s22, s12 + eps));
    glm::vec2 e1 = { -e0.y, e0.x };
    *eigen0 = e0;
    *eigen1 = e1;
}

inline float ss_sqrt(float x) {
    float y;
    _mm_store_ss(&y, _mm_sqrt_ss(_mm_load_ss(&x)));
    return y;
}

float min_of_abs(float x, float y, float z)
{
    return ss_min(ss_min(glm::abs(x), glm::abs(y)), glm::abs(z));
}
float max_of_abs(float x, float y)
{
    return ss_max(glm::abs(x), glm::abs(y));
}
void eigen_decomposition(glm::vec3 es[3], float lambdas[3], float A_00, float A_01, float A_02, float A_11, float A_12, float A_22)
{
    glm::vec3 wbasisXY[2] = { {1, 0, 0}, {0, 1, 0} };

    auto sincostan = [](float* s, float* c, float* t, float cos_2theta, float sin_2theta)
    {
        float X = cos_2theta;
        float Y = sin_2theta;
        float YY = Y * Y;
        float L = ss_sqrt(ss_max(X * X + YY, FLT_MIN) /* handle the singularity */);

        // The half vector
        float hx = X + sign_of(X) * L;
        float hy = Y;
        float hyhy = YY;
        float lh = ss_sqrt(hx * hx + hyhy);

        float tanTheta = hy / hx;
        float cosTheta = hx / lh;
        *t = tanTheta;
        *c = cosTheta;
        *s = cosTheta * tanTheta;
    };

    for (int i = 0; i < 32; i++)
    {
        {
            float b = A_12;
            float a = A_11;
            float d = A_22;

            // b is small enough not to affect the next rotation
            if (A_02 + b != A_02)
            {
                float c;
                float s;
                float t;
                sincostan(&s, &c, &t, 0.5f * (d - a), b);

                //glm::mat3 P = glm::mat3(
                //    1, 0, 0,
                //    0, c, -s,
                //    0, s, c
                //);

                for (int j = 0; j < 2; j++)
                {
                    float Y = wbasisXY[j].y;
                    float Z = wbasisXY[j].z;
                    wbasisXY[j].y = c * Y - s * Z;
                    wbasisXY[j].z = s * Y + c * Z;
                }

                float nA_01 = c * A_01 /* - s * A_02*/;
                float nA_02 = /* c * A_02 */ +s * A_01;

                if (i == 0)
                {
                    nA_01 += -s * A_02;
                    nA_02 += c * A_02;
                }

                // simplified via nA_12 == 0
                float nA_11 = A_11 - t * A_12;
                float nA_22 = A_22 + t * A_12;

                float nA_12 = 0.0f; // focus
                A_01 = nA_01;
                A_02 = nA_02;
                A_11 = nA_11;
                A_12 = nA_12;
                A_22 = nA_22;


                // Converged when the smallest diag can not get affected from the largest non-diag
                float minDiag = min_of_abs(A_00, A_11, A_22);
                if (minDiag + max_of_abs(A_02, A_01) == minDiag)
                {
                    break;
                }
            }
        }

        {
            float b = A_01;
            float a = A_00;
            float d = A_11;

            // b is small enough not to affect the next rotation
            if (A_12 + b != A_12)
            {
                float c;
                float s;
                float t;
                sincostan(&s, &c, &t, 0.5f * (d - a), b);

                //glm::mat3 P = glm::mat3(
                //    c, -s, 0,
                //    s, c, 0,
                //    0, 0, 1
                //);

                for (int j = 0; j < 2; j++)
                {
                    float X = wbasisXY[j].x;
                    float Y = wbasisXY[j].y;
                    wbasisXY[j].x = c * X - s * Y;
                    wbasisXY[j].y = s * X + c * Y;
                }

                // simplified via nA_12 == 0
                float nA_00 = A_00 - t * A_01;
                float nA_11 = A_11 + t * A_01;

                float nA_01 = 0.0f; // focus
                float nA_02 = c * A_02 /*- s * A_12*/;
                float nA_12 = /* c * A_12 + */ s * A_02;
                A_00 = nA_00;
                A_01 = nA_01;
                A_02 = nA_02;
                A_11 = nA_11;
                A_12 = nA_12;

                // Converged when the smallest diag can not get affected from the largest non-diag
                float minDiag = min_of_abs(A_00, A_11, A_22);
                if (minDiag + max_of_abs(A_02, A_12) == minDiag)
                {
                    break;
                }
            }
        }

        {
            float b = A_02;
            float a = A_00;
            float d = A_22;
            // b is small enough not to affect the next rotation
            if (A_01 + b != A_01)
            {
                float c;
                float s;
                float t;
                sincostan(&s, &c, &t, 0.5f * (d - a), b);

                //glm::mat3 P = glm::mat3(
                //    c, 0, -s,
                //    0, 1, 0,
                //    s, 0, c
                //);

                for (int j = 0; j < 2; j++)
                {
                    float X = wbasisXY[j].x;
                    float Z = wbasisXY[j].z;
                    wbasisXY[j].x = c * X - s * Z;
                    wbasisXY[j].z = s * X + c * Z;
                }

                // simplified via nA_12 == 0
                float nA_00 = A_00 - t * A_02;
                float nA_22 = A_22 + t * A_02;

                float nA_01 = /* c * A_01 */ -s * A_12;
                float nA_02 = 0.0f; // focus
                float nA_12 = c * A_12 /* + s * A_01*/;
                A_00 = nA_00;
                A_01 = nA_01;
                A_02 = nA_02;
                A_12 = nA_12;
                A_22 = nA_22;


                // Converged when the smallest diag can not get affected from the largest non-diag
                float minDiag = min_of_abs(A_00, A_11, A_22);
                if (minDiag + max_of_abs(A_12, A_01) == minDiag)
                {
                    break;
                }
            }
        }
    }

    glm::vec3 wbasisZ = glm::cross(wbasisXY[0], wbasisXY[1]);
    lambdas[0] = A_00;
    lambdas[1] = A_11;
    lambdas[2] = A_22;
    es[0] = { wbasisXY[0].x, wbasisXY[1].x, wbasisZ.x };
    es[1] = { wbasisXY[0].y, wbasisXY[1].y, wbasisZ.y };
    es[2] = { wbasisXY[0].z, wbasisXY[1].z, wbasisZ.z };
}

glm::vec3 plasma_quintic(float x)
{
    x = glm::clamp(x, 0.0f, 1.0f);
    glm::vec4 x1 = glm::vec4(1.0f, x, x * x, x * x * x); // 1 x x2 x3
    glm::vec4 x2 = x1 * x1.w * x; // x4 x5 x6 x7
    return glm::vec3(
        glm::dot(x1, glm::vec4(+0.063861086f, +1.992659096f, -1.023901152f, -0.490832805f)) + glm::dot(glm::vec2(x2), glm::vec2(+1.308442123f, -0.914547012f)),
        glm::dot(x1, glm::vec4(+0.049718590f, -0.791144343f, +2.892305078f, +0.811726816f)) + glm::dot(glm::vec2(x2), glm::vec2(-4.686502417f, +2.717794514f)),
        glm::dot(x1, glm::vec4(+0.513275779f, +1.580255060f, -5.164414457f, +4.559573646f)) + glm::dot(glm::vec2(x2), glm::vec2(-1.916810682f, +0.570638854f)));
}

class PLYFile
{
public:
    void load( const char* file )
    {
        FILE* fp = fopen(file, "rb");
        PR_ASSERT(fp);

        int elementVertex = 0;

        int head = 0;
        std::map<std::string, int> vAttribOffets;
        
        char line[256];
        bool header_ended = false;
        while (fgets(line, sizeof(line), fp))
        {
            if (strstr(line, "format"))
            {
                PR_ASSERT(strcmp(line, "format binary_little_endian 1.0\n") == 0 );
            }

            if(strstr(line, "element vertex"))
            {
                sscanf(line, "element vertex %d", &elementVertex);
            }
            if (strstr(line, "property"))
            {
                char type[256];
                char field[256];
                sscanf(line, "property %s %s", type, field);
                PR_ASSERT(strcmp(type, "float") == 0);

                vAttribOffets[field] = head;
                head += 4;
                printf("");
            }
            if (strcmp(line, "end_header\n") == 0)
            {
                header_ended = true;
                break;
            }
        }

        PR_ASSERT(header_ended);

        std::vector<uint8_t> binary( head * elementVertex );
        int nRead = fread(binary.data(), head * elementVertex, 1, fp);
        PR_ASSERT(nRead == 1);

        fclose(fp);

        m_vAttribOffets = std::move(vAttribOffets);
        m_size = elementVertex;
        m_binary = std::move(binary);
        m_stride = head;
    }
    int size() const
    {
        return m_size;
    }
    int attrib_offset(const char* field) const
    {
        return m_vAttribOffets.find(field)->second;
    }
    float value( int index, int offset ) const {
        float x;
        const uint8_t* ptr = m_binary.data() + index * m_stride + offset;
        memcpy(&x, ptr, sizeof(float));
        return x;
    }

    int m_size = 0;
    int m_stride = 0;
    std::map<std::string, int> m_vAttribOffets;
    std::vector<uint8_t> m_binary;
};
static float remap(float value, float inputMin, float inputMax, float outputMin, float outputMax)
{
    return (value - inputMin) * ((outputMax - outputMin) / (inputMax - inputMin)) + outputMin;
}
inline int32_t __float_as_int(float x)
{
    int i;
    memcpy(&i, &x, 4);
    return i;
}
inline uint32_t __float_as_uint(float x)
{
    uint32_t i;
    memcpy(&i, &x, 4);
    return i;
}
inline uint32_t getKeyBits( float x )
{
    if (x == 0.0f)
        x = 0.0f;

    uint32_t flip = uint32_t(__float_as_int(x) >> 31) | 0x80000000;
    return __float_as_uint(x) ^ flip;
}

float sigmoid(float x)
{
    return 1.0f / (1.0f + std::exp(-x));
}

glm::mat2 applyJacobian( float cov_00, float cov_01, float cov_02, float cov_11, float cov_12, float cov_22, float ka, float kb, float kc )
{
    float ka2 = ka * ka;
    float kb2 = kb * kb;
    float kc2 = kc * kc;

    float m00 = cov_00 * ka2 + cov_02 * 2.0f * ka * kb + cov_22 * kb2;
    float m11 = cov_11 * ka2 + cov_12 * 2.0f * ka * kc + cov_22 * kc2;
    float m01 = cov_01 * ka2 + cov_12 * ka * kb + cov_02 * ka * kc + cov_22 * kb * kc;
    return { m00, m01, m01, m11 };
}

float exp_approx(float x)
{
    // return expf( glm::clamp( x, -16.0f, 16.0f ) ); // use this for numerical varidation

    /*
    float L = 0.0f;
    float R = 1.0f;
    for (int i = 0 ; i < 1000 ; i++)
    {
        float m = ( L + R ) * 0.5f;
        float x = m;
        x *= x;
        x *= x;
        x *= x;
        if( x == 0.0f || fpclassify(x) == FP_SUBNORMAL )
        {
            L = m;
        }
        else
        {
            R = m;
        }
    }
    printf( "%.32f\n", R ); >> 0.00001814586175896693021059036255
    */
    x = 1.0f + x / 8.0f;
    if (x < 0.00001814586175896693021059036255f) // avoid subnormal
    {
        return 0.0f;
    }
    x *= x;
    x *= x;
    x *= x;
    return x;
}

int main() {
    using namespace pr;

    SetDataDir(ExecutableDir());

    PLYFile pointCould;
    pointCould.load(GetDataPath("models/bicycle/point_cloud/iteration_30000/point_cloud.ply").c_str());
    int ATTRIB_X = pointCould.attrib_offset("x");
    int ATTRIB_Y = pointCould.attrib_offset("y");
    int ATTRIB_Z = pointCould.attrib_offset("z");
    int ATTRIB_R = pointCould.attrib_offset("f_dc_0");
    int ATTRIB_G = pointCould.attrib_offset("f_dc_1");
    int ATTRIB_B = pointCould.attrib_offset("f_dc_2");
    int ATTRIB_QX = pointCould.attrib_offset("rot_0");
    int ATTRIB_QY = pointCould.attrib_offset("rot_1");
    int ATTRIB_QZ = pointCould.attrib_offset("rot_2");
    int ATTRIB_QW = pointCould.attrib_offset("rot_3");
    int ATTRIB_SX = pointCould.attrib_offset("scale_0");
    int ATTRIB_SY = pointCould.attrib_offset("scale_1");
    int ATTRIB_SZ = pointCould.attrib_offset("scale_2");
    int ATTRIB_O = pointCould.attrib_offset("opacity");

    Config config;
    config.ScreenWidth = 1920;
    config.ScreenHeight = 1080;
    config.SwapInterval = 1;
    Initialize(config);

    Camera3D camera;
    camera.origin = { 4, 4, 4 };
    // camera.origin = { 4.41570, 0.82581, -5.27443 };
    camera.lookat = { 0, 0, 0 };
    camera.zNear = 0.4f;
    camera.zFar = 1000.0f;
    // camera.fovy = glm::radians( 90.0f );

    //camera.origin *= 100.0f;
    //camera.fovy = 0.005f;

    double e = GetElapsedTime();

    ITexture* tex = CreateTexture();
    Image2DRGBA32 image;
    int stride = 4;

    while (pr::NextFrame() == false) {
        if (IsImGuiUsingMouse() == false) {
            UpdateCameraBlenderLike(&camera);
        }

        // ClearBackground(0.1f, 0.1f, 0.1f, 1);
        ClearBackground(tex);

        BeginCamera(camera);

        PushGraphicState();

        DrawGrid(GridAxis::XZ, 1.0f, 10, { 128, 128, 128 });

        static glm::vec3 rotation_handle = { -0.18264, -1.49141, -1.84652 };
        ManipulatePosition(camera, &rotation_handle, 0.2f);
        DrawArrow({ 0,0,0 }, rotation_handle, 0.01f, { 0,255,255 });
        glm::mat3 modelMat = glm::identity<glm::mat3>();
        {
            glm::vec3 xaxis, yaxis;
            GetOrthonormalBasis(glm::normalize(rotation_handle), &xaxis, &yaxis);
            modelMat = { xaxis, yaxis, glm::normalize(rotation_handle) };
        }

        // printf("%.5f, %.5f, %.5f\n", rotation_handle.x, rotation_handle.y, rotation_handle.z);
        printf("%.5f, %.5f, %.5f\n", camera.origin.x, camera.origin.y, camera.origin.z);

        DrawSphere({ 4,4,4 }, 0.01f, { 0, 255, 0 });

        SetObjectTransform(modelMat);

        //PrimBegin(PrimitiveMode::Points, 2);
        //for (int i = 0; i < pointCould.size(); i++)
        //{
        //    //if (i != 1783434)
        //    //    continue;

        //    glm::vec3 p = {
        //        pointCould.value(i, ATTRIB_X),
        //        pointCould.value(i, ATTRIB_Y),
        //        pointCould.value(i, ATTRIB_Z),
        //    };
        //    float SH_C0 = 0.28209479177387814;
        //    glm::vec3 splat_col = {
        //        0.5f + SH_C0 * pointCould.value(i, ATTRIB_R),
        //        0.5f + SH_C0 * pointCould.value(i, ATTRIB_G),
        //        0.5f + SH_C0 * pointCould.value(i, ATTRIB_B),
        //    };
        //    splat_col = glm::clamp(splat_col, glm::vec3(0), glm::vec3(1));

        //    PrimVertex(p, splat_col * 255.0f );
        //}
        //PrimEnd();

        //PrimBegin(PrimitiveMode::Lines, 2);
        //for (int i = 0; i < pointCould.size(); i++)
        //{
        //    if (i != 1783434)
        //        continue;
        //    glm::vec3 p = {
        //        pointCould.value(i, ATTRIB_X),
        //        pointCould.value(i, ATTRIB_Y),
        //        pointCould.value(i, ATTRIB_Z),
        //    };
        //    glm::vec3 s = {
        //        std::expf( pointCould.value(i, ATTRIB_SX) ),
        //        std::expf( pointCould.value(i, ATTRIB_SY) ),
        //        std::expf( pointCould.value(i, ATTRIB_SZ) ),
        //    };
        //    glm::quat q = {
        //        pointCould.value(i, ATTRIB_QX),
        //        pointCould.value(i, ATTRIB_QY),
        //        pointCould.value(i, ATTRIB_QZ),
        //        pointCould.value(i, ATTRIB_QW)
        //    };
        //    glm::mat3 R = glm::mat3_cast(q);
        //    glm::u8vec3 CR = { 255, 0, 0 };
        //    glm::u8vec3 CG = { 0, 255, 0 };
        //    glm::u8vec3 CB = { 0, 0, 255 };
        //    PrimVertex(p, CR);
        //    PrimVertex(p + R[0] * s.x, CR);
        //    PrimVertex(p, CG);
        //    PrimVertex(p + R[1] * s.y, CG);
        //    PrimVertex(p, CB);
        //    PrimVertex(p + R[2] * s.z, CB);
        //}
        //PrimEnd();

        //PrimBegin(PrimitiveMode::Lines, 1);
        //for (int i = 0; i < pointCould.size(); i++ )
        //{
        //    glm::vec3 p = {
        //        pointCould.value(i, ATTRIB_X),
        //        pointCould.value(i, ATTRIB_Y),
        //        pointCould.value(i, ATTRIB_Z),
        //    };
        //    glm::vec3 s = {
        //        std::expf(pointCould.value(i, ATTRIB_SX)),
        //        std::expf(pointCould.value(i, ATTRIB_SY)),
        //        std::expf(pointCould.value(i, ATTRIB_SZ)),
        //    };
        //    glm::quat q = {
        //        pointCould.value(i, ATTRIB_QX),
        //        pointCould.value(i, ATTRIB_QY),
        //        pointCould.value(i, ATTRIB_QZ),
        //        pointCould.value(i, ATTRIB_QW)
        //    };
        //    glm::mat3 R = glm::mat3_cast(q);

        //    glm::vec3 u = R[0] * s.x;
        //    glm::vec3 v = R[1] * s.y;
        //    glm::vec3 w = R[2] * s.z;
        //    
        //    glm::u8vec3 color = { 255, 255, 255 };
        //    int vertexCount = 8;
        //    {
        //        CircleGenerator circular(glm::pi<float>() * 2.0f / vertexCount);
        //        for (int i = 0; i <= vertexCount; i++)
        //        {
        //            PrimVertex(p + u * circular.sin() + v * circular.cos(), color);
        //            circular.step();
        //            PrimVertex(p + u * circular.sin() + v * circular.cos(), color);
        //        }
        //    }
        //    {
        //        CircleGenerator circular(glm::pi<float>() * 2.0f / vertexCount);
        //        for (int i = 0; i <= vertexCount; i++)
        //        {
        //            PrimVertex(p + v * circular.sin() + w * circular.cos(), color);
        //            circular.step();
        //            PrimVertex(p + v * circular.sin() + w * circular.cos(), color);
        //        }
        //    }
        //}
        //PrimEnd();
        SetObjectIdentify();

        static bool RT_mode = false;
        static bool showWire = true;


        // https://tavianator.com/2014/ellipsoid_bounding_boxes.html
        //float dx = std::sqrt(U.x * U.x + V.x * V.x + W.x * W.x);
        //float dy = std::sqrt(U.y * U.y + V.y * V.y + W.y * W.y);
        //float dz = std::sqrt(U.z * U.z + V.z * V.z + W.z * W.z);
        //DrawCube({ 0, 0, 0 }, { dx * 2, dy * 2, dz * 2 }, { 255,255,255 });

        //{
        //    float cov_00 = U.x * U.x + V.x * V.x + W.x * W.x;
        //    float cov_01 = U.x * U.y + V.x * V.y + W.x * W.y;
        //    float cov_02 = U.x * U.z + V.x * V.z + W.x * W.z;
        //    float cov_11 = U.y * U.y + V.y * V.y + W.y * W.y;
        //    float cov_12 = U.y * U.z + V.y * V.z + W.y * W.z;
        //    float cov_22 = U.z * U.z + V.z * V.z + W.z * W.z;
        //    glm::mat3 cov = {
        //        cov_00, cov_01, cov_02,
        //        cov_01, cov_11, cov_12,
        //        cov_02, cov_12, cov_22
        //    };

        //    glm::mat3 ellipsoidAffine = { U, V, W };

        //    float lambdas[3];
        //    glm::vec3 es[3];
        //    eigen_decomposition(es, lambdas, cov[0][0], cov[1][0], cov[2][0], cov[1][1], cov[2][1], cov[2][2]);
        //    DrawArrow({ 0,0,0 }, es[0], 0.01f, { 255,0,255 });
        //    DrawArrow({ 0,0,0 }, es[1], 0.01f, { 255,255,0 });
        //    DrawArrow({ 0,0,0 }, es[2], 0.01f, { 0,255,255 });
        //    DrawLine(-es[0] * std::sqrt(lambdas[0]), es[0] * std::sqrt(lambdas[0]), { 255,0,255 });
        //    DrawLine(-es[1] * std::sqrt(lambdas[1]), es[1] * std::sqrt(lambdas[1]), { 255,255,0 });
        //    DrawLine(-es[2] * std::sqrt(lambdas[2]), es[2] * std::sqrt(lambdas[2]), { 0,255,255 });

        //    if (showWire)
        //    {
        //        //SetObjectTransform(ellipsoidAffine);
        //        //DrawSphere({ 0,0,0 }, 1.0f, { 128 ,128 ,128 }, 32, 32);
        //        //SetObjectIdentify();
        //        glm::vec3 e0 = es[0] * std::sqrt(lambdas[0]);
        //        glm::vec3 e1 = es[1] * std::sqrt(lambdas[1]);
        //        glm::vec3 e2 = es[2] * std::sqrt(lambdas[2]);
        //        DrawEllipse({ 0,0,0 }, e0, e1, { 255,255,255 }, 64);
        //        DrawEllipse({ 0,0,0 }, e1, e2, { 255,255,255 }, 64);
        //        DrawEllipse({ 0,0,0 }, e2, e0, { 255,255,255 }, 64);
        //    }
        //}

        image.allocate(GetScreenWidth() / stride, GetScreenHeight() / stride);

        //glm::mat4 viewMat = GetCurrentViewMatrix();
        //glm::mat3 viewRot = glm::mat3(viewMat);

        glm::mat4 viewMat = GetCurrentViewMatrix() * glm::mat4(modelMat);
        glm::mat3 viewRot = glm::mat3(viewMat); 

        static std::vector<int>   pointIndices;
        static std::vector<float> pointDepthes;
        pointIndices.resize(pointCould.size());
        pointDepthes.resize(pointCould.size());
        for (int i = 0; i < pointIndices.size(); i++)
        {
            pointIndices[i] = i;

            glm::vec3 p = {
                pointCould.value(i, ATTRIB_X),
                pointCould.value(i, ATTRIB_Y),
                pointCould.value(i, ATTRIB_Z),
            };
            pointDepthes[i] = (viewMat * glm::vec4( p.x, p.y, p.z, 1 )).z;
        }

        // front to back
#if NDEBUG
        std::sort(pointIndices.begin(), pointIndices.end(), [](int a, int b) { return pointDepthes[a] > pointDepthes[b]; });
#else
        concurrency::parallel_radixsort(pointIndices.begin(), pointIndices.end(), [](int x) { return getKeyBits(pointDepthes[x]) ^ 0xFFFFFFFF; });
#endif

        float tanThetaY = std::tan(camera.fovy * 0.5f);
        float tanThetaX = tanThetaY / image.height() * image.width();

        std::fill(image.data(), image.data() + image.width() * image.height(), glm::vec4(0.0f, 0.0f, 0.0f, 1.0f));

        // for (int i = 0; i < pointCould.size(); i++)
        for( auto i : pointIndices )
        {
            //if (i != 1783434)
            //    continue;

            //if (i == 1783434)
            //    continue;

            glm::vec3 p = {
                pointCould.value(i, ATTRIB_X),
                pointCould.value(i, ATTRIB_Y),
                pointCould.value(i, ATTRIB_Z),
            };
            glm::vec3 s = {
                std::expf(pointCould.value(i, ATTRIB_SX)),
                std::expf(pointCould.value(i, ATTRIB_SY)),
                std::expf(pointCould.value(i, ATTRIB_SZ)),
            };
            glm::quat q = {
                pointCould.value(i, ATTRIB_QX),
                pointCould.value(i, ATTRIB_QY),
                pointCould.value(i, ATTRIB_QZ),
                pointCould.value(i, ATTRIB_QW)
            };
            glm::mat3 R = glm::mat3_cast(q);

            float SH_C0 = 0.28209479177387814;
            glm::vec3 splat_col = {
                0.5f + SH_C0 * pointCould.value(i, ATTRIB_R),
                0.5f + SH_C0 * pointCould.value(i, ATTRIB_G),
                0.5f + SH_C0 * pointCould.value(i, ATTRIB_B),
            };
            splat_col = glm::clamp(splat_col, glm::vec3(0), glm::vec3(1));
            float opacity = sigmoid ( pointCould.value(i, ATTRIB_O) );

            glm::vec3 u = viewRot * R[0] * s.x;
            glm::vec3 v = viewRot * R[1] * s.y;
            glm::vec3 w = viewRot * R[2] * s.z;

            float cov_00 = u.x * u.x + v.x * v.x + w.x * w.x;
            float cov_01 = u.x * u.y + v.x * v.y + w.x * w.y;
            float cov_02 = u.x * u.z + v.x * v.z + w.x * w.z;
            float cov_11 = u.y * u.y + v.y * v.y + w.y * w.y;
            float cov_12 = u.y * u.z + v.y * v.z + w.y * w.z;
            float cov_22 = u.z * u.z + v.z * v.z + w.z * w.z;

            glm::mat3 cov = {
                cov_00, cov_01, cov_02,
                cov_01, cov_11, cov_12,
                cov_02, cov_12, cov_22
            };
            glm::vec3 u_camera = viewMat * glm::vec4(p.x, p.y, p.z, 1);
            //glm::vec3 u_camera = viewMat * glm::vec4(0,0,0, 1);

            // culling back
            if( -0.5f < u_camera.z )
            {
                continue;
            }
            glm::vec2 x_rayspace = glm::vec2(u_camera.x / -u_camera.z, u_camera.y / -u_camera.z);

            // Apply a rough culling. 
            // This can also avoid unstable splat projection
            // aabb in viewspace
            // float dx2 = (u.x * u.x + v.x * v.x + w.x * w.x);
            // float dy2 = (u.y * u.y + v.y * v.y + w.y * w.y);
            // float dz2 = (u.z * u.z + v.z * v.z + w.z * w.z);
            float diagH_raySpace = sqrt( cov_00 + cov_11 + cov_22 ) / -u_camera.z; // half diagonal
            if(x_rayspace.x + diagH_raySpace < -tanThetaX || tanThetaX < x_rayspace.x - diagH_raySpace ||
                x_rayspace.y + diagH_raySpace < -tanThetaY || tanThetaY < x_rayspace.y - diagH_raySpace )
            {
                continue;
            }

            //glm::mat3 J;
            //J = {
            //    1.0f / u_camera.z, 0, 0,
            //    0, 1.0f / u_camera.z, 0,
            //    -u_camera.x / (u_camera.z * u_camera.z), -u_camera.y / (u_camera.z * u_camera.z), 1.0f
            //};

            //glm::mat3 covPrime = J * cov * glm::transpose(J);
            //glm::mat2 covPrime2d = glm::mat2(
            //    covPrime[0][0], covPrime[0][1],
            //    covPrime[1][0], covPrime[1][1]
            //);
            float pxSize = tanThetaX * 2.0f / image.width();
            glm::mat2 covPrime2d = applyJacobian(cov[0][0], cov[1][0], cov[2][0], cov[1][1], cov[2][1], cov[2][2], 1.0f / u_camera.z, -u_camera.x / (u_camera.z * u_camera.z), -u_camera.y / (u_camera.z * u_camera.z));
            float det_of_cov = glm::determinant(covPrime2d);
            if (glm::abs(det_of_cov) < sqr(pxSize / 32.0f * pxSize / 32.0f))
            {
                continue;
            }
            
            // dilation approach 
            //covPrime2d[0][0] += sqr( pxSize / 16.0f );
            //covPrime2d[1][1] += sqr( pxSize / 16.0f );

            glm::mat2 invCovPrime2d = glm::inverse(covPrime2d);



            //if( glm::abs(det_of_cov) < 0.0000001f )
            //    continue;

            // The exact bounding box from covariance matrix
            float hsize_invCovY = std::sqrt(invCovPrime2d[0][0] * det_of_cov) * SPLAT_BOUNDS;
            //float lowerY = x_rayspace.y + hsize_invCovY;
            //float upperY = x_rayspace.y - hsize_invCovY;
            //int begY = remap( lowerY, tanThetaY, -tanThetaY, 0, image.height() - 1 );
            //int endY = remap( upperY, tanThetaY, -tanThetaY, 0, image.height() - 1 );

            //int midY = remap(x_rayspace.y, tanThetaY, -tanThetaY, 0, image.height() - 1);
            //int begY = midY - 2;
            //int endY = midY + 2;
            int begY = remap(x_rayspace.y + hsize_invCovY, tanThetaY, -tanThetaY, 0, image.height() - 1);
            int endY = remap(x_rayspace.y - hsize_invCovY, tanThetaY, -tanThetaY, 0, image.height() - 1);
            if (endY < 0 || image.height() <= begY)
            {
                continue;
            }
            begY = ss_max(begY, 0);
            endY = ss_min(endY, image.height() - 1);

            for (int y = begY; y <= endY; y++)
            {
                if( y < 0 || image.height() <= y )
                    continue;

                // Minimum range of x
                float vy = remap( y, 0, image.height() - 1, tanThetaY, -tanThetaY ) - x_rayspace.y;
                float a = invCovPrime2d[0][0];
                float b = invCovPrime2d[1][0];
                float d = invCovPrime2d[1][1];
                float xs[2];
                int begX = -1;
                int endX = -1;
                if (solve_quadratic(xs, a, 2.0f * b * vy, d * vy * vy - SPLAT_BOUNDS * SPLAT_BOUNDS))
                {
                    //float lowerX = x_rayspace.x + xs[0];
                    //float upperX = x_rayspace.x + xs[1];
                    //begX = (lowerX + tanThetaX) * pxPerTanTheta;
                    //endX = (upperX + tanThetaX) * pxPerTanTheta;
                    begX = remap(x_rayspace.x + xs[0], -tanThetaX, tanThetaX, 0, image.width() - 1);
                    endX = remap(x_rayspace.x + xs[1], -tanThetaX, tanThetaX, 0, image.width() - 1);
                }
                if (endX < 0 || image.width() <= begX)
                {
                    continue;
                }
                begX = ss_max(begX, 0);
                endX = ss_min(endX, image.width() - 1);

                //if (endY - begY > 500)
                //{
                //    printf("");
                //}

                //float midX = remap( x_rayspace.x, -tanThetaX, tanThetaX, 0, image.width() - 1 );
                //int begX = midX - 2;
                //int endX = midX + 2;

                for (int x = begX; x <= endX; x++)
                {
                    if (x < 0 || image.width() <= x)
                        continue;

                    //image(x, y) = glm::vec4(1, 1, 1, 1);
                    // image(x, y) = glm::vec4(splat_col, 1);

                    // w as throughput
                    glm::vec4 color = image(x, y);
                    float T = color.w;

                    if (T < MIN_THROUGHPUT)
                        continue;

                    glm::vec2 p_rayspace = {
                        remap(x, 0, image.width() - 1, -tanThetaX, tanThetaX ),
                        remap(y, 0, image.height() - 1, tanThetaY, -tanThetaY )
                    };
                    glm::vec2 v = p_rayspace - x_rayspace;

                    float d2 = glm::dot(v, invCovPrime2d * v);
                    float alpha = exp_approx(-0.5f * ss_max(d2, 0.0f)) * opacity;

                    color.x += T * splat_col.x * alpha;
                    color.y += T * splat_col.y * alpha;
                    color.z += T * splat_col.z * alpha;

                    color.w *= (1.0f - alpha);

                    //if (isfinite(color.x) == false)
                    //{
                    //    printf("a");
                    //    exp_approx(-0.5f * d2);
                    //}
                    image(x, y) = color;
                    // image(x, y) = { 1, 0, 1, color.w };
                }
            }
        }

        {
            glm::vec4* ptr = image.data();
            for (int i = 0; i < image.width() * image.height(); i++)
            {
                //ptr[i].x = std::pow(  ptr[i].x, 1.0f / 2.2f );
                //ptr[i].y = std::pow(  ptr[i].y, 1.0f / 2.2f );
                //ptr[i].z = std::pow(  ptr[i].z, 1.0f / 2.2f );
                ptr[i].w = 1.0f;
            }
        }

        //CameraRayGenerator rayGenerator(GetCurrentViewMatrix(), GetCurrentProjMatrix(), image.width(), image.height());

        //glm::mat4 viewMat = GetCurrentViewMatrix();
        //glm::mat3 viewRot = glm::mat3(viewMat);

        //for (int j = 0; j < image.height(); ++j)
        //{
        //    for (int i = 0; i < image.width(); ++i)
        //    {
        //        if (RT_mode)
        //        {
        //            glm::vec3 ro, rd;
        //            rayGenerator.shoot(&ro, &rd, i, j, 0.5f, 0.5f);

        //            glm::vec3 Z = glm::normalize(rd);
        //            glm::vec3 X;
        //            glm::vec3 Y;
        //            GetOrthonormalBasis(Z, &X, &Y);

        //            glm::mat3 toLocal = {
        //                glm::vec3(X.x, Y.x, Z.x),
        //                glm::vec3(X.y, Y.y, Z.y),
        //                glm::vec3(X.z, Y.z, Z.z),
        //            };
        //            glm::vec3 u = toLocal * U;
        //            glm::vec3 v = toLocal * V;
        //            glm::vec3 w = toLocal * W;

        //            glm::vec2 v_rel = toLocal * ro;

        //            float cov_00 = u.x * u.x + v.x * v.x + w.x * w.x;
        //            float cov_01 = u.x * u.y + v.x * v.y + w.x * w.y;
        //            float cov_11 = u.y * u.y + v.y * v.y + w.y * w.y;
        //            glm::mat2 cov = { cov_00, cov_01, cov_01, cov_11 };
        //            glm::mat2 inv_cov = glm::inverse(cov);

        //            float d2 = glm::dot(v_rel, inv_cov * v_rel);
        //            float alpha = glm::exp(-0.5f * d2);
        //            glm::u8vec3 color = glm::u8vec3(glm::clamp(plasma_quintic(alpha) * 255.0f, glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(255.0f, 255.0f, 255.0f)));

        //            if (glm::abs(d2 - 1.0f) < 0.005f)
        //            {
        //                color = { 255,255,255 };
        //            }
        //            if (glm::abs(d2 - 4.0f) < 0.005f)
        //            {
        //                color = { 128,128,128 };
        //            }
        //            image(i, j) = glm::u8vec4(color, 255);
        //        }
        //        else
        //        {
        //            // viewRot
        //            glm::vec3 u = viewRot * U;
        //            glm::vec3 v = viewRot * V;
        //            glm::vec3 w = viewRot * W;

        //            float cov_00 = u.x * u.x + v.x * v.x + w.x * w.x;
        //            float cov_01 = u.x * u.y + v.x * v.y + w.x * w.y;
        //            float cov_02 = u.x * u.z + v.x * v.z + w.x * w.z;
        //            float cov_11 = u.y * u.y + v.y * v.y + w.y * w.y;
        //            float cov_12 = u.y * u.z + v.y * v.z + w.y * w.z;
        //            float cov_22 = u.z * u.z + v.z * v.z + w.z * w.z;
        //            glm::mat3 cov = {
        //                cov_00, cov_01, cov_02,
        //                cov_01, cov_11, cov_12,
        //                cov_02, cov_12, cov_22
        //            };
        //            // glm::vec3 u_ray = viewMat * glm::vec4( ( /* gaussian center */ - camera.origin ), 1.0f );
        //            glm::vec3 u_ray = viewMat * glm::vec4(0, 0, 0, 1);

        //            glm::mat3 J;
        //            J = {
        //                1.0f / u_ray.z, 0, 0,
        //                0, 1.0f / u_ray.z, 0,
        //                -u_ray.x / (u_ray.z * u_ray.z), -u_ray.y / (u_ray.z * u_ray.z), 1.0f
        //            };

        //            glm::mat3 covPrime = J * cov * glm::transpose(J);
        //            glm::mat2 covPrime2d = glm::mat2(
        //                covPrime[0][0], covPrime[0][1],
        //                covPrime[1][0], covPrime[1][1]
        //            );

        //            glm::mat2 inv_cov = glm::inverse(covPrime2d);

        //            float tanThetaY = std::tan(camera.fovy * 0.5f);
        //            float tanThetaX = tanThetaY / image.height() * image.width();

        //            glm::vec2 v_rel = glm::vec2(
        //                glm::mix(-tanThetaX, tanThetaX, (float)i / image.width()),
        //                glm::mix(tanThetaY, -tanThetaY, (float)j / image.height())
        //            ) - glm::vec2(u_ray.x / -u_ray.z, u_ray.y / -u_ray.z); // z is negative

        //            float d2 = glm::dot(v_rel, inv_cov * v_rel);
        //            float alpha = std::expf(-0.5f * d2);
        //            glm::u8vec3 color = glm::u8vec3(glm::clamp(plasma_quintic(alpha) * 255.0f, glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(255.0f, 255.0f, 255.0f)));

        //            if (glm::abs(d2 - 1.0f) < 0.005f)
        //            {
        //                color = { 255,255,255 };
        //            }
        //            if (glm::abs(d2 - 4.0f) < 0.005f)
        //            {
        //                color = { 128,128,128 };
        //            }
        //            image(i, j) = glm::u8vec4(color, 255);

        //            // image(i, j) = glm::u8vec4(0, 0, 0, 255);
        //        }
        //    }
        //}
        


        tex->upload(image);

        PopGraphicState();
        EndCamera();

        BeginImGui();

        ImGui::SetNextWindowSize({ 500, 800 }, ImGuiCond_Once);
        ImGui::Begin("Panel");
        ImGui::Text("fps = %f", GetFrameRate());
        ImGui::Text("%d points", pointCould.size());

        ImGui::SliderFloat("fov", &camera.fovy, 0, 0.1);
        ImGui::Text("cam d %f", glm::length(camera.origin));
        ImGui::Checkbox("RT_mode", &RT_mode);
        ImGui::Checkbox("showWire", &showWire);
        ImGui::End();

        EndImGui();
    }

    pr::CleanUp();
}

#pragma once

#include <cmath>
#include <xmmintrin.h>
#include <smmintrin.h>
#include <sstream>
#include <algorithm>

namespace svector {
    class float4 {
    public:
        union {
            struct {
                float w, x, y, z;
            };
            struct {
                float a, r, g, b;
            };
            __m128 xmm;
        };

        static const __m128 mm_one;
        static const __m128 mm_two;

        // constructors
        inline float4() : xmm(_mm_setzero_ps()) {}

        inline float4(const float4& t) : xmm(_mm_set_ps(t.z, t.y, t.x, t.w)) {};

        inline float4(const float4& t, float w) : xmm(_mm_set_ps(t.z, t.y, t.x, w)) {};

        inline float4(float x, float y, float z, float w = 0.0f)
            : xmm(_mm_set_ps(z, y, x, w)) {};

        inline float4(float v) : xmm(_mm_set1_ps(v)) {};

        inline float4(__m128 xmmx) : xmm(xmmx) {}

        // assigments
        inline void operator()(const float4& t) {
            xmm = _mm_set_ps(t.z, t.y, t.x, t.w);
        }

        inline void operator()(float x, float y, float z, float w) {
            xmm = _mm_set_ps(z, y, x, w);
        }

        inline void zero() { xmm = _mm_setzero_ps(); }

        // unary operators
        inline float4 operator+() { return *this; }
        inline float4 operator-() {
            float4 t(-x, -y, -z, -w);
            return t;
        }
        inline float4 operator++() {
            xmm = _mm_add_ps(xmm, mm_one);
            return *this;
        }
        inline float4 operator--() {
            xmm = _mm_sub_ps(xmm, mm_one);
            return *this;
        }
        inline float4 operator++(int d) {
            float4 v(*this);
            xmm = _mm_add_ps(xmm, mm_one);
            return v;
        }
        inline float4 operator--(int d) {
            float4 v(*this);
            xmm = _mm_sub_ps(xmm, mm_one);
            return v;
        }

        // xmm mutators
        inline void operator=(__m128 xmmx) { xmm = xmmx; }

        // scalar mutators
        inline void operator=(float v) { xmm = _mm_set1_ps(v); }
        inline float4& operator+=(float v) {
            float4 t(v);
            xmm = _mm_add_ps(xmm, t.xmm);
            return *this;
        }

        float4& operator-=(float v) {
            float4 t(v);
            xmm = _mm_sub_ps(xmm, t.xmm);
            return *this;
        }

        float4& operator*=(float v) {
            float4 t(v);
            xmm = _mm_mul_ps(xmm, t.xmm);
            return *this;
        }

        float4& operator/=(float v) {
            float4 t(v);
            xmm = _mm_div_ps(xmm, t.xmm);
            return *this;
        }

        // vector mutators
        inline void operator=(const float4& v) {
            xmm = _mm_set_ps(v.z, v.y, v.x, v.w);
        }

        inline float4& operator+=(const float4& v) {
            xmm = _mm_add_ps(xmm, v.xmm);
            return *this;
        }

        float4& operator-=(const float4& v) {
            xmm = _mm_sub_ps(xmm, v.xmm);
            return *this;
        }

        float4& operator*=(const float4& v) {
            xmm = _mm_mul_ps(xmm, v.xmm);
            return *this;
        }

        float4& operator/=(const float4& v) {
            xmm = _mm_div_ps(xmm, v.xmm);
            return *this;
        }

        // math
        inline float norm(const unsigned char mask = 0xF1) {
            //#ifdef _WIN32
            return _mm_cvtss_f32(_mm_sqrt_ss(_mm_dp_ps(xmm, xmm, 0xF1)));
            //#else
            //      return _mm_cvtss_f32(_mm_sqrt_ss(_mm_dp_ps(xmm, xmm, mask)));
            //#endif
        }

        inline float norm2(const unsigned char mask = 0xF1) {
#ifdef _WIN32
            return _mm_cvtss_f32(_mm_dp_ps(xmm, xmm, 0xF1));
#else
            return _mm_cvtss_f32(_mm_dp_ps(xmm, xmm, mask));
#endif
        }

        inline void normalize(const unsigned char mask = 0xF1) {
            float v = norm(mask);
            if (std::fabs(v) > 0.0001) {
                float4 t(v);
                xmm = _mm_div_ps(xmm, t.xmm);
            }
        }

        inline void clip(const float vmin, const float vmax) {
            x = std::min(std::max(x, vmin), vmax);
            y = std::min(std::max(y, vmin), vmax);
            z = std::min(std::max(z, vmin), vmax);
            w = std::min(std::max(w, vmin), vmax);
        }

        // quaternions functions
        inline float4& quaternion_mult(const float4& q) {
            __m128 q0 = _mm_set_ps(q.y, -q.z, q.w, q.x);
            __m128 q1 = _mm_set_ps(-q.x, q.w, q.z, q.y);
            __m128 q2 = _mm_set_ps(q.w, q.x, -q.y, q.z);
            __m128 q3 = _mm_set_ps(-q.z, -q.y, -q.x, q.w);

            float tx = _mm_cvtss_f32(_mm_dp_ps(q0, xmm, 0xf1));
            float ty = _mm_cvtss_f32(_mm_dp_ps(q1, xmm, 0xf1));
            float tz = _mm_cvtss_f32(_mm_dp_ps(q2, xmm, 0xf1));
            float tw = _mm_cvtss_f32(_mm_dp_ps(q3, xmm, 0xf1));

            xmm = _mm_set_ps(tz, ty, tx, tw);
            return *this;
        }

        inline void euler(float xe, float ye, float ze) {
            __m128 mmh = _mm_set_ps(ze, ye, xe, 0);
            mmh = _mm_div_ps(mmh, mm_two);
#ifdef _WIN32
            __declspec(align(16)) float h[4] = { 0 };
#else
            float h[4] = { 0 };
#endif
            _mm_store_ps(h, mmh);

            float c[3] = { float(cos(h[1])), float(cos(h[2])), float(cos(h[3])) };
            float s[3] = { float(sin(h[1])), float(sin(h[2])), float(sin(h[3])) };
            float tw = c[0] * c[1] * c[2] + s[0] * s[2] * s[1];
            float tx = s[0] * c[1] * c[2] - c[0] * s[2] * s[1];
            float ty = c[0] * s[1] * c[2] + s[0] * c[2] * s[1];
            float tz = c[0] * c[1] * s[2] - s[0] * s[2] * c[1];
            xmm = _mm_set_ps(tz, ty, tx, tw);
            normalize();
        }

        inline void rotationMatrix(float* Matrix) {
            float xx = x * x;
            float xy = x * y;
            float xz = x * z;
            float xw = w * x;
            float yy = y * y;
            float yz = y * z;
            float yw = w * y;
            float zz = z * z;
            float zw = w * z;

            Matrix[0] = 1 - 2 * (yy + zz);
            Matrix[1] = 2 * (xy - zw);
            Matrix[2] = 2 * (xz + yw);

            Matrix[4] = 2 * (xy + zw);
            Matrix[5] = 1 - 2 * (xx + zz);
            Matrix[6] = 2 * (yz - xw);

            Matrix[8] = 2 * (xz - yw);
            Matrix[9] = 2 * (yz + xw);
            Matrix[10] = 1 - 2 * (xx + yy);

            Matrix[3] = Matrix[7] = Matrix[11] = Matrix[12] = Matrix[13] = Matrix[14] =
                0;
            Matrix[15] = 1;
        }

        inline float4 axis() {
            float temp_angle = acos(w) * 2;
            float vnorm = norm(0xE1);
            __m128 t;
            if (std::fabs(vnorm) < 0.001) {
                t = _mm_set_ps(0, 0, 1, 0);
            }
            else {
                __m128 n = _mm_set1_ps(vnorm);
                t = _mm_div_ps(xmm, n);
            }
            float4 ret(t);
            ret.w = temp_angle;
            return ret;
        }

        // string
        const std::string str() {
            std::ostringstream output;
            output << "[" << x << ", " << y << ", " << z << ", " << w << "]";
            return output.str();
        }

    };

    const __m128 float4::mm_one = _mm_set1_ps(1.0f);
    const __m128 float4::mm_two = _mm_set1_ps(2.0f);

    // binary operators
    //================= (+) ====================
    inline float4 operator+(float k, const float4& v) {
        float4 w(k);
        w += v;
        return w;
    }

    inline float4 operator+(const float4& v, float k) {
        float4 w(k);
        w += v;
        return w;
    }

    inline float4 operator+(const float4& v, const float4& w) {
        float4 u;
        u.xmm = _mm_add_ps(v.xmm, w.xmm);
        return u;
    }

    //================= (-) ====================
    inline float4 operator-(float k, const float4& v) {
        float4 w(k);
        w -= v;
        return w;
    }

    inline float4 operator-(const float4& v, float k) {
        float4 w(k);
        return float4(_mm_sub_ps(v.xmm, w.xmm));
        ;
    }

    inline float4 operator-(const float4& v, const float4& w) {
        return float4(_mm_sub_ps(v.xmm, w.xmm));
    }

    //================= (*) ====================
    inline float4 operator*(float k, const float4& v) {
        float4 w(k);
        w *= v;
        return w;
    }

    inline float4 operator*(const float4& v, float k) {
        float4 w(k);
        return float4(_mm_mul_ps(v.xmm, w.xmm));
        ;
    }

    inline float4 operator*(const float4& v, const float4& w) {
        return float4(_mm_mul_ps(v.xmm, w.xmm));
    }

    //================= (/) ====================
    inline float4 operator/(float k, const float4& v) {
        float4 w(k);
        w /= v;
        return w;
    }

    inline float4 operator/(const float4& v, float k) {
        float4 w(k);
        return float4(_mm_div_ps(v.xmm, w.xmm));
        ;
    }

    inline float4 operator/(const float4& v, const float4& w) {
        return float4(_mm_div_ps(v.xmm, w.xmm));
    }


    // comparison Operators
    inline float4 operator>=(const float4& v, const float4& w) {
        return float4(_mm_cmpge_ps(v.xmm, w.xmm));
    }

    inline float4 operator<=(const float4& v, const float4& w) {
        return float4(_mm_cmple_ps(v.xmm, w.xmm));
    }

    inline float4 operator>(const float4& v, const float4& w) {
        return float4(_mm_cmpgt_ps(v.xmm, w.xmm));
    }

    inline float4 operator<(const float4& v, const float4& w) {
        return float4(_mm_cmplt_ps(v.xmm, w.xmm));
    }

    inline float4 operator==(const float4& v, const float4& w) {
        return float4(_mm_cmpeq_ps(v.xmm, w.xmm));
    }


    // other operators
    inline float dot3d(const float4& v, const float4& w) {
        return _mm_cvtss_f32(_mm_dp_ps(v.xmm, w.xmm, 0xE1));
    }

    inline float dot(const float4& v, const float4& w,
        const unsigned char mask = 0xF1) {
#ifdef _WIN32
        return _mm_cvtss_f32(_mm_dp_ps(v.xmm, w.xmm, 0xF1));
#else
        return _mm_cvtss_f32(_mm_dp_ps(v.xmm, w.xmm, 0xF1));
        //return _mm_cvtss_f32(_mm_dp_ps(v.xmm, w.xmm, mask));
#endif
    }

    inline float norm(const float4& v) {
        return _mm_cvtss_f32(_mm_sqrt_ss(_mm_dp_ps(v.xmm, v.xmm, 0xF1)));
    }

    inline float4 cross3d(const float4& v, const float4& w) {
        float4 x(v.y * w.z - v.z * w.y, v.z * w.x - v.x * w.z, v.x * w.y - v.y * w.x);
        return x;
    }

    inline float4 normal3d(const float4& v, const float4& w) {
        float4 x = cross3d(v, w);
        x.normalize();
        return x;
    }

    inline float4 sqrt(const float4& v) { return _mm_sqrt_ps(v.xmm); }

    inline float4 sqr(const float4& v) { return _mm_mul_ps(v.xmm, v.xmm); }

    inline float4 max(const float4& v1, const float4& v2) {
        float x = std::max(v1.x, v2.x);
        float y = std::max(v1.y, v2.y);
        float z = std::max(v1.z, v2.z);
        float w = std::max(v1.w, v2.w);
        return float4(x, y, z, w);
    }

    inline float4 min(const float4& v1, const float4& v2) {
        float x = std::min(v1.x, v2.x);
        float y = std::min(v1.y, v2.y);
        float z = std::min(v1.z, v2.z);
        float w = std::min(v1.w, v2.w);
        return float4(x, y, z, w);
    }

    inline float4 max(const float4& v1, const float4& v2, const float4& v3) {
        return max(max(v1, v2), v3);
    }

    inline float4 max(const float4& v1, const float4& v2, const float4& v3, const float4& v4) {
        return max(max(v1, v2, v3), v4);
    }

    inline float4 max(const float4& v1, const float4& v2, const float4& v3, const float4& v4, const float4& v5) {
        return max(max(v1, v2, v3, v4), v5);
    }

    inline float4 min(const float4& v1, const float4& v2, const float4& v3) {
        return min(min(v1, v2), v3);
    }

    inline float4 min(const float4& v1, const float4& v2, const float4& v3, const float4& v4) {
        return min(min(v1, v2, v3), v4);
    }

    inline float4 min(const float4& v1, const float4& v2, const float4& v3, const float4& v4, const float4& v5) {
        return min(min(v1, v2, v3, v4), v5);
    }

    inline float4 clamp(const float4& v, const float a, const float b) {
        float x = std::min(std::max(v.x, a), b);
        float y = std::min(std::max(v.y, a), b);
        float z = std::min(std::max(v.z, a), b);
        float w = std::min(std::max(v.w, a), b);

        return float4(x, y, z, w);
    }

    inline float4 saturate(const float4& v) {
        return clamp(v, 0.0f, 1.0f);
    }


    inline float4 fabs(const float4& v) {
        float x = std::fabs(v.x);
        float y = std::fabs(v.y);
        float z = std::fabs(v.z);
        float w = std::fabs(v.w);
        return float4(x, y, z, w);
    }

    inline float4 quaternion_mult(const float4& v, const float4& w) {
        float4 x(v);
        x.quaternion_mult(w);
        return x;
    }

}  // namespace svector
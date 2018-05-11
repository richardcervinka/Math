#pragma once

// Implementation switch.
//#define COMPILE_MATH_MATRIX_DIRECTX_MATH

// DirectX Math implementation helpers.
#ifdef COMPILE_MATH_MATRIX_DIRECTX_MATH

#include "Libs\DirectXMath\Inc\DirectXMath.h"

inline DirectX::XMMATRIX LoadXmMatrix(const Math::Matrix& m)
{
    return DirectX::XMLoadFloat4x4A(reinterpret_cast<const DirectX::XMFLOAT4X4A*>(&m));
}

inline void StoreXmMatrix(const DirectX::XMMATRIX& m, Math::Matrix& o)
{
    DirectX::XMStoreFloat4x4A(reinterpret_cast<DirectX::XMFLOAT4X4A*>(&o), m);
}

inline DirectX::XMVECTOR LoadXmVector(const Math::Vector& v)
{
    return DirectX::XMLoadFloat4A(reinterpret_cast<const DirectX::XMFLOAT4A*>(&v));
}

inline Math::Vector CreateVector(const DirectX::XMVECTOR& v)
{
    Math::Vector result;
    DirectX::XMStoreFloat4A(reinterpret_cast<DirectX::XMFLOAT4A*>(&result), v);
    return result;
}

#endif // COMPILE_MATH_MATRIX_DIRECTX_MATH

namespace Math
{

    inline float Matrix::Determinant(const Matrix& m)
    {
        #ifdef COMPILE_MATH_MATRIX_DIRECTX_MATH

        const DirectX::XMMATRIX xm = LoadXmMatrix(m);
        const DirectX::XMVECTOR xv = DirectX::XMMatrixDeterminant(xm);
        return DirectX::XMVectorGetX(xv);

        #else

        const float t0 = m.m[2][2] * m.m[3][3] - m.m[2][3] * m.m[3][2];
        const float t1 = m.m[2][3] * m.m[3][1] - m.m[2][1] * m.m[3][3];
        const float t2 = m.m[2][1] * m.m[3][2] - m.m[2][2] * m.m[3][1];

        float det00 = {
            m.m[1][1] * t0 +
            m.m[1][2] * t1 +
            m.m[1][3] * t2
        };

        float det10 = {
            m.m[0][1] * t0 +
            m.m[0][2] * t1 +
            m.m[0][3] * t2
        };

        float det20 = {
            m.m[0][1] * (m.m[1][2] * m.m[3][3] - m.m[1][3] * m.m[3][2]) +
            m.m[0][2] * (m.m[1][3] * m.m[3][1] - m.m[1][1] * m.m[3][3]) +
            m.m[0][3] * (m.m[1][1] * m.m[3][2] - m.m[1][2] * m.m[3][1])
        };

        float det30 = {
            m.m[0][1] * (m.m[1][2] * m.m[2][3] - m.m[1][3] * m.m[2][2]) +
            m.m[0][2] * (m.m[1][3] * m.m[2][1] - m.m[1][1] * m.m[2][3]) +
            m.m[0][3] * (m.m[1][1] * m.m[2][2] - m.m[1][2] * m.m[2][1])
        };

        return m.m[0][0] * det00 - m.m[1][0] * det10 + m.m[2][0] * det20 - m.m[3][0] * det30;

        #endif
    }

    inline void Matrix::MatrixTranspose(Matrix& m)
    {
        #ifdef COMPILE_MATH_MATRIX_DIRECTX_MATH

        DirectX::XMMATRIX xm = LoadXmMatrix(m);
        xm = DirectX::XMMatrixTranspose(xm);
        StoreXmMatrix(xm, m);

        #else

        std::swap(m.m[1][0], m.m[0][1]);
        std::swap(m.m[2][0], m.m[0][2]);
        std::swap(m.m[3][0], m.m[0][3]);
        std::swap(m.m[2][1], m.m[1][2]);
        std::swap(m.m[3][1], m.m[1][3]);
        std::swap(m.m[3][2], m.m[2][3]);

        #endif
    }

    inline void Matrix::Transpose(const Matrix& m, Matrix& o)
    {
        #ifdef COMPILE_MATH_MATRIX_DIRECTX_MATH

        DirectX::XMMATRIX xm = LoadXmMatrix(m);
        xm = DirectX::XMMatrixTranspose(xm);
        StoreXmMatrix(xm, o);

        #else

        o.m[0][0] = m.m[0][0]; o.m[0][1] = m.m[1][0]; o.m[0][2] = m.m[2][0]; o.m[0][3] = m.m[3][0];
        o.m[1][0] = m.m[0][1]; o.m[1][1] = m.m[1][1]; o.m[1][2] = m.m[2][1]; o.m[1][3] = m.m[3][1];
        o.m[2][0] = m.m[0][2]; o.m[2][1] = m.m[1][2]; o.m[2][2] = m.m[2][2]; o.m[2][3] = m.m[3][2];
        o.m[3][0] = m.m[0][3]; o.m[3][1] = m.m[1][3]; o.m[3][2] = m.m[2][3]; o.m[3][3] = m.m[3][3];

        #endif
    }

    inline bool Matrix::Invert(const Matrix& m, Matrix& o)
    {
        #ifdef COMPILE_MATH_MATRIX_DIRECTX_MATH

        DirectX::XMMATRIX xm = LoadXmMatrix(m);
        xm = DirectX::XMMatrixInverse(nullptr, xm);
        if (DirectX::XMMatrixIsNaN(xm))
        {
            return false;
        }
        StoreXmMatrix(xm, o);
        return true;

        #else

        // Precomputationed values.
        const float t00 = m.m[2][2] * m.m[3][3] - m.m[2][3] * m.m[3][2];
        const float t01 = m.m[2][3] * m.m[3][0] - m.m[2][0] * m.m[3][3];
        const float t02 = m.m[2][1] * m.m[3][2] - m.m[2][2] * m.m[3][1];
        const float t03 = m.m[2][3] * m.m[3][1] - m.m[2][1] * m.m[3][3];
        const float t04 = m.m[2][0] * m.m[3][1] - m.m[2][1] * m.m[3][0];
        const float t05 = m.m[2][0] * m.m[3][2] - m.m[2][2] * m.m[3][0];
        const float t06 = m.m[2][1] * m.m[3][3] - m.m[2][3] * m.m[3][1];
        const float t07 = m.m[2][2] * m.m[3][0] - m.m[2][0] * m.m[3][2];
        const float t08 = m.m[1][2] * m.m[3][3] - m.m[1][3] * m.m[3][2];
        const float t09 = m.m[1][1] * m.m[3][2] - m.m[1][2] * m.m[3][1];
        const float t10 = m.m[1][0] * m.m[3][1] - m.m[1][1] * m.m[3][0];
        const float t11 = m.m[1][2] * m.m[2][3] - m.m[1][3] * m.m[2][2];
        const float t12 = m.m[1][1] * m.m[2][2] - m.m[1][2] * m.m[2][1];
        const float t13 = m.m[1][3] * m.m[2][0] - m.m[1][0] * m.m[2][3];
        const float t14 = m.m[1][0] * m.m[2][1] - m.m[1][1] * m.m[2][0];
        const float t15 = m.m[1][3] * m.m[3][0] - m.m[1][0] * m.m[3][3];

        // Minors.
        float det00 = m.m[1][1] * t00 + m.m[1][2] * t03 + m.m[1][3] * t02;
        float det01 = m.m[1][0] * t00 + m.m[1][2] * t01 + m.m[1][3] * t05;
        float det02 = m.m[1][0] * t06 + m.m[1][1] * t01 + m.m[1][3] * t04;
        float det03 = m.m[1][0] * t02 + m.m[1][1] * t07 + m.m[1][2] * t04;
        float det10 = m.m[0][1] * t00 + m.m[0][2] * t03 + m.m[0][3] * t02;
        float det11 = m.m[0][0] * t00 + m.m[0][2] * t01 + m.m[0][3] * t05;
        float det12 = m.m[0][0] * t06 + m.m[0][1] * t01 + m.m[0][3] * t04;
        float det13 = m.m[0][0] * t02 + m.m[0][1] * t07 + m.m[0][2] * t04;
        float det20 = m.m[0][1] * t08 + m.m[0][3] * t09 + m.m[0][2] * (m.m[1][3] * m.m[3][1] - m.m[1][1] * m.m[3][3]);
        float det21 = m.m[0][0] * t08 + m.m[0][2] * t15 + m.m[0][3] * (m.m[1][0] * m.m[3][2] - m.m[1][2] * m.m[3][0]);
        float det22 = m.m[0][1] * t15 + m.m[0][3] * t10 + m.m[0][0] * (m.m[1][1] * m.m[3][3] - m.m[1][3] * m.m[3][1]);
        float det23 = m.m[0][0] * t09 + m.m[0][2] * t10 + m.m[0][1] * (m.m[1][2] * m.m[3][0] - m.m[1][0] * m.m[3][2]);
        float det30 = m.m[0][1] * t11 + m.m[0][3] * t12 + m.m[0][2] * (m.m[1][3] * m.m[2][1] - m.m[1][1] * m.m[2][3]);
        float det31 = m.m[0][0] * t11 + m.m[0][2] * t13 + m.m[0][3] * (m.m[1][0] * m.m[2][2] - m.m[1][2] * m.m[2][0]);
        float det32 = m.m[0][1] * t13 + m.m[0][3] * t14 + m.m[0][0] * (m.m[1][1] * m.m[2][3] - m.m[1][3] * m.m[2][1]);
        float det33 = m.m[0][0] * t12 + m.m[0][2] * t14 + m.m[0][1] * (m.m[1][2] * m.m[2][0] - m.m[1][0] * m.m[2][2]);

        // The matrix determinant.
        float det = m.m[0][0] * det00 - m.m[1][0] * det10 + m.m[2][0] * det20 - m.m[3][0] * det30;

        // The matrix is singular (not invertible).
        if (det == 0.0f)
        {
            return false;
        }

        // Set inverse entries.
        o.m[0][0] = det00 / det;
        o.m[0][1] = -det10 / det;
        o.m[0][2] = det20 / det;
        o.m[0][3] = -det30 / det;
        o.m[1][0] = -det01 / det;
        o.m[1][1] = det11 / det;
        o.m[1][2] = -det21 / det;
        o.m[1][3] = det31 / det;
        o.m[2][0] = det02 / det;
        o.m[2][1] = -det12 / det;
        o.m[2][2] = det22 / det;
        o.m[2][3] = -det32 / det;
        o.m[3][0] = -det03 / det;
        o.m[3][1] = det13 / det;
        o.m[3][2] = -det23 / det;
        o.m[3][3] = det33 / det;

        return true;

        #endif
    }

    inline Matrix Matrix::Multiply(const Matrix& l, const Matrix& r)
    {
        #ifdef COMPILE_MATH_MATRIX_DIRECTX_MATH

        const DirectX::XMMATRIX xm = DirectX::XMMatrixMultiply(LoadXmMatrix(l), LoadXmMatrix(r));
        Matrix m;
        StoreXmMatrix(xm, m);
        return m;

        #else

        Matrix result;
        result.Zero();

        for (int row = 0; row < 4; ++row)
        {
            for (int col = 0; col < 4; ++col)
            {
                for (int i = 0; i < 4; ++i)
                {
                    result.m[row][col] += l.m[row][i] * r.m[i][col];
                }
            }
        }

        return result;

        #endif
    }

    inline Vector Matrix::Multiply(const Matrix& m, const Vector& v)
    {
        #ifdef COMPILE_MATH_MATRIX_DIRECTX_MATH

        DirectX::XMMATRIX xm = LoadXmMatrix(m);
        xm = DirectX::XMMatrixTranspose(xm);
        DirectX::XMVECTOR xv1 = LoadXmVector(v);
        DirectX::XMVECTOR xv2 = DirectX::XMVector4Transform(xv1, xm);
        return CreateVector(xv2);

        #else

        return Vector(
            m.m[0][0] * v.x + m.m[0][1] * v.y + m.m[0][2] * v.z + m.m[0][3] * v.w,
            m.m[1][0] * v.x + m.m[1][1] * v.y + m.m[1][2] * v.z + m.m[1][3] * v.w,
            m.m[2][0] * v.x + m.m[2][1] * v.y + m.m[2][2] * v.z + m.m[2][3] * v.w,
            m.m[3][0] * v.x + m.m[3][1] * v.y + m.m[3][2] * v.z + m.m[3][3] * v.w
        );

        #endif
    }

    inline Vector Matrix::Multiply(const Vector& v, const Matrix& m)
    {
        #ifdef COMPILE_MATH_MATRIX_DIRECTX_MATH

        DirectX::XMMATRIX xm = LoadXmMatrix(m);
        DirectX::XMVECTOR xv1 = LoadXmVector(v);
        DirectX::XMVECTOR xv2 = DirectX::XMVector4Transform(xv1, xm);
        return CreateVector(xv2);

        #else

        return Vector(
            v.x * m.m[0][0] + v.y * m.m[1][0] + v.z * m.m[2][0] + v.w * m.m[3][0],
            v.x * m.m[0][1] + v.y * m.m[1][1] + v.z * m.m[2][1] + v.w * m.m[3][1],
            v.x * m.m[0][2] + v.y * m.m[1][2] + v.z * m.m[2][2] + v.w * m.m[3][2],
            v.x * m.m[0][3] + v.y * m.m[1][3] + v.z * m.m[2][3] + v.w * m.m[3][3]
        );

        #endif
    }

    inline void Matrix::Add(const Matrix& l, const Matrix& r, Matrix& result)
    {
        for (int row = 0; row < 4; ++row)
        {
            for (int col = 0; col < 4; ++col)
            {
                result.m[row][col] = l.m[row][col] + r.m[row][col];
            }
        }
    }

    // Matrix * column-major vector (vector transformation).
    inline Vector operator*(const Matrix& m, const Vector& v)
    {
        return Matrix::Multiply(m, v);
    }

    // Row-major vector * matrix (vector transformation).
    inline Vector operator*(const Vector& v, const Matrix& m)
    {
        return Matrix::Multiply(v, m);
    }

    inline Matrix Matrix::Identity()
    {
        Matrix m;
        m.m[0][0] = 1.0f; m.m[0][1] = 0.0f; m.m[0][2] = 0.0f; m.m[0][3] = 0.0f;
        m.m[1][0] = 0.0f; m.m[1][1] = 1.0f; m.m[1][2] = 0.0f; m.m[1][3] = 0.0f;
        m.m[2][0] = 0.0f; m.m[2][1] = 0.0f; m.m[2][2] = 1.0f; m.m[2][3] = 0.0f;
        m.m[3][0] = 0.0f; m.m[3][1] = 0.0f; m.m[3][2] = 0.0f; m.m[3][3] = 1.0f;
        return m;
    }

    inline Matrix Matrix::Transpose(const Matrix& m)
    {
        Matrix n;
        Transpose(m, n);
        return n;
    }

    inline Matrix Matrix::Inverse(const Matrix& m)
    {
        Matrix n;
        Invert(m, n);
        return n;
    }

    inline Matrix& Matrix::operator=(const Matrix& matrix)
    {
        memcpy(m, matrix.m, sizeof(m));
        return *this;
    }

    inline Matrix Matrix::Zero()
    {
        Matrix m;
        memset(&m, 0, sizeof(m));
        return m;
    }

    inline void Matrix::Transpose()
    {
        MatrixTranspose(*this);
    }

    inline bool Matrix::operator==(const Matrix& r) const
    {
        return std::memcmp(m, r.m, sizeof(m)) == 0;
    }

    inline bool Matrix::operator!=(const Matrix& r) const
    {
        return std::memcmp(m, r.m, sizeof(m)) != 0;
    }

    inline Matrix Matrix::operator*(const Matrix& r) const
    {
        return Multiply(*this, r);
    }

    inline Matrix Matrix::operator+(const Matrix& r) const
    {
        Matrix result;
        Add(*this, r, result);
        return result;
    }

    inline Matrix& Matrix::operator+=(const Matrix& r)
    {
        Add(*this, r, *this);
        return *this;
    }

    inline void Matrix::SwapRows(const int r1, const int r2)
    {
        std::swap(m[r1][0], m[r2][0]);
        std::swap(m[r1][1], m[r2][1]);
        std::swap(m[r1][2], m[r2][2]);
        std::swap(m[r1][3], m[r2][3]);
    }

    inline void Matrix::AddRow(const int srcRow, const int destRow, const float multiplier)
    {
        m[destRow][0] += (m[srcRow][0] * multiplier);
        m[destRow][1] += (m[srcRow][1] * multiplier);
        m[destRow][2] += (m[srcRow][2] * multiplier);
        m[destRow][3] += (m[srcRow][3] * multiplier);
    }

    inline void Matrix::MulRow(const int row, const float value)
    {
        m[row][0] *= value;
        m[row][1] *= value;
        m[row][2] *= value;
        m[row][3] *= value;
    }

    inline bool Matrix::Invert()
    {
        return Invert(*this, *this);
    }

    inline float Matrix::Determinant() const
    {
        return Determinant(*this);
    }

    inline Vector Matrix::Transform(const Vector& v) const
    {
        return *this * v;
    }

    inline void Matrix::AppendTransformations(const Matrix& t)
    {
        *this = Multiply(t, *this);
    }

    inline void Matrix::PrependTransformations(const Matrix& t)
    {
        *this = Multiply(*this, t);
    }

    inline void Matrix::AppendTranslation(const float x, const float y, const float z)
    {
        m[0][0] = x * m[3][0] + m[0][0];
        m[0][1] = x * m[3][1] + m[0][1];
        m[0][2] = x * m[3][2] + m[0][2];
        m[0][3] = x * m[3][3] + m[0][3];
        m[1][0] = y * m[3][0] + m[1][0];
        m[1][1] = y * m[3][1] + m[1][1];
        m[1][2] = y * m[3][2] + m[1][2];
        m[1][3] = y * m[3][3] + m[1][3];
        m[2][0] = z * m[3][0] + m[2][0];
        m[2][1] = z * m[3][1] + m[2][1];
        m[2][2] = z * m[3][2] + m[2][2];
        m[2][3] = z * m[3][3] + m[2][3];
    }

    inline void Matrix::PrependTranslation(const float x, const float y, const float z)
    {
        m[0][3] = x * m[0][0] + y * m[0][1] + z * m[0][2] + m[0][3];
        m[1][3] = x * m[1][0] + y * m[1][1] + z * m[1][2] + m[1][3];
        m[2][3] = x * m[2][0] + y * m[2][1] + z * m[2][2] + m[2][3];
        m[3][3] = x * m[3][0] + y * m[3][1] + z * m[3][2] + m[3][3];
    }

    inline void Matrix::AppendScaling(const float x, const float y, const float z)
    {
        m[0][0] *= x;
        m[0][1] *= x;
        m[0][2] *= x;
        m[0][3] *= x;
        m[1][0] *= y;
        m[1][1] *= y;
        m[1][2] *= y;
        m[1][3] *= y;
        m[2][0] *= z;
        m[2][1] *= z;
        m[2][2] *= z;
        m[2][3] *= z;
    }

    inline void Matrix::PrependScaling(const float x, const float y, const float z)
    {
        m[0][0] *= x;
        m[0][1] *= y;
        m[0][2] *= z;
        m[1][0] *= x;
        m[1][1] *= y;
        m[1][2] *= z;
        m[2][0] *= x;
        m[2][1] *= y;
        m[2][2] *= z;
        m[3][0] *= x;
        m[3][1] *= y;
        m[3][2] *= z;
    }

    inline Matrix Matrix::RotationX(const float rad)
    {
        const float sin = std::sinf(rad);
        const float cos = std::cosf(rad);
        Matrix m;
        m.m[0][0] = 1.0f; m.m[0][1] = 0.0f; m.m[0][2] = 0.0f; m.m[0][3] = 0.0f;
        m.m[1][0] = 0.0f; m.m[1][1] = cos;  m.m[1][2] = -sin; m.m[1][3] = 0.0f;
        m.m[2][0] = 0.0f; m.m[2][1] = sin;  m.m[2][2] = cos;  m.m[2][3] = 0.0f;
        m.m[3][0] = 0.0f; m.m[3][1] = 0.0f; m.m[3][2] = 0.0f; m.m[3][3] = 1.0f;
        return m;
    }

    inline Matrix Matrix::RotationY(const float rad)
    {
        const float sin = std::sinf(rad);
        const float cos = std::cosf(rad);
        Matrix m;
        m.m[0][0] = cos;  m.m[0][1] = 0.0f; m.m[0][2] = sin;  m.m[0][3] = 0.0f;
        m.m[1][0] = 0.0f; m.m[1][1] = 1.0f; m.m[1][2] = 0.0f; m.m[1][3] = 0.0f;
        m.m[2][0] = -sin; m.m[2][1] = 0.0f; m.m[2][2] = cos;  m.m[2][3] = 0.0f;
        m.m[3][0] = 0.0f; m.m[3][1] = 0.0f; m.m[3][2] = 0.0f; m.m[3][3] = 1.0f;
        return m;
    }

    inline Matrix Matrix::RotationZ(const float rad)
    {
        const float sin = std::sinf(rad);
        const float cos = std::cosf(rad);
        Matrix m;
        m.m[0][0] = cos;  m.m[0][1] = -sin; m.m[0][2] = 0.0f; m.m[0][3] = 0.0f;
        m.m[1][0] = sin;  m.m[1][1] = cos;  m.m[1][2] = 0.0f; m.m[1][3] = 0.0f;
        m.m[2][0] = 0.0f; m.m[2][1] = 0.0f; m.m[2][2] = 1.0f; m.m[2][3] = 0.0f;
        m.m[3][0] = 0.0f; m.m[3][1] = 0.0f; m.m[3][2] = 0.0f; m.m[3][3] = 1.0f;
        return m;
    }

    inline Matrix Matrix::Rotation(const float x, const float y, const float z)
    {
        const float sx = std::sinf(x);
        const float cx = std::cosf(x);
        const float sy = std::sinf(y);
        const float cy = std::cosf(y);
        const float sz = std::sinf(z);
        const float cz = std::cosf(z);

        Matrix m;
        m.m[0][0] = sx * sy * sz + cy * cz;
        m.m[0][1] = sx * sy * cz - sz * cy;
        m.m[0][2] = sy * cx;
        m.m[0][3] = 0.0f;
        m.m[1][0] = sz * cx;
        m.m[1][1] = cx * cz;
        m.m[1][2] = -sx;
        m.m[1][3] = 0.0f;
        m.m[2][0] = sx * sz * cy - sy * cz;
        m.m[2][1] = sx * cy * cz + sy * sz;
        m.m[2][2] = cx * cy;
        m.m[2][3] = 0.0f;
        m.m[3][0] = 0.0f;
        m.m[3][1] = 0.0f;
        m.m[3][2] = 0.0f;
        m.m[3][3] = 1.0f;
        return m;
    }

    // Rotation decomposition
    // float rz = std::atan2f(rotations.m[1][0], rotations.m[1][1]);
    // float ry = std::atan2f(rotations.m[0][2], rotations.m[2][2]);

    inline Matrix Matrix::Scale(const float x, const float y, const float z)
    {
        Matrix m;
        m.m[0][0] = x;    m.m[0][1] = 0.0f; m.m[0][2] = 0.0f; m.m[0][3] = 0.0f;
        m.m[1][0] = 0.0f; m.m[1][1] = y;    m.m[1][2] = 0.0f; m.m[1][3] = 0.0f;
        m.m[2][0] = 0.0f; m.m[2][1] = 0.0f; m.m[2][2] = z;    m.m[2][3] = 0.0f;
        m.m[3][0] = 0.0f; m.m[3][1] = 0.0f; m.m[3][2] = 0.0f; m.m[3][3] = 1.0f;
        return m;
    }

    inline Matrix Matrix::Translation(const float x, const float y, const float z)
    {
        Matrix m;
        m.m[0][0] = 1.0f; m.m[0][1] = 0.0f; m.m[0][2] = 0.0f; m.m[0][3] = x;
        m.m[1][0] = 0.0f; m.m[1][1] = 1.0f; m.m[1][2] = 0.0f; m.m[1][3] = y;
        m.m[2][0] = 0.0f; m.m[2][1] = 0.0f; m.m[2][2] = 1.0f; m.m[2][3] = z;
        m.m[3][0] = 0.0f; m.m[3][1] = 0.0f; m.m[3][2] = 0.0f; m.m[3][3] = 1.0f;
        return m;
    }

} // namespace Math
#pragma once

#include <algorithm>
#include "Vector.h"

namespace Math
{
    // 4x4 matrix.
    class alignas(16) Matrix
    {
    public:
        // Column major entries.
        float m[4][4];

        Matrix() = default;
        Matrix(const Matrix&) = default;
        Matrix& operator=(const Matrix&);

        static Matrix MakeIdentity();
        static Matrix MakeTransposed(const Matrix& m);
        static Matrix MakeInverted(const Matrix& m);

        //static Matrix MakeRotationX(const float rad);
        //static Matrix MakeRotationY(const float rad);
        //static Matrix MakeRotationZ(const float rad);
        //static Matrix MakeScale(const float scaleX, const float scaleY, const float scaleZ);
        //static Matrix MakeTranslation(const float x, const float y, const float z);

        void Identity();
        void Zero();
        void Transpose();
        bool Invert();

        float Determinant() const;

        // Equality.
        bool operator==(const Matrix&) const;
        bool operator!=(const Matrix&) const;

        // Add.
        Matrix operator+(const Matrix&) const;
        Matrix& operator+=(const Matrix&);

        // Multiply.
        Matrix operator*(const Matrix&) const;

        // Swap rows.
        void SwapRows(const int r1, const int r2);

        // Add a multiplied src row to a dest row.
        void AddRow(const int srcRow, const int destRow, const float multiplier);

        // Multiply all entries of a row.
        void MulRow(const int row, const float value);

        // Apply transformations to column-major vector.
        Vector Transform(const Vector& v) const;

        void AppendTransformations(const Matrix&);
        void PrependTransformations(const Matrix&);

        void AppendTranslation(const float x, const float y, const float z);
        void PrependTranslation(const float x, const float y, const float z);

        void AppendScaling(const float x, const float y, const float z);
        void PrependScaling(const float x, const float y, const float z);
    };
}

// Include implementation file.
#include "Matrix.inl"
//#include "MatrixDirectXMath.inl"

namespace Math
{
    // Matrix * column-major vector (vector transformation).
    inline Vector operator*(const Matrix& m, const Vector& v)
    {
        return MatrixMultiply(m, v);
    }

    // Row-major vector * matrix (vector transformation).
    inline Vector operator*(const Vector& v, const Matrix& m)
    {
        return MatrixMultiply(v, m);
    }

    inline Matrix Matrix::MakeIdentity()
    {
        Matrix m;
        m.m[0][0] = 1.0f; m.m[0][1] = 0.0f; m.m[0][2] = 0.0f; m.m[0][3] = 0.0f;
        m.m[1][0] = 0.0f; m.m[1][1] = 1.0f; m.m[1][2] = 0.0f; m.m[1][3] = 0.0f;
        m.m[2][0] = 0.0f; m.m[2][1] = 0.0f; m.m[2][2] = 1.0f; m.m[2][3] = 0.0f;
        m.m[3][0] = 0.0f; m.m[3][1] = 0.0f; m.m[3][2] = 0.0f; m.m[3][3] = 1.0f;
        return m;
    }

    inline Matrix Matrix::MakeTransposed(const Matrix& m)
    {
        Matrix n;
        MatrixTranspose(m, n);
        return n;
    }
    
    inline Matrix Matrix::MakeInverted(const Matrix& m)
    {
        Matrix n;
        MatrixInvert(m, n);
        return n;
    }

    inline Matrix& Matrix::operator=(const Matrix& matrix)
    {
        memcpy(m, matrix.m, sizeof(m));
        return *this;
    }

    inline void Matrix::Zero()
    {
        memset(m, 0, sizeof(m));
    }

    inline void Matrix::Identity()
    {
        m[0][0] = 1; m[0][1] = 0; m[0][2] = 0; m[0][3] = 0;
        m[1][0] = 0; m[1][1] = 1; m[1][2] = 0; m[1][3] = 0;
        m[2][0] = 0; m[2][1] = 0; m[2][2] = 1; m[2][3] = 0;
        m[3][0] = 0; m[3][1] = 0; m[3][2] = 0; m[3][3] = 1;
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
        return MatrixMultiply(*this, r);
    }

    inline Matrix Matrix::operator+(const Matrix& r) const
    {
        Matrix result;
        MatrixAdd(*this, r, result);
        return result;
    }

    inline Matrix& Matrix::operator+=(const Matrix& r)
    {
        MatrixAdd(*this, r, *this);
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
        MatrixMulRow(*this, row, value);
    }

    inline bool Matrix::Invert()
    {
        return MatrixInvert(*this, *this);
    }

    inline float Matrix::Determinant() const
    {
        return MatrixDeterminant(*this);
    }

    inline Vector Matrix::Transform(const Vector& v) const
    {
        return v * (*this);
    }

    inline void Matrix::AppendTransformations(const Matrix& t)
    {
        *this = MatrixMultiply(t, *this);
    }
    
    inline void Matrix::PrependTransformations(const Matrix& t)
    {
        *this = MatrixMultiply(*this, t);
    }

    inline void Matrix::AppendTranslation(const float x, const float y, const float z)
    {
        MatrixAppendTranslation(*this, x, y, z);
    }

    inline void Matrix::PrependTranslation(const float x, const float y, const float z)
    {
        MatrixPrependTranslation(*this, x, y, z);
    }

    inline void Matrix::AppendScaling(const float x, const float y, const float z)
    {
        MatrixAppendScaling(*this, x, y, z);
    }

    inline void Matrix::PrependScaling(const float x, const float y, const float z)
    {
        MatrixPrependScaling(*this, x, y, z);
    }
}
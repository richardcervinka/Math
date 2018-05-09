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

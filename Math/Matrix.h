#pragma once

#include <algorithm>
#include "Vector.h"

namespace Math
{
    // 4x4 matrix.
    class alignas(16) Matrix
    {
    public:
        // Row major entries: m[row][column].
        float m[4][4];

        Matrix() = default;
        Matrix(const Matrix&) = default;
        Matrix& operator=(const Matrix&);

        static Matrix Zero();
        static Matrix Identity();
        static Matrix Transpose(const Matrix& m);
        static Matrix Inverse(const Matrix& m);

        //static Matrix RotationX(const float rad);
        //static Matrix RotationY(const float rad);
        //static Matrix RotationZ(const float rad);
        //static Matrix Scale(const float scaleX, const float scaleY, const float scaleZ);
        //static Matrix Translation(const float x, const float y, const float z);

        void Transpose();
        bool Invert();

        // Determinant.
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

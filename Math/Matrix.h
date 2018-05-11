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

        Matrix(const Matrix&) = default;
        Matrix& operator=(const Matrix&);

        // Create a zero matrix.
        static Matrix Zero();

        // Create an identity matrix.
        static Matrix Identity();

        // Create transpose matrix.
        static Matrix Transpose(const Matrix& m);

        // Create inverse matrix.
        static Matrix Inverse(const Matrix& m);

        // Create rotation matrices around the world axes.
        static Matrix RotationX(const float rad);
        static Matrix RotationY(const float rad);
        static Matrix RotationZ(const float rad);

        // Create a rotation matrix around all worl axes in order Z, X, Y.
        static Matrix Rotation(const float x, const float y, const float z);

        // Create a rotation matrix around the axis defined as unit vector.
        //static Matrix RotationAxes(const Vector& v);

        // Create a scale matrix.
        static Matrix Scale(const float x, const float y, const float z);

        // Create a translation matrix.
        static Matrix Translation(const float x, const float y, const float z);

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

        // Apply transformations (this * v).
        Vector Transform(const Vector& v) const;

        void AppendTransformations(const Matrix&);
        void PrependTransformations(const Matrix&);

        void AppendTranslation(const float x, const float y, const float z);
        void PrependTranslation(const float x, const float y, const float z);

        void AppendScaling(const float x, const float y, const float z);
        void PrependScaling(const float x, const float y, const float z);

    private:
        // Create a uninitialized matrix. Only internal use.
        Matrix() = default;

        // Implementation.
        static float Determinant(const Matrix& m);
        static void MatrixTranspose(Matrix& m);
        static void Transpose(const Matrix& m, Matrix& o);
        static bool Invert(const Matrix& m, Matrix& o);
        static Matrix Multiply(const Matrix& l, const Matrix& r);
        static Vector Multiply(const Matrix& m, const Vector& v);
        static Vector Multiply(const Vector& v, const Matrix& m);
        static void Add(const Matrix& l, const Matrix& r, Matrix& result);

        // Friends that needs access to the Multiply function.
        friend Vector operator*(const Matrix& m, const Vector& v);
        friend Vector operator*(const Vector& v, const Matrix& m);
    };
}

// Implementation file.
#include "Matrix.inl"

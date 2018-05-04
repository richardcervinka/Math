#pragma once

inline float MatrixDeterminant(const Math::Matrix& m)
{
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
}

inline void MatrixTranspose(Math::Matrix& m)
{
    std::swap(m.m[1][0], m.m[0][1]);
    std::swap(m.m[2][0], m.m[0][2]);
    std::swap(m.m[3][0], m.m[0][3]);
    std::swap(m.m[2][1], m.m[1][2]);
    std::swap(m.m[3][1], m.m[1][3]);
    std::swap(m.m[3][2], m.m[2][3]);
}

inline void MatrixTranspose(const Math::Matrix& m, Math::Matrix& o)
{
    o.m[0][0] = m.m[0][0]; o.m[0][1] = m.m[1][0]; o.m[0][2] = m.m[2][0]; o.m[0][3] = m.m[3][0];
    o.m[1][0] = m.m[0][1]; o.m[1][1] = m.m[1][1]; o.m[1][2] = m.m[2][1]; o.m[1][3] = m.m[3][1];
    o.m[2][0] = m.m[0][2]; o.m[2][1] = m.m[1][2]; o.m[2][2] = m.m[2][2]; o.m[2][3] = m.m[3][2];
    o.m[3][0] = m.m[0][3]; o.m[3][1] = m.m[1][3]; o.m[3][2] = m.m[2][3]; o.m[3][3] = m.m[3][3];
}

inline bool MatrixInvert(const Math::Matrix& m, Math::Matrix& o)
{
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
}

inline Math::Matrix MatrixMultiply(const Math::Matrix& l, const Math::Matrix& r)
{
    Math::Matrix result;
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
}

// Matrix * column-major vector (vector transformation).
inline Math::Vector MatrixMultiply(const Math::Matrix& m, const Math::Vector& v)
{
    return Math::Vector(
        m.m[0][0] * v.x + m.m[0][1] * v.y + m.m[0][2] * v.z + m.m[0][3] * v.w,
        m.m[1][0] * v.x + m.m[1][1] * v.y + m.m[1][2] * v.z + m.m[1][3] * v.w,
        m.m[2][0] * v.x + m.m[2][1] * v.y + m.m[2][2] * v.z + m.m[2][3] * v.w,
        m.m[3][0] * v.x + m.m[3][1] * v.y + m.m[3][2] * v.z + m.m[3][3] * v.w
    );
}

// Row-major vector * matrix (vector transformation).
inline Math::Vector MatrixMultiply(const Math::Vector& v, const Math::Matrix& m)
{
    return Math::Vector(
        v.x * m.m[0][0] + v.y * m.m[1][0] + v.z * m.m[2][0] + v.w * m.m[3][0],
        v.x * m.m[0][1] + v.y * m.m[1][1] + v.z * m.m[2][1] + v.w * m.m[3][1],
        v.x * m.m[0][2] + v.y * m.m[1][2] + v.z * m.m[2][2] + v.w * m.m[3][2],
        v.x * m.m[0][3] + v.y * m.m[1][3] + v.z * m.m[2][3] + v.w * m.m[3][3]
    );
}

inline void MatrixAdd(const Math::Matrix& l, const Math::Matrix& r, Math::Matrix& result)
{
    for (int row = 0; row < 4; ++row)
    {
        for (int col = 0; col < 4; ++col)
        {
            result.m[row][col] = l.m[row][col] + r.m[row][col];
        }
    }
}

inline void MatrixMulRow(Math::Matrix& m, const int row, const float value)
{
    m.m[row][0] *= value;
    m.m[row][1] *= value;
    m.m[row][2] *= value;
    m.m[row][3] *= value;
}

inline void MatrixAppendTranslation(Math::Matrix& m, const float x, const float y, const float z)
{
    m.m[0][0] = x * m.m[3][0] + m.m[0][0];
    m.m[0][1] = x * m.m[3][1] + m.m[0][1];
    m.m[0][2] = x * m.m[3][2] + m.m[0][2];
    m.m[0][3] = x * m.m[3][3] + m.m[0][3];
    m.m[1][0] = y * m.m[3][0] + m.m[1][0];
    m.m[1][1] = y * m.m[3][1] + m.m[1][1];
    m.m[1][2] = y * m.m[3][2] + m.m[1][2];
    m.m[1][3] = y * m.m[3][3] + m.m[1][3];
    m.m[2][0] = z * m.m[3][0] + m.m[2][0];
    m.m[2][1] = z * m.m[3][1] + m.m[2][1];
    m.m[2][2] = z * m.m[3][2] + m.m[2][2];
    m.m[2][3] = z * m.m[3][3] + m.m[2][3];
}

inline void MatrixPrependTranslation(Math::Matrix& m, const float x, const float y, const float z)
{
    m.m[0][3] = x * m.m[0][0] + y * m.m[0][1] + z * m.m[0][2] + m.m[0][3];
    m.m[1][3] = x * m.m[1][0] + y * m.m[1][1] + z * m.m[1][2] + m.m[1][3];
    m.m[2][3] = x * m.m[2][0] + y * m.m[2][1] + z * m.m[2][2] + m.m[2][3];
    m.m[3][3] = x * m.m[3][0] + y * m.m[3][1] + z * m.m[3][2] + m.m[3][3];
}

inline void MatrixAppendScaling(Math::Matrix& m, const float x, const float y, const float z)
{
    m.m[0][0] *= x;
    m.m[0][1] *= x;
    m.m[0][2] *= x;
    m.m[0][3] *= x;
    m.m[1][0] *= y;
    m.m[1][1] *= y;
    m.m[1][2] *= y;
    m.m[1][3] *= y;
    m.m[2][0] *= z;
    m.m[2][1] *= z;
    m.m[2][2] *= z;
    m.m[2][3] *= z;
}

inline void MatrixPrependScaling(Math::Matrix& m, const float x, const float y, const float z)
{
    m.m[0][0] *= x;
    m.m[0][1] *= y;
    m.m[0][2] *= z;
    m.m[1][0] *= x;
    m.m[1][1] *= y;
    m.m[1][2] *= z;
    m.m[2][0] *= x;
    m.m[2][1] *= y;
    m.m[2][2] *= z;
    m.m[3][0] *= x;
    m.m[3][1] *= y;
    m.m[3][2] *= z;
}

inline void MatrixAppendRotationX(Math::Matrix& m, const float rad)
{
    float s = std::sinf(rad);
    float c = std::cosf(rad);

    // To do ...
}
#pragma once

#include <cmath>

// Column vectors.
// Row major matrix.
// Right-hand coordinate system.

namespace Math
{
    // Math constants.
    const float Pi = 3.1415926535898f;
    const float Pi2 = Pi * 2.0f;
    const float PiHalf = Pi / 2.0f;
}

#include "Vector.h"
#include "Matrix.h"

namespace Math
{
    // Linear interpolation.
    inline float Interpolate(const float from, const float to, const float ratio)
    {
        return from + (to - from) * ratio;
    }

    // Linear interpolation.
    inline Vector Interpolate(const Vector& from, const Vector& to, const float ratio)
    {
        return Vector(
            from.x + (to.x - from.x) * ratio,
            from.y + (to.y - from.y) * ratio,
            from.z + (to.z - from.z) * ratio,
            from.w + (to.w - from.w) * ratio
        );
    }
}
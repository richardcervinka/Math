#pragma once

#include <cmath>

namespace Math
{
    class alignas(16) Vector
    {
    public:
        float x = 0.0f;
        float y = 0.0f;
        float z = 0.0f;
        float w = 0.0f;

        // Construction.
        Vector() = default;
        Vector(const float x_, const float y_, const float z_, const float w_): x(x_), y(y_), z(z_), w(w_) { }
        Vector(const Vector&) = default;

        // Default destructor.
        ~Vector() = default;

        // Assignment.
        Vector& operator=(const Vector&) = default;
        void Set(const float x_, const float y_, const float z_, const float w_) noexcept;

        // Returns normalized vector.
        static Vector Normalized(const Vector&) noexcept;
        static Vector Normalized(Vector&&) noexcept;

        // Returns normal vector.
        static Vector Normal(const Vector& v1, const Vector& v2) noexcept;

        // Returns orthogonal projection v1 to the v2.
        static Vector Projection(const Vector& v1, const Vector& v2) noexcept;

        // Returns length of the result of the orthogonal projection v1 to the v2.
        static float ProjectionLength(const Vector& v1, const Vector& v2) noexcept;

        // Angle between two vectors.
        static float Angle(const Vector&, const Vector&) noexcept;

        // Cosinus of the angle between two vectors.
        static float CosAngle(const Vector&, const Vector&) noexcept;

        // Dot product.
        static float Dot(const Vector&, const Vector&) noexcept;

        // Cross product.
        static Vector Cross(const Vector&, const Vector&) noexcept;

        // Length of the vector.
        static float Length(const Vector&) noexcept;

        // Squared length of the vector.
        static float LengthSq(const Vector&) noexcept;

        // Comparsion two vectors.
        bool operator==(const Vector&) const noexcept;
        bool operator!=(const Vector&) const noexcept;

        // Comparsion length of the vector and scalar.
        bool operator==(const float) const noexcept;
        bool operator!=(const float) const noexcept;
        bool operator>(const float) const noexcept;
        bool operator>=(const float) const noexcept;
        bool operator<(const float) const noexcept;
        bool operator<=(const float) const noexcept;

        // Comparsion length of the two vectors.
        bool operator>(const Vector&) const noexcept;
        bool operator>=(const Vector&) const noexcept;
        bool operator<(const Vector&) const noexcept;
        bool operator<=(const Vector&) const noexcept;

        // Add.
        Vector operator+(const Vector&) const noexcept;
        Vector& operator+=(const Vector&) noexcept;

        // Subtract.
        Vector operator-(const Vector&) const noexcept;
        Vector& operator-=(const Vector&) noexcept;

        // Multiply.
        Vector operator*(const float) const noexcept;
        Vector& operator*=(const float) noexcept;

        // Divide.
        Vector operator/(const float) const noexcept;
        Vector& operator/=(const float) noexcept;

        // Negation.
        Vector operator-() const noexcept;

        // Dot product.
        float operator*(const Vector&) const noexcept;

        // Get length.
        float Length() const noexcept;

        // Get length square.
        float LengthSq() const noexcept;

        // Set length. If initial length is 0, does nothing.
        void SetLength(const float) noexcept;

        // Set length to 1.
        void Normalize() noexcept;

        // Vector negation.
        void Invert() noexcept;
    };

    // Vector members implementation.

    inline void Vector::Set(const float x_, const float y_, const float z_, const float w_) noexcept
    {
        x = x_;
        y = y_;
        z = z_;
        w = w_;
    }

    inline float Vector::operator*(const Vector& r) const noexcept
    {
        return (x * r.x) + (y * r.y) + (z * r.z) + (w * r.w);
    }

    inline Vector& Vector::operator+=(const Vector& r) noexcept
    {
        Set(x + r.x, y + r.y, z + r.z, w + r.w);
        return *this;
    }

    inline Vector Vector::operator+(const Vector& r) const noexcept
    {
        return Vector(x + r.x, y + r.y, z + r.z, w + r.w);
    }

    inline Vector& Vector::operator-=(const Vector& r) noexcept
    {
        Set(x - r.x, y - r.y, z - r.z, w - r.w);
        return *this;
    }

    inline Vector Vector::operator-(const Vector& r) const noexcept
    {
        return Vector(x - r.x, y - r.y, z - r.z, w - r.w);
    }

    inline Vector& Vector::operator*=(const float value) noexcept
    {
        Set(x * value, y * value, z * value, w * value);
        return *this;
    }

    inline Vector Vector::operator*(const float value) const noexcept
    {
        return Vector(x * value, y * value, z * value, w * value);
    }

    inline Vector& Vector::operator/=(const float value) noexcept
    {
        Set(x / value, y / value, z / value, w / value);
        return *this;
    }

    inline Vector Vector::operator/(const float value) const noexcept
    {
        return Vector(x / value, y / value, z / value, w / value);
    }

    inline Vector Vector::operator-() const noexcept
    {
        return Vector(-x, -y, -z, -w);
    }

    inline bool Vector::operator==(const Vector& r) const noexcept
    {
        return (x == r.x) && (y == r.y) && (z == r.z) && (w == r.w);
    }

    inline bool Vector::operator!=(const Vector& r) const noexcept
    {
        return (x != r.x) || (y != r.y) || (z != r.z) || (w != r.w);
    }

    inline bool Vector::operator==(const float length) const noexcept
    {
        return LengthSq() == (length * length);
    }

    inline bool Vector::operator!=(const float length) const noexcept
    {
        return LengthSq() != (length * length);
    }

    inline bool Vector::operator>(const float length) const noexcept
    {
        return LengthSq() > (length * length);
    }

    inline bool Vector::operator>=(const float length) const noexcept
    {
        return LengthSq() >= (length * length);
    }

    inline bool Vector::operator<(const float length) const noexcept
    {
        return LengthSq() < (length * length);
    }

    inline bool Vector::operator<=(const float length) const noexcept
    {
        return LengthSq() <= (length * length);
    }

    inline bool Vector::operator>(const Vector& r) const noexcept
    {
        return LengthSq() > r.LengthSq();
    }

    inline bool Vector::operator>=(const Vector& r) const noexcept
    {
        return LengthSq() >= r.LengthSq();
    }

    inline bool Vector::operator<(const Vector& r) const noexcept
    {
        return LengthSq() < r.LengthSq();
    }

    inline bool Vector::operator<=(const Vector& r) const noexcept
    {
        return LengthSq() <= r.LengthSq();
    }

    inline float Vector::Length() const noexcept
    {
        return sqrtf(x * x + y * y + z * z + w * w);
    }

    inline float Vector::LengthSq() const noexcept
    {
        return x * x + y * y + z * z + w * w;
    }

    inline void Vector::SetLength(const float value) noexcept
    {
        float l = Length();
        if (l != 0)
        {
            float n = value / l;
            Set(x * n, y * n, z * n, w * n);
        }
    }

    inline void Vector::Normalize() noexcept
    {
        float l = Length();
        if (l != 0)
        {
            Set(x / l, y / l, z / l, w / l);
        }
    }

    inline void Vector::Invert() noexcept
    {
        Set(-x, -y, -z, -w);
    }

    inline Vector Vector::Normalized(Vector&& v) noexcept
    {
        v.Normalize();
        return v;
    }

    inline Vector Vector::Normalized(const Vector& v) noexcept
    {
        return Normalized(Vector(v));
    }

    inline Vector Vector::Normal(const Vector& v1, const Vector& v2) noexcept
    {
        return Normalized(Cross(v1, v2));
    }

    inline Vector Vector::Projection(const Vector& v1, const Vector& v2) noexcept
    {
        return (v2 * (v1 * v2)) / v2.LengthSq();
    }

    inline float Vector::ProjectionLength(const Vector& v1, const Vector& v2) noexcept
    {
        return (v1 * v2) / v2.Length();
    }

    inline float Vector::Angle(const Vector& v1, const Vector& v2) noexcept
    {
        return acosf(CosAngle(v1, v2));
    }

    inline float Vector::CosAngle(const Vector& v1, const Vector& v2) noexcept
    {
        return (v1 * v2) / (v1.Length() * v2.Length());
    }

    inline float Vector::Dot(const Vector& l, const Vector& r) noexcept
    {
        return l * r;
    }

    inline Vector Vector::Cross(const Vector& l, const Vector& r) noexcept
    {
        return Vector(
            l.y * r.z - l.z * r.y,
            l.z * r.x - l.x * r.z,
            l.x * r.y - l.y * r.x,
            0
        );
    }

    inline float Vector::Length(const Vector& v) noexcept
    {
        return v.Length();
    }

    inline float Vector::LengthSq(const Vector& v) noexcept
    {
        return v.LengthSq();
    }

    // Non-member operators.

    inline bool operator==(const float length, const Vector& v) noexcept
    {
        return v == length;
    }

    inline bool operator!=(const float length, const Vector& v) noexcept
    {
        return v != length;
    }

    inline bool operator>(const float length, const Vector& v) noexcept
    {
        return v < length;
    }

    inline bool operator>=(const float length, const Vector& v) noexcept
    {
        return v <= length;
    }

    inline bool operator<(const float length, const Vector& v) noexcept
    {
        return v > length;
    }

    inline bool operator<=(const float length, const Vector& v) noexcept
    {
        return v >= length;
    }
}
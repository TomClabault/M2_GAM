#include "point.h"
#include "vector.h"

#include <cmath>

Vector operator-(const Point& a, const Point& b)
{
    return Vector(a.x - b.x, a.y - b.y, a.z - b.z);
}

Vector operator/(const Vector& vec, double k)
{
    return Vector(vec.x / k, vec.y / k, vec.z / k);
}

double length(const Vector& vec)
{
    return std::sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
}

Vector normalize(const Vector& vec)
{
    return vec / length(vec);
}

double dot (const Vector& a, const Vector& b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vector cross(const Vector& a, const Vector& b)
{
    return Vector(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}


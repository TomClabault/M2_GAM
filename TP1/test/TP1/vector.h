#ifndef VECTOR_H
#define VECTOR_H

#include "point.h"

#include <cmath>

class Point;

class Vector
{
public:
    Vector() : x(0), y(0), z(0) {}
    Vector(double x, double y, double z) : x(x), y(y), z(z) {}

    Vector operator+=(const Vector& other);

    double x, y, z;
};

Vector operator-(const Point& a, const Point& b);
Vector operator/(const Vector& vec, double k);
Vector operator*(const Vector& vec, double k);
Vector operator*(double k, const Vector& vec);
Vector operator+(const Vector& a, double k);

double length(const Vector& vec);
Vector normalize(const Vector& vec);
double dot (const Vector& a, const Vector& b);
Vector cross(const Vector& a, const Vector& b);

#endif // VECTOR_H

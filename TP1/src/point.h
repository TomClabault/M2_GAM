#ifndef POINT_H
#define POINT_H

#include <ostream>

#include "TP1/vector.h"

class Vector;

class Point
{
public:
    double x;
    double y;
    double z;

    Point():x(),y(),z() {}
    Point(float x_, float y_, float z_):x(x_),y(y_),z(z_) {}

    /**
     * @return A positive number if the points are oriented CCW, negative if CW
     */
    static float orientation_test(const Point& a, const Point& b, const Point& c);
    static int is_point_in_triangle(const Point& point_to_test, const Point& a, const Point& b, const Point& c);

    double operator[] (int index);

    Point& operator/=(double k);
    friend bool operator==(const Point& a, const Point& b);
    friend std::ostream& operator<<(std::ostream& os, const Point& point);
};

inline Point min(const Point& a, const Point& b)
{
    return Point(std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z));
}

inline Point max(const Point& a, const Point& b)
{
    return Point(std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z));
}

inline Point operator*(const Point& p, double k)
{
    return Point(p.x * k, p.y * k, p.z * k);
}

inline Point operator*(double k, const Point& p)
{
    return p * k;
}

inline Point operator*(const Point& p, const Point& p2)
{
    return Point(p.x + p2.x, p.y + p2.y, p.z + p2.z);
}

inline Point operator+(const Point& p, const Vector& v)
{
    return Point(p.x + v.x, p.y + v.y, p.z + v.z);
}

#endif // POINT_H

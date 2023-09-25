#include "point.h"

bool operator==(const Point& a, const Point& b)
{
    return (std::abs(a.x - b.x) < 1.0e-6f && std::abs(a.y - b.y) < 1.0e-6f && std::abs(a.z - b.z) < 1.0e-6f);
}

std::ostream& operator<<(std::ostream& os, const Point& point)
{
    os << "Point[" << point.x << ", " << point.y << ", " << point.z <<  "]";

    return os;
}

float Point::orientation_test(const Point &a, const Point &b, const Point &c)
{
    Vector ab = b - a;
    Vector ac = c - a;

    Vector crossed = cross(ab, ac);

    return crossed.z;
}

int Point::is_point_in_triangle(const Point &point_to_test, const Point &a, const Point &b, const Point &c)
{
    return (Point::orientation_test(a, b, point_to_test) > 0
         && Point::orientation_test(b, c, point_to_test) > 0
         && Point::orientation_test(c, a, point_to_test) > 0) * 2 - 1;
}

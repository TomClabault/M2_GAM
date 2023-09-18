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


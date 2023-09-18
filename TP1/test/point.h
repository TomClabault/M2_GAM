#ifndef POINT_H
#define POINT_H

#include <ostream>

class Point
{
public:
    double x;
    double y;
    double z;

    Point():x(),y(),z() {}
    Point(float x_, float y_, float z_):x(x_),y(y_),z(z_) {}

    friend bool operator==(const Point& a, const Point& b);
    friend std::ostream& operator<<(std::ostream& os, const Point& point);
};

#endif // POINT_H

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

    friend bool operator==(const Point& a, const Point& b);
    friend std::ostream& operator<<(std::ostream& os, const Point& point);
};

#endif // POINT_H

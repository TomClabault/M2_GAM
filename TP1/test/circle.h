#ifndef CIRCLE_H
#define CIRCLE_H

#include "point.h"

class Circle
{
public:
    Circle(const Point& center, float radius);

    bool contains_point(const Point& point) const;

private:
    Point m_center;
    float m_radius;
};

#endif // CIRCLE_H

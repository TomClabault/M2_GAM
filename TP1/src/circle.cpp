#include "circle.h"

Circle::Circle(const Point& center, float radius) : m_center(center), m_radius(radius) {}

bool Circle::contains_point(const Point &point) const
{
    return length(m_center - point) < m_radius;
}

const Point &Circle::get_center() const
{
    return m_center;
}

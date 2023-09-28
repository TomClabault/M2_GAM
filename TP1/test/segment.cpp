#include "segment.h"

Segment::Segment(const Point& a, const Point& b) : a(a), b(b) { }

bool Segment::intersect(const Segment &other) const
{
    int orientation_first_test = Point::orientation_test(a, other.a, other.b) > 0 ? -1 : 1;
    int orientation_second_test = Point::orientation_test(b, other.a, other.b) > 0 ? -1 : 1;
    int orientation_third_test = Point::orientation_test(other.a, a, b) > 0 ? -1 : 1;
    int orientation_fourth_test = Point::orientation_test(other.b, a, b) > 0 ? -1 : 1;

    return (orientation_first_test != orientation_second_test)
        && (orientation_third_test != orientation_fourth_test);
}

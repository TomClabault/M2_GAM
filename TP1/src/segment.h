#ifndef SEGMENT_H
#define SEGMENT_H

#include "point.h"

class Segment
{
public:
    Segment(const Point& a, const Point& b);

    bool intersect(const Segment& other) const;

    friend bool operator==(const Segment& a, const Segment& b);

    Point a, b;
};

#endif // SEGMENT_H

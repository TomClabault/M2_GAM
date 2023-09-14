#ifndef POINT_H
#define POINT_H

class Point
{
public:
    double x;
    double y;
    double z;

    Point():x(),y(),z() {}
    Point(float x_, float y_, float z_):x(x_),y(y_),z(z_) {}
};

#endif // POINT_H

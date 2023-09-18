#ifndef GEOMETRICWORLD_H
#define GEOMETRICWORLD_H

#include "mesh.h"
#include "point.h"

class GeometricWorld //Here used to create a singleton instance
{
    QVector<Point> _bBox;  // Bounding box
public :
    void load_off(const char* filepath);

    void draw();
    void drawWireFrame();

    // ** TP Can be extended with further elements;
    Mesh _mesh;
};

#endif // GEOMETRICWORLD_H

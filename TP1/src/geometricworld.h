#ifndef GEOMETRICWORLD_H
#define GEOMETRICWORLD_H

#include "mesh.h"
#include "point.h"

class GeometricWorld //Here used to create a singleton instance
{
    QVector<Point> _bBox;  // Bounding box
public :
    void load_off(const char* filepath);
    void precompute_mesh_curvature();

    void draw();
    void drawWireFrame();

    // ** TP Can be extended with further elements;
    Mesh _mesh;
private:
    std::vector<Vector> m_precomputed_curvatures;
};

#endif // GEOMETRICWORLD_H

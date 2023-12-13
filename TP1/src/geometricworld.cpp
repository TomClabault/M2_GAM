#include "geometricworld.h"
#include "segment.h"
#include "TP1/offreader.h"

#include <iostream>

// Draw a Point
void glPointDraw(const Point & p) {
    glVertex3f(p.x, p.y, p.z);
}

void GeometricWorld::load_off(const char* filepath)
{
    // --------------- Demo terrain --------------- //
    /*

    //Sans Delaunay incrémental
    //_mesh.insert_point_cloud("../src/data/alpes_random_2.txt", 8, false);

    //Avec Delaunay incrémental
    _mesh.insert_point_cloud("../src/data/alpes_random_2.txt", 8, true);

    _mesh.scale_to_min_max_points(Point(-3, -3, -3), Point(3, 3, 3));

    */



    // ---------- Demo edge flip et face split ---------- //
    /*

    std::vector<Point> points;
    points.push_back(Point(0, 0, 0));
    points.push_back(Point(1, 0, 0));
    points.push_back(Point(2, 0.5, 0));
    points.push_back(Point(1.5, 1.5, 0));
    points.push_back(Point(0.5, 2, 0));
    points.push_back(Point(-0.5, 2, 0));
    points.push_back(Point(-0.5, 1, 0));

    for (const Point& point : points)
        _mesh.insert_point_2D(point, true);

    //_mesh.edge_flip(std::make_pair(1, 4));
    _mesh.face_split(0, Point(0.5, 0.5, 0));

    */





    // ---------- Demo Ruppert (cassé) ---------- //
    /*

    std::vector<Point> points;
    points.push_back(Point(0, 0, 0));
    points.push_back(Point(1, 0, 0));
    points.push_back(Point(2, 0.5, 0));
    points.push_back(Point(1.5, 1.5, 0));
    points.push_back(Point(0.5, 2, 0));
    points.push_back(Point(-0.5, 2, 0));
    points.push_back(Point(-0.5, 1, 0));

    std::vector<std::pair<int, int>> constraint_segments;
    for (int i = 0; i < points.size() - 1; i++)
        constraint_segments.push_back(std::make_pair(i, i + 1));
    constraint_segments.push_back(std::make_pair(0, points.size() - 1));

    _mesh.ruppert(constraint_segments, 20);

    */
}

void GeometricWorld::precompute_mesh_curvature()
{
    m_precomputed_curvatures.resize(_mesh.m_vertices.size());

    Vector max_curvature = Vector(std::numeric_limits<double>::min());
    for (int i = 0; i < _mesh.m_vertices.size(); i++)
    {
        Vector curvature = length(_mesh.laplacian_mean_curvature(i));
        max_curvature = max(max_curvature, curvature);

        m_precomputed_curvatures[i] = curvature;
    }

    for (Vector& curvature : m_precomputed_curvatures)
        curvature /= max_curvature;
}

#define COLOR_CURVATURE 0

//Example with a bBox
void GeometricWorld::draw()
{
    for (int i = 0; i < _mesh.m_faces.size(); i++)
    {
        const Face& face = _mesh.m_faces.at(i);
        Point& vertex_a = _mesh.m_vertices.at(face.m_a).get_point();
        Point& vertex_b = _mesh.m_vertices.at(face.m_b).get_point();
        Point& vertex_c = _mesh.m_vertices.at(face.m_c).get_point();

        glBegin(GL_TRIANGLES);
#if COLOR_CURVATURE
        Vector& curvature_a = m_precomputed_curvatures[face.m_a];
        glColor3d(curvature_a.x, curvature_a.y, curvature_a.z);
        glPointDraw(vertex_a);

        Vector& curvature_b = m_precomputed_curvatures[face.m_b];
        glColor3d(curvature_b.x, curvature_b.y, curvature_b.z);
        glPointDraw(vertex_b);

        Vector& curvature_c = m_precomputed_curvatures[face.m_c];
        glColor3d(curvature_c.x, curvature_c.y, curvature_c.z);
        glPointDraw(vertex_c);
#else
        glColor3d(1, 0, 0);
        glPointDraw(vertex_a);

        glColor3d(0, 1, 0);
        glPointDraw(vertex_b);

        glColor3d(0, 0, 1);
        glPointDraw(vertex_c);
#endif

        glEnd();
    }
}

//Example with a wireframe bBox
void GeometricWorld::drawWireFrame() {
    glColor3d(0,1,0);
    glBegin(GL_LINE_STRIP);
    glPointDraw(_bBox[0]);
    glPointDraw(_bBox[1]);
    glEnd();
    glColor3d(0,0,1);
    glBegin(GL_LINE_STRIP);
    glPointDraw(_bBox[0]);
    glPointDraw(_bBox[2]);
    glEnd();
    glColor3d(1,0,0);
    glBegin(GL_LINE_STRIP);
    glPointDraw(_bBox[0]);
    glPointDraw(_bBox[3]);
    glEnd();
}

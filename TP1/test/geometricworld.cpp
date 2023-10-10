#include "geometricworld.h"
#include "TP1/offreader.h"

#include <iostream>

// Draw a Point
void glPointDraw(const Point & p) {
    glVertex3f(p.x, p.y, p.z);
}

void GeometricWorld::load_off(const char* filepath)
{
    //_mesh = OffReader::read_off(filepath);
    //std::cout << std::filesystem::current_path() << std::endl;
    //_mesh = OffReader::read_off("cube_maillage_triangles.off");
    //_mesh = OffReader::read_off("queen.off");

    _mesh.add_vertex(Vertex(0, Point(0, 0, 0)));
    _mesh.add_vertex(Vertex(0, Point(1, 0, 0)));
    _mesh.add_vertex(Vertex(0, Point(0.5, 0.5, 0)));
    _mesh.add_vertex(Vertex(1, Point(1.5, 0.25, 0)));
    _mesh.add_vertex(Vertex(2, Point(-0.5, 0.25, 0)));

    _mesh.add_face(Face(0, 1, 2, 1, 2, -1));
    _mesh.add_face(Face(1, 3, 2, -1, 0, -1));
    _mesh.add_face(Face(0, 2, 4, -1, -1, 0));

    auto start = _mesh.incident_faces(0);
    auto end = _mesh.incident_faces_past_the_end();

    for (; start != end; ++start)
    {
        std::cout << *start << std::endl;
    }

    _mesh.insert_outside_convex_hull_2D(Point(0.5, -0.5, 0));
    _mesh.insert_outside_convex_hull_2D(Point(-0.5, -0.6, 0));
    _mesh.insert_outside_convex_hull_2D(Point(1.0, -0.6, 0));

    //_mesh.delaunayize_lawson();
    //precompute_mesh_curvature();
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

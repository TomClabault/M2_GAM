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
    _mesh = OffReader::read_off("cube_maillage_triangles.off");


    std::cout << Point::orientation_test(Point(0, 0, 0), Point(1, 0, 0), Point(0.5, 0.5, 0)) << std::endl;
    std::cout << Point::orientation_test(Point(0, 0, 0), Point(0.5, 0.5, 0), Point(1, 0, 0)) << std::endl;

    std::cout << Point::is_point_in_triangle(Point(0.25, 0.25, 0), Point(0, 0, 0), Point(1, 0, 0), Point(0.5, 0.5, 0)) << std::endl;

//    _mesh.add_vertex(Vertex(0, Point(0, 0, 0)));
//    _mesh.add_vertex(Vertex(0, Point(1, 0, 0)));
//    _mesh.add_vertex(Vertex(0, Point(0.5, 0.5, 0)));

//    _mesh.add_vertex(Vertex(1, Point(0.5, -0.5, 0)));
//    _mesh.add_vertex(Vertex(2, Point(1.5, 0.25, 0)));
//    _mesh.add_vertex(Vertex(3, Point(-0.5, 0.25, 0)));

//    _mesh.add_face(Face(0, 1, 2, 2, 3, 1));
//    _mesh.add_face(Face(0, 3, 1, -1, 0, -1));
//    _mesh.add_face(Face(1, 4, 2, -1, 0, -1));
//    _mesh.add_face(Face(2, 5, 0, -1, 0, -1));

    _mesh.laplacian_mean_curvature(0);
}

//Example with a bBox
void GeometricWorld::draw() {

    for (int i = 0; i < _mesh.m_faces.size(); i++)
    {
        const Face& face = _mesh.m_faces.at(i);
        Point vertex_a = _mesh.m_vertices.at(face.m_a).get_point();
        Point vertex_b = _mesh.m_vertices.at(face.m_b).get_point();
        Point vertex_c = _mesh.m_vertices.at(face.m_c).get_point();

        glBegin(GL_TRIANGLES);
        Vector curvature_a = (_mesh.laplacian_mean_curvature(face.m_a) + 1) * 0.5;
        glColor3d(curvature_a.x, curvature_a.y, curvature_a.z);
        //glColor3d(1, 0, 0);
        glPointDraw(vertex_a);

        //glColor3d(0, 1, 0);
        Vector curvature_b = (_mesh.laplacian_mean_curvature(face.m_b) + 1) * 0.5;
        glColor3d(curvature_b.x, curvature_b.y, curvature_b.z);
        glPointDraw(vertex_b);

        //glColor3d(0, 0, 1);
        Vector curvature_c = (_mesh.laplacian_mean_curvature(face.m_c) + 1) * 0.5;
        glColor3d(curvature_c.x, curvature_c.y, curvature_c.z);
        glPointDraw(vertex_c);

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

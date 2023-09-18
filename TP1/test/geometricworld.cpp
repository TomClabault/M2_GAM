#include "geometricworld.h"
#include "TP1/offreader.h"

#include <iostream>

// Draw a Point
void glPointDraw(const Point & p) {
    glVertex3f(p.x, p.y, p.z);
}

void GeometricWorld::load_off(const char* filepath)
{
    _mesh = OffReader::read_off(filepath);

    Mesh::Iterator_on_vertices its;
    Mesh::Circulator_on_faces cf;
    for (its = _mesh.vertices_begin();
         its != _mesh.vertices_past_the_end();
         ++its)
    {
        Mesh::Circulator_on_faces cfbegin
            = _mesh.incident_faces(*its) ;
        int cmpt = 0;
        for (cf=cfbegin, ++cf;
             cf!=cfbegin;
             cf++)
            cmpt++ ;
        std::cout<< "valence of the vertex" << cmpt << std::endl;
    }

//    Mesh::Circulator_on_faces ite = _mesh.incident_faces(0);
//    while (ite != _mesh.incident_faces_past_the_end())
//    {
//        std::cout << *ite << std::endl;
//        ite++;
//    }
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
        glColor3d(1, 0, 0);
        glPointDraw(vertex_a);
        glColor3d(0, 1, 0);
        glPointDraw(vertex_b);
        glColor3d(0, 0, 1);
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

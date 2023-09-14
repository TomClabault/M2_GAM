#include "gldisplaywidget.h"
#include "mesh.h"
#include "TP1/color.h"
#include "TP1/offreader.h"
#include "TP1/vector.h"

// The following functions could be displaced into a module OpenGLDisplayGeometricWorld that would include mesh.h

// Draw a Point
void glPointDraw(const Point & p) {
    glVertex3f(p.x, p.y, p.z);
}

void GeometricWorld::load_off(const char* filepath)
{
    _mesh = OffReader::read_off(filepath);
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
//    glColor3d(1,0,0);
//    glBegin(GL_TRIANGLES);
//    glPointDraw(_bBox[0]);
//    glPointDraw(_bBox[1]);
//    glPointDraw(_bBox[2]);
//    glEnd();

//    glColor3d(0,1,0);
//    glBegin(GL_TRIANGLES);
//    glPointDraw(_bBox[0]);
//    glPointDraw(_bBox[2]);
//    glPointDraw(_bBox[3]);
//    glEnd();

//    glColor3d(0,0,1);
//    glBegin(GL_TRIANGLES);
//    glPointDraw(_bBox[0]);
//    glPointDraw(_bBox[3]);
//    glPointDraw(_bBox[1]);
//    glEnd();

    //glColor3d(1,1,0);
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


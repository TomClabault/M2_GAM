#ifndef MESH_H
#define MESH_H

#include "point.h"

#include <QOpenGLWidget>
#include <vector>

class Face
{
public:
    /**
     * @brief Face
     * @param a Index of the first vertex of the face
     * @param b Index of the second vertex of the face
     * @param fa Index of the face that is opposing the vertex "a"
     * @param fb Index of the face that is opposing the vertex "b"
     * @param fc Index of the face that is opposing the vertex "c"
     */
    Face(int a, int b, int c, int fa, int fb, int fc) : m_a(a), m_b(b), m_c(c), m_fa(fa), m_fb(fb), m_fc(fc) { }

    int& m_opposing_faces(int index) {return *(&m_fa + index);}
    int& m_vertices(int index) {return *(&m_a + index);}

    //Vertices of the face
    int m_a, m_b, m_c;

    //Adjacent faces. The faces are opposing the corresponding vertex:
    //fa is the face opposing the vertex 'a' (so it's the face that shares the 'bc' line)
    //...
    int m_fa, m_fb, m_fc;
};

class Vertex
{
public:
    Vertex(int adjacent_face_index, Point geometric_coordinates) : m_adjacent_face_index(adjacent_face_index), m_geometric_point(geometric_coordinates) {}

    int get_adjacent_face_index() { return m_adjacent_face_index; }
    void set_adjacent_face_index(int adjacent_face_index) { m_adjacent_face_index = adjacent_face_index; }

    Point& get_point() { return m_geometric_point; }

private:
    //Index of a face that contains the current vertex
    int m_adjacent_face_index;

    Point m_geometric_point;
};

//** TP : TO MODIFY

class Mesh
{
public:
    Mesh() {}
    Mesh(std::vector<Face>& faces, std::vector<Vertex>& vertices) : m_faces(faces), m_vertices(vertices) {}

    void add_vertex(const Vertex& vertex)
    {
        m_vertices.push_back(vertex);
    }

    std::vector<Face> m_faces;
    std::vector<Vertex> m_vertices;
    //void drawMesh();
    //void drawMeshWireFrame();
};

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


#endif // MESH_H

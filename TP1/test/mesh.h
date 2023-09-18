#ifndef MESH_H
#define MESH_H

#include "face.h"
#include "vertex.h"

#include <QOpenGLWidget>
#include <vector>

class Mesh
{
public:
    struct Iterator_on_faces
    {
    public:
        Iterator_on_faces() { m_mesh = nullptr; };
        Iterator_on_faces(Mesh& mesh, bool past_the_end);
        Iterator_on_faces(Mesh& mesh) : m_mesh(&mesh), m_current_face_index(0) , m_past_the_end(false) {}

        Face& operator*();

        friend Mesh::Iterator_on_faces& operator++(Mesh::Iterator_on_faces& operand);
        friend Mesh::Iterator_on_faces operator++(Mesh::Iterator_on_faces& operand, int dummy);

        friend bool operator ==(const Mesh::Iterator_on_faces& a, const Mesh::Iterator_on_faces& b);
        friend bool operator !=(const Mesh::Iterator_on_faces& a, const Mesh::Iterator_on_faces& b);

    private:

        Mesh* m_mesh;
        unsigned long long int m_current_face_index;

        bool m_past_the_end;
    };

    struct Iterator_on_vertices
    {
    public:
        Iterator_on_vertices() { m_mesh = nullptr; };
        Iterator_on_vertices(Mesh& mesh, bool past_the_end);
        Iterator_on_vertices(Mesh& mesh) : m_mesh(&mesh), m_current_vertex_index(0) , m_past_the_end(false) {}

        Vertex& operator*();

        friend Mesh::Iterator_on_vertices& operator++(Mesh::Iterator_on_vertices& operand);
        friend Mesh::Iterator_on_vertices operator++(Mesh::Iterator_on_vertices& operand, int dummy);

        friend bool operator ==(const Mesh::Iterator_on_vertices& a, const Mesh::Iterator_on_vertices& b);
        friend bool operator !=(const Mesh::Iterator_on_vertices& a, const Mesh::Iterator_on_vertices& b);

    private:
        Mesh* m_mesh;
        unsigned long long int m_current_vertex_index;

        bool m_past_the_end;
    };

    struct Circulator_on_faces
    {
    public:
        Circulator_on_faces() { m_mesh = nullptr; };
        Circulator_on_faces(Mesh& mesh, int current_face_index) : m_mesh(&mesh), m_current_face_index(current_face_index) {}
        Circulator_on_faces(Mesh& mesh, Vertex& vertex) : m_mesh(&mesh), m_vertex_circulating_around(&vertex), m_current_face_index(vertex.get_adjacent_face_index()) {}

        Face& operator*();

        friend Mesh::Circulator_on_faces& operator++(Mesh::Circulator_on_faces& operand);
        friend Mesh::Circulator_on_faces operator++(Mesh::Circulator_on_faces& operand, int dummy);

        friend bool operator ==(const Mesh::Circulator_on_faces& a, const Mesh::Circulator_on_faces& b);
        friend bool operator !=(const Mesh::Circulator_on_faces& a, const Mesh::Circulator_on_faces& b);

    private:
        Mesh* m_mesh;
        Vertex* m_vertex_circulating_around;

        int m_current_face_index;
    };

    Mesh() {}
    Mesh(std::vector<Face>& faces, std::vector<Vertex>& vertices) : m_faces(faces), m_vertices(vertices) {}

    void add_vertex(const Vertex& vertex);

    Iterator_on_faces faces_begin();
    Iterator_on_faces faces_past_the_end();

    Iterator_on_vertices vertices_begin();
    Iterator_on_vertices vertices_past_the_end();

    Circulator_on_faces incident_faces(int vertex_index);
    Circulator_on_faces incident_faces(Vertex& vertex);
    Circulator_on_faces incident_faces_past_the_end();

    std::vector<Face> m_faces;
    std::vector<Vertex> m_vertices;
    //void drawMesh();
    //void drawMeshWireFrame();
};

#endif // MESH_H

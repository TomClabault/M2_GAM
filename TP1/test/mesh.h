#ifndef MESH_H
#define MESH_H

#include "circle.h"
#include "face.h"
#include "vertex.h"

#include <QOpenGLWidget>
#include <unordered_set>
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

    struct Iterator_on_edges
    {
    public:
        Iterator_on_edges() { m_mesh = nullptr; }
        Iterator_on_edges(Mesh& mesh, bool past_the_end) { m_mesh = &mesh; m_past_the_end = past_the_end; };
        Iterator_on_edges(Mesh& mesh);

        std::pair<int, int> operator*();

        friend Mesh::Iterator_on_edges& operator++(Mesh::Iterator_on_edges& operand);
        friend Mesh::Iterator_on_edges operator++(Mesh::Iterator_on_edges& operand, int dummy);

        friend bool operator ==(const Mesh::Iterator_on_edges& a, const Mesh::Iterator_on_edges& b);
        friend bool operator !=(const Mesh::Iterator_on_edges& a, const Mesh::Iterator_on_edges& b);

        int get_current_face_index() { return m_current_face_index; }
        int get_opposite_face_index();

    private:
        //From: https://stackoverflow.com/questions/15160889/how-can-i-make-an-unordered-set-of-pairs-of-integers-in-c
        struct pair_hash {
            inline std::size_t operator()(const std::pair<int,int> & v) const {
                return v.first*31+v.second;
            }
        };

        Mesh* m_mesh;
        std::pair<int, int> m_current_edge;
        int m_current_face_index;
        int m_current_edge_in_current_face;

        bool m_past_the_end = false;

        std::unordered_set<std::pair<int, int>, pair_hash> already_visited_edges;
    };

    struct Circulator_on_faces
    {
    public:
        Circulator_on_faces() { m_mesh = nullptr; };
        Circulator_on_faces(Mesh& mesh, int current_face_index) : m_mesh(&mesh), m_current_face_index(current_face_index), m_start_face_index(current_face_index) {}
        Circulator_on_faces(Mesh& mesh, Vertex& vertex) : m_mesh(&mesh), m_vertex_circulating_around(&vertex), m_current_face_index(vertex.get_adjacent_face_index()), m_start_face_index(vertex.get_adjacent_face_index()) {}

        Face& operator*();
        int get_current_face_index() { return m_current_face_index; }

        friend Mesh::Circulator_on_faces& operator++(Mesh::Circulator_on_faces& operand);
        friend Mesh::Circulator_on_faces operator++(Mesh::Circulator_on_faces& operand, int dummy);

        friend bool operator ==(const Mesh::Circulator_on_faces& a, const Mesh::Circulator_on_faces& b);
        friend bool operator !=(const Mesh::Circulator_on_faces& a, const Mesh::Circulator_on_faces& b);

    private:
        Mesh* m_mesh;
        Vertex* m_vertex_circulating_around;

        int m_current_face_index;
        int m_start_face_index;

        bool m_circulating_backwards = false;
    };

    Mesh() {}
    Mesh(std::vector<Face>& faces, std::vector<Vertex>& vertices) : m_faces(faces), m_vertices(vertices) {}

    void add_vertex(const Vertex& vertex);
    void add_face(const Face &face);
    void push_convex_hull_edge(int index_vertex1, int index_vertex2);
    void push_convex_hull_edge_face(int face_index);

    Iterator_on_faces faces_begin();
    Iterator_on_faces faces_past_the_end();

    Iterator_on_vertices vertices_begin();
    Iterator_on_vertices vertices_past_the_end();

    Iterator_on_edges edges_begin();
    Iterator_on_edges edges_past_the_end();

    Circulator_on_faces incident_faces(int vertex_index);
    Circulator_on_faces incident_faces(Vertex& vertex);
    Circulator_on_faces incident_faces_past_the_end();

    double face_area(const Face& face);
    Point barycenter_of_face(const Face& face) const;
    Circle get_circumscribed_circle_of_face(const Face& face) const;

    Vector laplacian_mean_curvature(const int vertex_index);

    std::pair<int, int> get_faces_indices_of_edge(std::pair<int, int> two_vertices_indices);
    void face_split(const int face_index, const Point& new_point);
    void edge_flip(const std::pair<int, int>& vertex_index_pair);
    void edge_flip(const int face_index_1, const int face_index_2);
    void insert_point_2D(const Point& point);
    bool is_edge_locally_delaunay(int face1_index, int face2_index);
    bool is_edge_locally_delaunay(const std::pair<int, int> two_vertex_indices);

    void delaunayize_lawson();

    void compute_convex_hull_edges();
    void insert_outside_convex_hull_2D(const Point& point);

    //Edges that are on the convex hull of our triangulation
    std::list<std::pair<int, int>> m_convex_hull_edges;
    //Face asociated with the edge of 'm_convex_hull_edges'
    std::list<int> m_convex_hull_edges_faces;
    std::vector<Face> m_faces;
    std::vector<Vertex> m_vertices;
};

#endif // MESH_H

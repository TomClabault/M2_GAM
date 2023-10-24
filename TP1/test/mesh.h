#ifndef MESH_H
#define MESH_H

#include "circle.h"
#include "face.h"
#include "segment.h"
#include "vertex.h"

#include <QOpenGLWidget>
#include <unordered_set>
#include <vector>

class Mesh
{
public:
    //From: https://stackoverflow.com/questions/15160889/how-can-i-make-an-unordered-set-of-pairs-of-integers-in-c
    struct pair_hash {
        inline std::size_t operator()(const std::pair<int,int> & v) const {
            return v.first*31+v.second;
        }
    };

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

    struct Iterator_on_convex_hull_edges_in_order
    {
        struct Convex_hull_edge
        {
            int face_index;
            int local_index_of_vertex_opposite_to_edge;

            friend bool operator==(const Convex_hull_edge& a, const Convex_hull_edge& b);
        };

        Iterator_on_convex_hull_edges_in_order(Mesh& mesh);
        Iterator_on_convex_hull_edges_in_order(Mesh& mesh, Convex_hull_edge start_edge) : m_mesh(&mesh), m_current_edge(start_edge), m_start_edge(start_edge) {}
        Iterator_on_convex_hull_edges_in_order(bool past_the_end) : m_past_the_end(past_the_end) {}

        const std::pair<int, int> operator*();

        friend Mesh::Iterator_on_convex_hull_edges_in_order& operator++(Mesh::Iterator_on_convex_hull_edges_in_order& operand);
        friend Mesh::Iterator_on_convex_hull_edges_in_order operator++(Mesh::Iterator_on_convex_hull_edges_in_order& operand, int dummy);

        const Convex_hull_edge& get_current_edge();

        friend bool operator ==(const Mesh::Iterator_on_convex_hull_edges_in_order& a, const Mesh::Iterator_on_convex_hull_edges_in_order& b);
        friend bool operator !=(const Mesh::Iterator_on_convex_hull_edges_in_order& a, const Mesh::Iterator_on_convex_hull_edges_in_order& b);

        Mesh* m_mesh;

        Convex_hull_edge m_current_edge;
        Convex_hull_edge m_start_edge;

        bool m_past_the_end = false;
    };

    Mesh() {}
    Mesh(std::vector<Face>& faces, std::vector<Vertex>& vertices) : m_faces(faces), m_vertices(vertices) {}

    void add_vertex(const Vertex& vertex);
    void add_face(const Face &face);
    int get_one_convex_hull_face();

    Iterator_on_faces faces_begin();
    Iterator_on_faces faces_past_the_end();

    Iterator_on_vertices vertices_begin();
    Iterator_on_vertices vertices_past_the_end();

    Iterator_on_edges edges_begin();
    Iterator_on_edges edges_past_the_end();

    Circulator_on_faces incident_faces(int vertex_index);
    Circulator_on_faces incident_faces(Vertex& vertex);
    Circulator_on_faces incident_faces_past_the_end();

    Iterator_on_convex_hull_edges_in_order convex_hull_edges_in_order_begin();
    Iterator_on_convex_hull_edges_in_order convex_hull_edges_in_order_begin(Mesh::Iterator_on_convex_hull_edges_in_order::Convex_hull_edge& current_edge);
    Iterator_on_convex_hull_edges_in_order convex_hull_edges_in_order_past_the_end();

    void save_as_obj(const std::string& filepath) const;
    void save_as_off(const std::string& filepath) const;
    void insert_point_cloud(const std::string &filepath);

    double face_area(const Face& face);
    Point barycenter_of_face(int face_index) const;
    Point barycenter_of_face(const Face& face) const;
    std::vector<std::pair<int, int>> get_edges_of_face(const Face& face);
    Circle get_circumscribed_circle_of_face(const Face& face) const;
    std::pair<int, int> get_faces_indices_of_edge(std::pair<int, int> two_vertices_indices);

    Vector laplacian_mean_curvature(const int vertex_index);

    void edge_split(const std::pair<int, int>& two_vertices);
    void face_split(const int face_index, const Point& new_point, bool apply_lawson_after_insertion = true);
    void check_and_flip_multiple_edges(std::vector<std::pair<int, int>>& edges_to_flip);
    std::vector<std::pair<int, int>> edge_flip(const std::pair<int, int>& vertex_index_pair);
    std::vector<std::pair<int, int>> edge_flip(const int face_index_1, const int face_index_2);
    void insert_point_2D(const Point& point, bool apply_lawson_after_insertion);
    bool point_inside_triangulation_brute_force(const Point& point);
    bool is_edge_locally_delaunay(int face1_index, int face2_index);
    bool is_edge_locally_delaunay(const std::pair<int, int> two_vertex_indices);

    void delaunayize_lawson();
    void ruppert(const std::vector<Segment> &constraint_segments);

    void compute_convex_hull_edges();
    void insert_outside_convex_hull_2D(const Point& point, bool apply_lawson_after_insertion = true);

    void scale_to_min_max_points(const Point& min, const Point& max);

    std::vector<Face> m_faces;
    std::vector<Vertex> m_vertices;
};

#endif // MESH_H

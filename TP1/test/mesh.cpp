#include "mesh.h"
#include "segment.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <queue>

Mesh::Iterator_on_faces::Iterator_on_faces(Mesh& mesh, bool past_the_end)
{
    m_mesh = &mesh;
    m_past_the_end = past_the_end;
}

Face& Mesh::Iterator_on_faces::operator*()
{
    return m_mesh->m_faces.at(m_current_face_index);
}

//Prefix
Mesh::Iterator_on_faces& operator++(Mesh::Iterator_on_faces& operand)
{
    if (operand.m_current_face_index >= (operand.m_mesh->m_faces.size() - 1))
        operand.m_past_the_end = true;
    else
        operand.m_current_face_index++;

    return operand;
}

//Postfix
Mesh::Iterator_on_faces operator++(Mesh::Iterator_on_faces& operand, int dummy)
{
    Mesh::Iterator_on_faces copy = operand;

    ++operand;

    return copy;
}

bool operator ==(const Mesh::Iterator_on_faces& a, const Mesh::Iterator_on_faces& b)
{
    if (a.m_past_the_end != b.m_past_the_end)
        return false;
    else if (a.m_past_the_end == b.m_past_the_end && a.m_mesh == b.m_mesh)
        return true;

    if (a.m_current_face_index != b.m_current_face_index)
        return false;

    if (a.m_mesh != b.m_mesh)
        return false;

    return true;
}

bool operator !=(const Mesh::Iterator_on_faces& a, const Mesh::Iterator_on_faces& b)
{
    return !(a == b);
}

Mesh::Iterator_on_faces Mesh::faces_begin()
{
    return Mesh::Iterator_on_faces(*this);
}

Mesh::Iterator_on_faces Mesh::faces_past_the_end()
{
    return Mesh::Iterator_on_faces(*this, true);
}











Mesh::Iterator_on_vertices::Iterator_on_vertices(Mesh& mesh, bool past_the_end)
{
    m_mesh = &mesh;
    m_past_the_end = past_the_end;
}

Vertex& Mesh::Iterator_on_vertices::operator*()
{
    return m_mesh->m_vertices.at(m_current_vertex_index);
}

//Prefix
Mesh::Iterator_on_vertices& operator++(Mesh::Iterator_on_vertices& operand)
{
    if (operand.m_current_vertex_index >= (operand.m_mesh->m_vertices.size() - 1))
        operand.m_past_the_end = true;
    else
        operand.m_current_vertex_index++;

    return operand;
}

//Postfix
Mesh::Iterator_on_vertices operator++(Mesh::Iterator_on_vertices& operand, int dummy)
{
    Mesh::Iterator_on_vertices copy = operand;

    ++operand;

    return copy;
}

bool operator ==(const Mesh::Iterator_on_vertices& a, const Mesh::Iterator_on_vertices& b)
{
    if (a.m_past_the_end != b.m_past_the_end)
        return false;
    else if (a.m_past_the_end == b.m_past_the_end && a.m_mesh == b.m_mesh)
        return true;

    if (a.m_current_vertex_index != b.m_current_vertex_index)
        return false;

    if (a.m_mesh != b.m_mesh)
        return false;

    return true;
}

bool operator !=(const Mesh::Iterator_on_vertices& a, const Mesh::Iterator_on_vertices& b)
{
    return !(a == b);
}

Mesh::Iterator_on_edges::Iterator_on_edges(Mesh& mesh) : m_mesh(&mesh), m_current_face_index(0), m_current_edge_in_current_face(0)
{
    Face& current_face = m_mesh->m_faces[m_current_face_index];

    m_current_edge = std::make_pair(current_face.m_a, current_face.m_b);
}

std::pair<int, int> Mesh::Iterator_on_edges::operator*()
{
    return m_current_edge;
}

int Mesh::Iterator_on_edges::get_opposite_face_index()
{
    if (!m_past_the_end)
    {
        return m_mesh->m_faces[m_current_face_index].opposing_face((m_current_edge_in_current_face + 2) % 3);
    }
}

Mesh::Iterator_on_edges& operator++(Mesh::Iterator_on_edges& operand)
{
    if (!operand.m_past_the_end)
    {
        Mesh* mesh = operand.m_mesh;

        Face& current_face = mesh->m_faces[operand.m_current_face_index];

        int edge_point_1_global_index = current_face.global_index_of_local_vertex_index(operand.m_current_edge_in_current_face);
        int edge_point_2_global_index = current_face.global_index_of_local_vertex_index((operand.m_current_edge_in_current_face + 1) % 3);

        bool edge_p1_smallest = edge_point_1_global_index < edge_point_2_global_index;

        std::pair<int, int> new_edge = std::make_pair(edge_p1_smallest ? edge_point_1_global_index : edge_point_2_global_index,
                                                      edge_p1_smallest ? edge_point_2_global_index : edge_point_1_global_index);

        bool edge_already_visited = operand.already_visited_edges.find(new_edge) == operand.already_visited_edges.end();
        bool need_to_increment = false;
        if (!edge_already_visited)
        {
            operand.m_current_edge = new_edge;
            operand.already_visited_edges.insert(operand.m_current_edge);
        }
        else
            //We're going to have to skip to the next edge because we already visited this edge
            need_to_increment = true;


        if (operand.m_current_edge_in_current_face == 2)
        {
            operand.m_current_edge_in_current_face = 0;
            operand.m_current_face_index++;
            if ((unsigned long long int)operand.m_current_face_index == mesh->m_faces.size())
                operand.m_past_the_end = true;
        }
        else
            operand.m_current_edge_in_current_face++;

        if (need_to_increment)
            ++operand;
    }

    return operand;
}

Mesh::Iterator_on_edges operator++(Mesh::Iterator_on_edges& operand, int dummy)
{
    Mesh::Iterator_on_edges copy = operand;

    ++operand;

    return copy;
}

bool operator ==(const Mesh::Iterator_on_edges& a, const Mesh::Iterator_on_edges& b)
{
    return (a.m_past_the_end == b.m_past_the_end) || ((a.m_current_edge == b.m_current_edge) && (a.m_current_face_index == b.m_current_face_index));
}

bool operator !=(const Mesh::Iterator_on_edges& a, const Mesh::Iterator_on_edges& b)
{
    return !(a == b);
}

Mesh::Iterator_on_vertices Mesh::vertices_begin()
{
    return Mesh::Iterator_on_vertices(*this);
}

Mesh::Iterator_on_vertices Mesh::vertices_past_the_end()
{
    return Mesh::Iterator_on_vertices(*this, true);
}

Mesh::Iterator_on_edges Mesh::edges_begin()
{
    return Iterator_on_edges(*this);
}

Mesh::Iterator_on_edges Mesh::edges_past_the_end()
{
    return Iterator_on_edges(*this, true);
}

Mesh::Circulator_on_faces Mesh::incident_faces(int vertex_index)
{
    return Mesh::Circulator_on_faces(*this, m_vertices[vertex_index]);
}

Mesh::Circulator_on_faces Mesh::incident_faces(Vertex& vertex)
{
    return Mesh::Circulator_on_faces(*this, vertex);
}

Mesh::Circulator_on_faces Mesh::incident_faces_past_the_end()
{
    return Mesh::Circulator_on_faces(*this, -1);
}

Mesh::Iterator_on_convex_hull_edges_in_order Mesh::convex_hull_edges_in_order_begin()
{
    return Mesh::Iterator_on_convex_hull_edges_in_order(*this);
}

Mesh::Iterator_on_convex_hull_edges_in_order Mesh::convex_hull_edges_in_order_begin(Mesh::Iterator_on_convex_hull_edges_in_order::Convex_hull_edge& start_edge)
{
    return Mesh::Iterator_on_convex_hull_edges_in_order(*this, start_edge);
}

Mesh::Iterator_on_convex_hull_edges_in_order Mesh::convex_hull_edges_in_order_past_the_end()
{
    return Mesh::Iterator_on_convex_hull_edges_in_order(true);
}

void Mesh::save_as_obj(const std::string& filepath) const
{
    std::ofstream output_file(filepath);

    for (const Vertex& vertex : m_vertices)
    {
        Point p = vertex.get_point();
        output_file << "v " << p.x / 1000 << " " << p.y / 1000 << " " << p.z / 1000 << std::endl;
    }

    for (const Face& face : m_faces)
        output_file << "f " << face.m_a + 1 << " " << face.m_b + 1 << " " << face.m_c + 1 << std::endl;
}

void Mesh::save_as_off(const std::string& filepath) const
{
    std::ofstream output_file(filepath);
    output_file << "OFF" << std::endl;
    output_file << m_vertices.size() << " " << m_faces.size() << " 0" << std::endl;

    for (const Vertex& vertex : m_vertices)
    {
        Point p = vertex.get_point();
        output_file << p.x / 100 << " " << p.y / 100 << " " << p.z / 100 << std::endl;
    }

    for (const Face& face : m_faces)
        output_file << "3 " << face.m_a << " " << face.m_b << " " << face.m_c << std::endl;
}

void Mesh::insert_point_cloud(const std::string& filepath)
{
    std::ifstream input_file(filepath);
    if (!input_file.is_open())
    {
        std::cerr << "Error while opening the point cloud file " << filepath << std::endl;
        return;
    }

    int nb_vertices;

    input_file >> nb_vertices;

    //Reading all the vertices of the file
    std::vector<double> z_coordinates(nb_vertices);
    for (int i = 0; i < nb_vertices; i++)
    {
        double x, y, z;

        input_file >> x;
        input_file >> y;
        input_file >> z;
        z_coordinates[i] = z;
        z = 0.0;//We don't want the z coordinate for our usecase

        Point point(x, y, z);
        insert_point_2D(point, true);
    }

    //Restoring the z coordinates
    for (int i = 0; i < m_vertices.size(); i++)
        m_vertices[i].get_point().z = z_coordinates[i];
}

double Mesh::face_area(const Face &face)
{
    Vector face_ab = m_vertices[face.m_b] - m_vertices[face.m_a];
    Vector face_ac = m_vertices[face.m_c] - m_vertices[face.m_a];

    return length(cross(face_ab, face_ac)) / 2;
}

Vector Mesh::laplacian_mean_curvature(const int vertex_index)
{
    Vertex& vertex = m_vertices[vertex_index];
    const Point& vertex_point = vertex.get_point();
    Mesh::Circulator_on_faces circulator = incident_faces(vertex);
    Mesh::Circulator_on_faces circulator_end =  incident_faces_past_the_end();

    Vector laplacien_sum(0, 0, 0);
    double neighbor_faces_area_sum = 0.0;
    do
    {
        Face& current_circulator_face = *circulator;
        neighbor_faces_area_sum += this->face_area(current_circulator_face) / 3;

        int local_index_in_face_of_vertex_index = current_circulator_face.local_index_of_global_vertex_index(vertex_index);
        int neighbor_point_global_index = current_circulator_face.global_index_of_local_vertex_index((local_index_in_face_of_vertex_index + 2) % 3);
        Vertex& neighbor_vertex = m_vertices[neighbor_point_global_index];
        Point& neighbor_point = neighbor_vertex.get_point();

        Face& neighbor_face2 = m_faces[current_circulator_face.opposing_face((local_index_in_face_of_vertex_index + 1) % 3)];
        int alpha_point_global_index = current_circulator_face.global_index_of_local_vertex_index((local_index_in_face_of_vertex_index + 1) % 3);
        int beta_point_global_index = neighbor_face2.global_index_of_local_vertex_index((neighbor_face2.local_index_of_global_vertex_index(vertex_index) + 2) % 3);

        Point alpha_point = m_vertices[alpha_point_global_index].get_point();
        Point beta_point = m_vertices[beta_point_global_index].get_point();

        Vector alpha_edge_1 = normalize(neighbor_point - alpha_point);
        Vector alpha_edge_2 = normalize(vertex_point - alpha_point);
        double alpha_cos = dot(alpha_edge_1, alpha_edge_2);
        //Vector alpha_cross = cross(alpha_edge_1, alpha_edge_2);
        //Vector alpha_face_normal = cross(alpha_edge_1, alpha_edge_2);
        //float alpha_cross_sign = dot(alpha_cross, alpha_face_normal);
        //double alpha_sin = length(alpha_cross) * alpha_cross_sign;
        //double alpha_sin = length(alpha_cross);
        double alpha_sin = length(cross(alpha_edge_1, alpha_edge_2));
        double cotangent_alpha = alpha_cos / alpha_sin;

        Vector beta_edge_1 = normalize(vertex_point - beta_point);
        Vector beta_edge_2 = normalize(neighbor_point - beta_point);
        double beta_cos = dot(beta_edge_1, beta_edge_2);
        //Vector beta_cross = cross(beta_edge_1, beta_edge_2);
        //Vector beta_face_normal = cross()
        //float beta_sign = dot(beta_cross, beta_face_normal);
        //double beta_sin = length(beta_cross) * beta_sign;
        double beta_sin = length(cross(beta_edge_1, beta_edge_2));
        double cotangent_beta = beta_cos / beta_sin;

        Vector difference = neighbor_point - vertex_point;

        Vector laplacien = (cotangent_alpha + cotangent_beta) * difference;

        laplacien_sum += laplacien;

        circulator++;
    } while (circulator != circulator_end);

    return laplacien_sum / (2 * neighbor_faces_area_sum);
}

void Mesh::edge_split(const std::pair<int, int>& two_vertices)
{
    std::pair<int, int> faces = get_faces_indices_of_edge(two_vertices);

    int face0_index = faces.first;
    int face1_index = faces.second;
    Face& face0 = m_faces[face0_index];
    Face& face1 = m_faces[face1_index];
    const Face& face0_copy = m_faces[face0_index];
    const Face& face1_copy = m_faces[face1_index];

    int global_vertex_on_face0_opposing_to_face1 = face0_copy.find_global_vertex_index_with_opposing_face(face1_index);
    int global_vertex_on_face1_opposing_to_face0 = face1_copy.find_global_vertex_index_with_opposing_face(face0_index);
    int local_vertex_on_face0_opposing_to_face_1 = face0.local_index_of_global_vertex_index(global_vertex_on_face0_opposing_to_face1);
    int local_vertex_on_face1_opposing_to_face_0 = face1.local_index_of_global_vertex_index(global_vertex_on_face1_opposing_to_face0);

    //Changing the vertices of the already existing faces
    int new_vertex_index = m_vertices.size();
    face0.m_a = new_vertex_index;
    face1.m_a = new_vertex_index;

    const Point& first_edge_point = m_vertices[two_vertices.first].get_point();
    const Point& second_edge_point = m_vertices[two_vertices.second].get_point();
    Vertex new_vertex(face0_index, (first_edge_point + second_edge_point) / 2);

    //Creating the 2 new faces
    int new_face_2_index = m_faces.size();
    int new_face_3_index = new_face_2_index + 1;
    Face new_face_2 = Face(two_vertices.first, global_vertex_on_face1_opposing_to_face0, new_vertex_index,
                           face1_index, new_face_3_index, face1_copy.opposing_face((local_vertex_on_face1_opposing_to_face_0 + 1) % 3));
    Face new_face_3 = Face(two_vertices.first, new_vertex_index, global_vertex_on_face0_opposing_to_face1,
                           face0_index, face0_copy.opposing_face((local_vertex_on_face0_opposing_to_face_1 + 2) % 3), new_face_2_index);

    //Update opposing faces
    face0.opposing_face(0) = face0_copy.opposing_face(0);
    face0.opposing_face(1) = new_face_3_index;
    face0.opposing_face(2) = face1_index;

    face1.opposing_face(0) = face1_copy.opposing_face(0);
    face1.opposing_face(1) = face0_index;
    face1.opposing_face(2) = new_face_2_index;

    //Update adjacencies
    m_vertices[two_vertices.first].set_adjacent_face_index(new_face_2_index);
    m_vertices[two_vertices.second].set_adjacent_face_index(face0_index);

    m_vertices.push_back(new_vertex);
    m_faces.push_back(new_face_2);
    m_faces.push_back(new_face_3);
}

std::pair<int, int> Mesh::get_faces_indices_of_edge(std::pair<int, int> vertex_index_pair)
{
    int first_face_index = -1, second_face_index = -1;

    auto faces_circulator = incident_faces(vertex_index_pair.first);
    auto faces_circulator_end = incident_faces_past_the_end();

    for (; faces_circulator != faces_circulator_end; faces_circulator++)
    {
        Face& current_face = *faces_circulator;
            if (current_face.contains_global_vertex_index(vertex_index_pair.first)
         && current_face.contains_global_vertex_index(vertex_index_pair.second))
        {
            first_face_index = faces_circulator.get_current_face_index();

            int local_index_of_first_vertex = current_face.local_index_of_global_vertex_index(vertex_index_pair.first);
            int local_index_of_second_vertex = current_face.local_index_of_global_vertex_index(vertex_index_pair.second);

            second_face_index = current_face.opposing_face(current_face.get_local_vertex_index_opposing_to_edge(local_index_of_first_vertex, local_index_of_second_vertex));

            break;
        }
    }

    return std::make_pair(first_face_index, second_face_index);
}

Point Mesh::barycenter_of_face(int face_index) const
{
    return barycenter_of_face(m_faces[face_index]);
}

Point Mesh::barycenter_of_face(const Face &face) const
{
    return (m_vertices[face.m_a].get_point() + m_vertices[face.m_b].get_point() + m_vertices[face.m_c].get_point()) / 3;
}

std::vector<std::pair<int, int>> Mesh::get_edges_of_face(const Face& face)
{
    return std::vector<std::pair<int, int>> {std::make_pair(face.m_a, face.m_b),
                                             std::make_pair(face.m_b, face.m_c),
                                             std::make_pair(face.m_c, face.m_a)};
}

Circle Mesh::get_circumscribed_circle_of_face(const Face &face) const
{
    Point point_a, point_b, point_c;
    point_a = m_vertices[face.m_a].get_point();
    point_b = m_vertices[face.m_b].get_point();
    point_c = m_vertices[face.m_c].get_point();

    Vector AB = point_b - point_a;
    Vector BC = point_c - point_b;
    Vector CA = point_a - point_c;

    double sin_A, sin_B, sin_C;

    sin_A = cross(AB, -CA).z;
    sin_B = cross(BC, -AB).z;
    sin_C = cross(CA, -BC).z;

    double dot_A = dot(AB, -CA);
    double dot_B = dot(BC, -AB);
    double dot_C = dot(CA, -BC);
    double prod_A = dot_B * dot_C;
    double prod_B = dot_A * dot_C;
    double prod_C = dot_A * dot_B;

    Point barycentric_coords = Point(sin_B * prod_B + sin_C * prod_C, sin_C * prod_C + sin_A * prod_A, sin_A * prod_A + sin_B * prod_B);
    barycentric_coords /= (barycentric_coords[0] + barycentric_coords[1] + barycentric_coords[2]);

    Point circle_center = Point(barycentric_coords[0] * point_a + barycentric_coords[1] * point_b + barycentric_coords[2] * point_c);
    return Circle(circle_center, length(circle_center - point_a));
}

void Mesh::face_split(const int face_index, const Point& new_point, bool apply_lawson_after_insertion)
{
    //TODO test que le apply lawson after insertion marche bien
    Face face_0_copy = m_faces[face_index];
    Face& new_face0 = m_faces[face_index];
    Face new_face1, new_face2;

    int new_face1_index = m_faces.size();
    int new_face2_index = m_faces.size() + 1;

    int face_touching_new_face_0 = face_0_copy.opposing_face(1);
    int face_touching_new_face_1 = face_0_copy.opposing_face(2);
    int face_touching_new_face_2 = face_0_copy.opposing_face(0);

    m_vertices.push_back(Vertex(face_index, new_point));

    //The current face we're inserting into is going to become 'face 0' of the
    //3 new faces created
    //The new point is always the vertex 0 of the new faces
    new_face0.m_a = m_vertices.size() - 1;
    new_face1.m_a = m_vertices.size() - 1;
    new_face2.m_a = m_vertices.size() - 1;

    new_face0.m_b = face_0_copy.m_c;
    new_face1.m_b = face_0_copy.m_a;
    new_face2.m_b = face_0_copy.m_b;

    new_face0.m_c = face_0_copy.m_a;
    new_face1.m_c = face_0_copy.m_b;
    new_face2.m_c = face_0_copy.m_c;



    new_face0.m_fa = face_touching_new_face_0;
    new_face1.m_fa = face_touching_new_face_1;
    new_face2.m_fa = face_touching_new_face_2;

    new_face0.m_fb = new_face1_index;
    new_face1.m_fb = new_face2_index;
    new_face2.m_fb = face_index;

    new_face0.m_fc = new_face2_index;
    new_face1.m_fc = face_index;
    new_face2.m_fc = new_face1_index;



    int vertex_opposing_new_face_0 = -1, vertex_opposing_new_face_1 = -1, vertex_opposing_new_face_2 = -1;

    if (face_touching_new_face_0 != -1)
        vertex_opposing_new_face_0 = m_faces[face_touching_new_face_0].find_local_vertex_index_with_opposing_face(face_index);
    if (face_touching_new_face_1 != -1)
        vertex_opposing_new_face_1 = m_faces[face_touching_new_face_1].find_local_vertex_index_with_opposing_face(face_index);
    if (face_touching_new_face_2 != -1)
        vertex_opposing_new_face_2 = m_faces[face_touching_new_face_2].find_local_vertex_index_with_opposing_face(face_index);

    if (vertex_opposing_new_face_0 != -1)
        m_faces[face_touching_new_face_0].opposing_face(vertex_opposing_new_face_0) = face_index;
    if (vertex_opposing_new_face_1 != -1)
        m_faces[face_touching_new_face_1].opposing_face(vertex_opposing_new_face_1) = new_face1_index;
    if (vertex_opposing_new_face_2 != -1)
        m_faces[face_touching_new_face_2].opposing_face(vertex_opposing_new_face_2) = new_face2_index;

    //If this vertex (the vertex shared by the new face1 and new face2 and not the point we just)
    //inserted was pointing to face 0, the face 0 is now too far away so this vertex is going to
    //have to point to another adjacent face (we'll give face1_index but it could also have been
    //face2 index)
    Vertex & vertex = m_vertices[face_0_copy.m_b];
    if (vertex.get_adjacent_face_index() == face_index)
        vertex.set_adjacent_face_index(new_face1_index);

    m_faces.push_back(new_face1);
    m_faces.push_back(new_face2);

    if (apply_lawson_after_insertion)
    {
        std::vector<std::pair<int, int>> edges_to_flip;
        edges_to_flip.push_back(std::make_pair(face_0_copy.m_a, face_0_copy.m_b));
        edges_to_flip.push_back(std::make_pair(face_0_copy.m_b, face_0_copy.m_c));
        edges_to_flip.push_back(std::make_pair(face_0_copy.m_c, face_0_copy.m_a));

        check_and_flip_multiple_edges(edges_to_flip);
    }
}

void Mesh::check_and_flip_multiple_edges(std::vector<std::pair<int, int>>& edges_to_flip)
{
    std::queue<std::pair<int, int>> to_flip_queue;
    //Used to make sure that we don't have the same edge to
    //flip multiple times in the queue (because this would mean that the
    //second time we encounter the same edge in the queue, we would be flipping
    //an edge that doesn't exist anymore since it was already flipped once
    //[and so the two vertices indices don't refer to a valid edge anymore])
    std::unordered_set<std::pair<int, int>, Mesh::pair_hash> to_flip_unique;

    //Inserting the edges to flip in the queue
    for (std::pair<int, int>& edge_to_flip : edges_to_flip)
    {
        std::pair<int, int> ordered_edge = std::make_pair(std::min(edge_to_flip.first, edge_to_flip.second),
                                                          std::max(edge_to_flip.first, edge_to_flip.second));

        //If the edge isn't already in the queue
        if (to_flip_unique.find(ordered_edge) == to_flip_unique.end())
        {
            to_flip_queue.push(ordered_edge);
            to_flip_unique.insert(ordered_edge);
        }
    }

    //Flipping the edges
    while (!to_flip_queue.empty())
    {
        std::pair<int, int> edge_to_flip = to_flip_queue.front();
        to_flip_queue.pop();
        to_flip_unique.erase(to_flip_unique.find(edge_to_flip));

        if (!is_edge_locally_delaunay(edge_to_flip))
        {
            std::vector<std::pair<int, int>> affected_edges = edge_flip(edge_to_flip);
            for (std::pair<int, int> affected_edge : affected_edges)
            {
                std::pair<int, int> ordered_edge = std::make_pair(std::min(affected_edge.first, affected_edge.second),
                                                                  std::max(affected_edge.first, affected_edge.second));
                if (to_flip_unique.find(ordered_edge) == to_flip_unique.end())
                {
                    to_flip_queue.push(ordered_edge);
                    to_flip_unique.insert(ordered_edge);
                }
            }
        }
    }
}

std::vector<std::pair<int, int>> Mesh::edge_flip(const std::pair<int, int>& vertex_index_pair)
{
    std::pair<int, int> faces_indices = get_faces_indices_of_edge(vertex_index_pair);

    return edge_flip(faces_indices.first, faces_indices.second);
}

std::vector<std::pair<int, int>> Mesh::edge_flip(const int face_0_index, const int face_1_index)
{
    Face face_0_copy = m_faces[face_0_index];
    Face face_1_copy = m_faces[face_1_index];

    Face& face_0 = m_faces[face_0_index];
    Face& face_1 = m_faces[face_1_index];

    int local_vertex_on_face_0_opposing_to_face_1 = face_0_copy.find_local_vertex_index_with_opposing_face(face_1_index);
    int local_vertex_on_face_1_opposing_to_face_0 = face_1_copy.find_local_vertex_index_with_opposing_face(face_0_index);

    //The edges that are going to be impacted by the flipping of the edge given in parameter
    std::vector<std::pair<int, int>> impacted_edges;
    impacted_edges.push_back(std::make_pair(face_0_copy.global_index_of_local_vertex_index((local_vertex_on_face_0_opposing_to_face_1 + 0) % 3),
                                            face_0_copy.global_index_of_local_vertex_index((local_vertex_on_face_0_opposing_to_face_1 + 1) % 3)));
    impacted_edges.push_back(std::make_pair(face_0_copy.global_index_of_local_vertex_index((local_vertex_on_face_0_opposing_to_face_1 + 2) % 3),
                                            face_0_copy.global_index_of_local_vertex_index((local_vertex_on_face_0_opposing_to_face_1 + 0) % 3)));

    impacted_edges.push_back(std::make_pair(face_1.global_index_of_local_vertex_index((local_vertex_on_face_1_opposing_to_face_0 + 0) % 3),
                                            face_1.global_index_of_local_vertex_index((local_vertex_on_face_1_opposing_to_face_0 + 1) % 3)));
    impacted_edges.push_back(std::make_pair(face_1.global_index_of_local_vertex_index((local_vertex_on_face_1_opposing_to_face_0 + 2) % 3),
                                            face_1.global_index_of_local_vertex_index((local_vertex_on_face_1_opposing_to_face_0 + 0) % 3)));

    if (local_vertex_on_face_0_opposing_to_face_1 == -1
     || local_vertex_on_face_1_opposing_to_face_0 == -1)
    {
        //The two faces are not adjacent, they do not form an edge

        std::cerr << "The two faces [" << face_0_index << ", " << face_1_index << "] given to edge_flip do not form an edge" << std::endl;
        return std::vector<std::pair<int, int>>();
    }

    int local_vertex1_on_face_0 = (local_vertex_on_face_0_opposing_to_face_1 + 1) % 3;
    int local_vertex2_on_face_0 = (local_vertex1_on_face_0 + 1) % 3;
    int local_vertex2_on_face_1 = (local_vertex_on_face_1_opposing_to_face_0 + 1) % 3;
    int local_vertex1_on_face_1 = (local_vertex2_on_face_1 + 1) % 3;

    face_0.global_index_of_local_vertex_index(local_vertex1_on_face_0) = face_1_copy.global_index_of_local_vertex_index(local_vertex_on_face_1_opposing_to_face_0);
    face_1.global_index_of_local_vertex_index(local_vertex2_on_face_1) = face_0_copy.global_index_of_local_vertex_index(local_vertex_on_face_0_opposing_to_face_1);

    int local_vertex0_on_face_0 = (local_vertex1_on_face_0 + 2) % 3;

    int face_2_index = face_0_copy.opposing_face(local_vertex1_on_face_0);
    int face_3_index = face_1_copy.opposing_face(local_vertex1_on_face_1);
    face_0.opposing_face(local_vertex0_on_face_0) = face_3_index;
    face_0.opposing_face((local_vertex0_on_face_0 + 1) % 3) = face_2_index;
    face_0.opposing_face((local_vertex0_on_face_0 + 2) % 3) = face_1_index;

    int face_4_index = face_1_copy.opposing_face(local_vertex2_on_face_1);
    int face_5_index = face_0_copy.opposing_face(local_vertex2_on_face_0);
    face_1.opposing_face(local_vertex1_on_face_1) = face_0_index;
    face_1.opposing_face((local_vertex1_on_face_1 + 1) % 3) = face_5_index;
    face_1.opposing_face((local_vertex1_on_face_1 + 2) % 3) = face_4_index;

    if (face_3_index != -1)
    {
        Face& face_3 = m_faces[face_3_index];
        int local_vertex4_on_face_3 = face_3.find_local_vertex_index_with_opposing_face(face_1_index);
        face_3.opposing_face(local_vertex4_on_face_3) = face_0_index;
    }

    if (face_5_index != -1)
    {
        Face& face_5 = m_faces[face_5_index];
        int local_vertex6_on_face_5 = face_5.find_local_vertex_index_with_opposing_face(face_0_index);
        face_5.opposing_face(local_vertex6_on_face_5) = face_1_index;
    }

    int global_index_of_vertex2_on_face_0 = face_0_copy.global_index_of_local_vertex_index(local_vertex2_on_face_0);
    int global_index_of_vertex0_on_face_1 = face_1_copy.global_index_of_local_vertex_index(local_vertex1_on_face_1);
    Vertex& vertex2_on_face_0 = m_vertices[global_index_of_vertex2_on_face_0];
    Vertex& vertex0_on_face_1 = m_vertices[global_index_of_vertex0_on_face_1];

    //If the first vertex on the edge between the faces that have been flipped was pointing to one of the two
    //faces, we're going to have to check whether the face that it was pointing to is still adjacent to the vertex
    bool was_pointing_to_face_0 = false;
    bool was_pointing_to_face_1 = false;
    bool needs_update = true;
    if (vertex2_on_face_0.get_adjacent_face_index() == face_0_index)
    {
        was_pointing_to_face_0 = true;
        //We're going to update the adjacent face pointed to by the vertices that composED the edge before the flip
        //before now that the faces are flipped, the face pointed to by the vertex may not be adjacent anymore

        //Looking at the vertices of the face 0
        for (int i = 0; i < 3; i++)
        {
            if(face_0.global_index_of_local_vertex_index(i) == global_index_of_vertex2_on_face_0)
            {
                //We found the vertex on the flipped face so there's nothing to update, the vertex is still
                //adjacent to the flipped face
                needs_update = false;

                break;
            }
        }

        //We didn't find the vertex on the flipped face. This means that the
        //vertex was pointing to the face that is now innaccessible so
        //we're going to have to update the vertex's adjacent face
        //needs_update = true; //Redundant but just for clarity
    }
    else if (vertex2_on_face_0.get_adjacent_face_index() == face_1_index)
    {
        was_pointing_to_face_1 = true;

        for (int i = 0; i < 3; i++)
        {
            if(face_1.global_index_of_local_vertex_index(i) == global_index_of_vertex2_on_face_0)
            {
                needs_update = false;

                break;
            }
        }
    }

    if (was_pointing_to_face_0 && needs_update)
        //Because the face 0 is no longer adjacent the vertex, this means that the vertex has to be pointing
        //to the other face
        vertex2_on_face_0.set_adjacent_face_index(face_1_index);
    else if (was_pointing_to_face_1 && needs_update)
        vertex2_on_face_0.set_adjacent_face_index(face_0_index);



    //Same check for the other vertex of the edge that was flipped
    was_pointing_to_face_0 = false;
    was_pointing_to_face_1 = false;
    needs_update = true;
    if (vertex0_on_face_1.get_adjacent_face_index() == face_0_index)
    {
        was_pointing_to_face_0 = true;
        //We're going to update the adjacent face pointed to by the vertices that composED the edge before the flip
        //before now that the faces are flipped, the face pointed to by the vertex may not be adjacent anymore

        //Looking at the vertices of the face 0
        for (int i = 0; i < 3; i++)
        {
            if(face_0.global_index_of_local_vertex_index(i) == global_index_of_vertex0_on_face_1)
                //We found the vertex on the flipped face so there's nothing to update, the vertex is still
                //adjacent to the flipped face
                return impacted_edges;
        }

        //We didn't find the vertex on the flipped face. This means that the
        //vertex was pointing to the face that is now innaccessible so
        //we're going to have to update the vertex's adjacent face
        //needs_update = true; //Redundant but just for clarity
    }
    else if (vertex0_on_face_1.get_adjacent_face_index() == face_1_index)
    {
        was_pointing_to_face_1 = true;

        for (int i = 0; i < 3; i++)
            if(face_1.global_index_of_local_vertex_index(i) == global_index_of_vertex0_on_face_1)
                return impacted_edges;
    }

    if (was_pointing_to_face_0 && needs_update)
        //Because the face 0 is no longer adjacent the vertex, this means that the vertex has to be pointing
        //to the other face
        vertex0_on_face_1.set_adjacent_face_index(face_1_index);
    else if (was_pointing_to_face_1 && needs_update)
        vertex0_on_face_1.set_adjacent_face_index(face_0_index);

    return impacted_edges;
}

void Mesh::insert_point_2D(const Point &point, bool apply_lawson_after_insertion)
{
    if (m_vertices.size() < 3)
    {
        m_vertices.push_back(Vertex(-1, point));

        //We just inserted the third point of our mesh, we can create the first face
        if (m_vertices.size() == 3)
        {
            m_vertices[0].set_adjacent_face_index(0);
            m_vertices[1].set_adjacent_face_index(0);
            m_vertices[2].set_adjacent_face_index(0);

            Point a, b, c;
            a = m_vertices[0].get_point();
            b = m_vertices[1].get_point();
            c = m_vertices[2].get_point();

            //Creating the face with the right orientation
            if (Point::orientation_test(a, b, c) > 0)
                m_faces.push_back(Face(0, 1, 2, -1, -1, -1));
            else if (Point::orientation_test(a, c, b) > 0)
                m_faces.push_back(Face(0, 2, 1, -1, -1, -1));
        }

        return;
    }

    //Determining whether the point that we want to insert is inside the triangulation or not
    int current_face_index = 0;
    Face current_face = m_faces[0];
    bool point_inside_triangulation = point_inside_triangulation_brute_force(point);
    if (point_inside_triangulation)
    {
        for (int i = 0; i < 3; i++)
        {
            Point current_face_barycenter = barycenter_of_face(current_face);
            Segment from_to_segment(current_face_barycenter, point);

            Point segment_point_a = m_vertices[current_face.global_index_of_local_vertex_index(i)].get_point();
            Point segment_point_b = m_vertices[current_face.global_index_of_local_vertex_index((i + 1) % 3)].get_point();

            Segment edge = Segment(segment_point_a, segment_point_b);
            if (edge.intersect(from_to_segment))
            {
                current_face_index = current_face.opposing_face((i + 2) % 3);
                if (current_face_index == -1)
                    break;

                current_face = m_faces[current_face_index];

                //We're going to start looping over the edges of the new face again so
                //we position i to -1 so that it gets incremented to 0 to start over
                i = -1;
            }
        }
    }

    if (point_inside_triangulation)
    {
        Point a, b, c;
        a = m_vertices[current_face.m_a].get_point();
        b = m_vertices[current_face.m_b].get_point();
        c = m_vertices[current_face.m_c].get_point();

        if (Point::is_point_in_triangle(point, a, b, c))
            face_split(current_face_index, point, apply_lawson_after_insertion);
    }
    else
        insert_outside_convex_hull_2D(point, apply_lawson_after_insertion);
}

bool Mesh::point_inside_triangulation_brute_force(const Point& point)
{
    int debug_index = 0;
    for (Face& face : m_faces)
    {
        debug_index++;
        Point A, B, C;
        A = m_vertices[face.m_a].get_point();
        B = m_vertices[face.m_b].get_point();
        C = m_vertices[face.m_c].get_point();

        if (Point::is_point_in_triangle(point, A, B, C) > 0)
            return true;
    }

    return false;
}

bool Mesh::is_edge_locally_delaunay(const std::pair<int, int> two_vertex_indices)
{
    std::pair<int, int> faces_indices = get_faces_indices_of_edge(two_vertex_indices);

    return is_edge_locally_delaunay(faces_indices.first, faces_indices.second);
}

bool Mesh::is_edge_locally_delaunay(int face1_index, int face2_index)
{
    //Edge on the convex hull
    if (face1_index == -1 || face2_index == -1)
        return true;

    const Face& face1 = m_faces[face1_index];
    const Face& face2 = m_faces[face2_index];

    //Finding the points opposite to the edge on the face 1 and face 2
    int vertex_index_on_face1_opposite_to_face2 = face1.find_local_vertex_index_with_opposing_face(face2_index);
    int vertex_index_on_face2_opposite_to_face1 = face2.find_local_vertex_index_with_opposing_face(face1_index);

    Point opposite_vertex_on_face1 = m_vertices[face1.global_index_of_local_vertex_index(vertex_index_on_face1_opposite_to_face2)].get_point();
    Point opposite_vertex_on_face2 = m_vertices[face2.global_index_of_local_vertex_index(vertex_index_on_face2_opposite_to_face1)].get_point();

    return !get_circumscribed_circle_of_face(face1).contains_point(opposite_vertex_on_face2)
        && !get_circumscribed_circle_of_face(face2).contains_point(opposite_vertex_on_face1);
}

void Mesh::delaunayize_lawson()
{
    std::queue<std::pair<int, int>> to_flip;
    //Set to check whether an edge is already in the queue to be flipped
    std::unordered_set<std::pair<int, int>, Mesh::pair_hash> to_flip_unique;

    for (Face& face : m_faces)
    {
        std::vector<std::pair<int, int>> edges_of_face = get_edges_of_face(face);
        for (std::pair<int, int>& edge : edges_of_face)
        {
            std::pair<int, int> edge_reordered = std::make_pair(std::min(edge.first, edge.second), std::max(edge.first, edge.second));

            //This edge isn't alredy in the queue to be flipped
            if (to_flip_unique.find(edge_reordered) == to_flip_unique.end())
            {
                to_flip.push(edge_reordered);
                to_flip_unique.insert(edge_reordered);
            }
        }
    }

    int debug_index = 0;
    while (!to_flip.empty())
    {
        debug_index++;
        std::pair<int, int> edge_to_flip = to_flip.front();
        to_flip_unique.erase(to_flip_unique.find(edge_to_flip));
        to_flip.pop();

        if (debug_index > 1565)
            save_as_off("flipped_debug_" + std::to_string(debug_index) + ".off");

        //If the edge still needs to be flipped (because it may have been put in
        //the queue before but since then, other edges have been flipped which has
        //made this edge locally delaunay and thus doesn't need to be flipped anymore)
        if (!is_edge_locally_delaunay(edge_to_flip))
        {
            std::vector<std::pair<int, int>> impacted_edges = edge_flip(edge_to_flip);

            for (std::pair<int, int>& new_edge_to_flip : impacted_edges)
            {
                std::pair<int, int> edge_reordered = std::make_pair(std::min(new_edge_to_flip.first, new_edge_to_flip.second), std::max(new_edge_to_flip.first, new_edge_to_flip.second));

                //This edge isn't alredy in the queue to be flipped
                if (to_flip_unique.find(edge_reordered) == to_flip_unique.end())
                {
                    to_flip.push(edge_reordered);
                    to_flip_unique.insert(edge_reordered);
                }
            }
        }
    }
}

void Mesh::ruppert(const std::vector<Segment>& constraint_segments)
{
    std::unordered_set<int> constraint_segments_found;
    int face_index = 0;
    for (const Face& face : m_faces)
    {
        std::pair<int, int> edge1 = std::make_pair(face.m_a, face.m_b);
        std::pair<int, int> edge2 = std::make_pair(face.m_b, face.m_c);
        std::pair<int, int> edge3 = std::make_pair(face.m_c, face.m_a);

        Segment face_segment1 = Segment(m_vertices[edge1.first].get_point(), m_vertices[edge1.second].get_point());
        Segment face_segment2 = Segment(m_vertices[edge2.first].get_point(), m_vertices[edge2.second].get_point());
        Segment face_segment3 = Segment(m_vertices[edge3.first].get_point(), m_vertices[edge3.second].get_point());
        //TODO mauvaise complexité, ça serait mieux d'avoir un set des segments
        //pour pouvoir instant lookup plutot que de faire une grosse boucle a chaque
        //fois

        //Finding the given constraint segments that are not already in the triangulation
        for (int segment_index = 0; segment_index < constraint_segments.size(); segment_index++)
        {
            const Segment& segment = constraint_segments[segment_index];
            if (segment == face_segment1 || segment == face_segment2 || segment == face_segment3)
                constraint_segments_found.insert(segment_index);
        }
        face_index++;
    }

    //Making a queue from the segments that are not in the triangulation
    std::queue<Segment> constraint_segments_queue;
    for (int segment_index = 0; segment_index < constraint_segments.size(); segment_index++)
        if (constraint_segments_found.find(segment_index) == constraint_segments_found.end())
            constraint_segments_queue.push(constraint_segments[segment_index]);

    while (!constraint_segments_queue.empty())
    {
        const Segment& segment = constraint_segments_queue.front();
        constraint_segments_queue.pop();

        Point middle_point = (segment.a + segment.b) / 2;
        insert_point_2D(middle_point, true);

        constraint_segments_queue.push(Segment(segment.a, middle_point));
        constraint_segments_queue.push(Segment(middle_point, segment.b));
    }
}

void Mesh::insert_outside_convex_hull_2D(const Point &point, bool apply_lawson_after_insertion)
{
    Vertex new_vertex = Vertex(m_faces.size(), point);
    int new_vertex_index = m_vertices.size();
    m_vertices.push_back(new_vertex);

    //[global vertex index of new edge 1, global vertex index of new edge 2] -> [global vertex index to update, face index of vertex to update]
    std::map<std::pair<int, int>, std::pair<int, int>> new_edges_to_vertex_to_update_and_face;
    //List of the opposing faces of vertices that will need to be updated after the creation of
    //the new faces
    //This is needed because when we're creating new faces, we will have vertices that
    //had -1 for opposing face that will now have a new face instead so we will
    //need to update those vertices. We're doing it afterwards instead of during
    //the algorithm because we don't want to mess with the convex_hull_edges iterator
    //[face_index, [local_vertex_index_in_the_face, new_opposing_face]]
    std::vector<std::pair<int, std::pair<int, int>>> opposing_faces_to_update;





    // Looking for an edge that is visible from the point we want to insert
    auto convex_hull_edge_ite = convex_hull_edges_in_order_begin();
    bool non_visible_found = false;
    bool visible = false;
    do
    {
        std::pair<int, int> convex_hull_edge = *convex_hull_edge_ite;

        Point edge_point1 = m_vertices[convex_hull_edge.first].get_point();
        Point edge_point2 = m_vertices[convex_hull_edge.second].get_point();

        visible = Point::orientation_test(point, edge_point2, edge_point1) > 0;
        if (!visible)
            non_visible_found = true;

        if (visible && non_visible_found)
            //If we found our edge, we don't to increment the iterator any further
            //so we're getting out of the loop
            break;
        else
            ++convex_hull_edge_ite;
    } while (!visible || !non_visible_found);

    auto first_visible_edge = convex_hull_edge_ite.get_current_edge();




    convex_hull_edge_ite = convex_hull_edges_in_order_begin(first_visible_edge);
    auto convex_hull_edge_ite_end = convex_hull_edges_in_order_past_the_end();
    int edge_index = 0;
    std::vector<Face> new_faces_to_add;
    std::vector<std::pair<int, int>> edges_to_check; //Edges that will need to be flipped
    //if they are not locally delaunay
    for (; convex_hull_edge_ite != convex_hull_edge_ite_end;
           convex_hull_edge_ite++, edge_index++)
    {
        Iterator_on_convex_hull_edges_in_order::Convex_hull_edge convex_hull_edge_struct = convex_hull_edge_ite.get_current_edge();
        std::pair<int, int> convex_hull_edge = *convex_hull_edge_ite;
        int convex_hull_edge_face_index = convex_hull_edge_struct.face_index;

        Point edge_point1 = m_vertices[convex_hull_edge.first].get_point();
        Point edge_point2 = m_vertices[convex_hull_edge.second].get_point();

        if (Point::orientation_test(point, edge_point2, edge_point1) > 0)
        {
            //If the edge is visible from the point we want to insert

            Face new_face = Face(new_vertex_index, convex_hull_edge.second, convex_hull_edge.first, convex_hull_edge_face_index, -1, -1);
            edges_to_check.push_back(std::make_pair(new_vertex_index, convex_hull_edge.second));
            edges_to_check.push_back(std::make_pair(new_vertex_index, convex_hull_edge.first));
            edges_to_check.push_back(std::make_pair(convex_hull_edge.second, convex_hull_edge.first));

            int new_face_index = m_faces.size() + new_faces_to_add.size();

            //We just created a new face where there was nothing before so we're going to
            //update the opposing face of the vertex that was pointing to -1 before
            //and now is pointing to the new face
            //We're adding it to a list for later update (because if we update it now, it's
            //going to mess up with the convex_hull_edges iterator)
            opposing_faces_to_update.push_back(std::make_pair(convex_hull_edge_face_index,
                                               std::make_pair(convex_hull_edge_struct.local_index_of_vertex_opposite_to_edge, new_face_index)));

            std::pair<int, int> new_edge_1, new_edge_2;
            new_edge_1 = std::make_pair(new_vertex_index < convex_hull_edge.first ? new_vertex_index : convex_hull_edge.first,
                                        new_vertex_index < convex_hull_edge.first ? convex_hull_edge.first : new_vertex_index);
            new_edge_2 = std::make_pair(new_vertex_index < convex_hull_edge.second ? new_vertex_index : convex_hull_edge.second,
                                        new_vertex_index < convex_hull_edge.second ? convex_hull_edge.second : new_vertex_index);

            auto edge_1_find = new_edges_to_vertex_to_update_and_face.find(new_edge_1);
            if (edge_1_find != new_edges_to_vertex_to_update_and_face.end())
            {
                //This edge has already been created before and we just created a face along
                //that edge so we will set to opposite face of the vertex to update
                int face_index = edge_1_find->second.second;
                //Because the new faces haven't been added yet, we're
                //getting them from the next_faces_to_add array
                Face& face_of_opposing_vertex = new_faces_to_add[face_index - m_faces.size()];

                //////////////
                //TODO on ne positionne pas les opposing face comme il faut, la face 4, vertex local 1
                //devrait pointer sur -1 et pas 4
                //////////////
                face_of_opposing_vertex.opposing_face(face_of_opposing_vertex.local_index_of_global_vertex_index(edge_1_find->second.first)) = new_face_index;
                new_face.opposing_face(new_face.local_index_of_global_vertex_index(convex_hull_edge.second)) = face_index;
            }
            else
                new_edges_to_vertex_to_update_and_face.insert(std::make_pair(new_edge_1, std::make_pair(convex_hull_edge.second, new_face_index)));

            //Same for the other new edge
            auto edge_2_find = new_edges_to_vertex_to_update_and_face.find(new_edge_2);
            if (edge_2_find != new_edges_to_vertex_to_update_and_face.end())
            {
                int face_index = edge_2_find->second.second;
                Face& opposing_vertex_face = new_faces_to_add[face_index - m_faces.size()];

                opposing_vertex_face.opposing_face(opposing_vertex_face.local_index_of_global_vertex_index(edge_2_find->second.first)) = new_face_index;
                new_face.opposing_face(new_face.local_index_of_global_vertex_index(convex_hull_edge.first)) = face_index;
            }
            else
                new_edges_to_vertex_to_update_and_face.insert(std::make_pair(new_edge_2, std::make_pair(convex_hull_edge.first, new_face_index)));

            new_faces_to_add.push_back(new_face);
        }
        else
            break;
    }

    for (Face& face : new_faces_to_add)
        m_faces.push_back(face);

    for (std::pair<int, std::pair<int, int>>& to_update : opposing_faces_to_update)
        m_faces[to_update.first].opposing_face(to_update.second.first) = to_update.second.second;

    if (apply_lawson_after_insertion)
        check_and_flip_multiple_edges(edges_to_check);
}

void Mesh::scale_to_min_max_points(const Point &target_min, const Point &target_max)
{
    Point mesh_min(std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max()),
          mesh_max(std::numeric_limits<double>::min(), std::numeric_limits<double>::min(), std::numeric_limits<double>::min());

    for (const Vertex& vertex : m_vertices)
    {
        mesh_min = min(mesh_min, vertex.get_point());
        mesh_max = max(mesh_max, vertex.get_point());
    }

    double target_diagonal = length(target_min - target_max);
    double mesh_diagonal = length(mesh_min - mesh_max);
    double scale_needed = target_diagonal / mesh_diagonal;
    Point scaled_mesh_min = mesh_min * scale_needed;
    Vector translation_needed = target_min - scaled_mesh_min;

    for (Vertex& vertex : m_vertices)
        vertex.get_point() = vertex.get_point() * scale_needed + translation_needed;

}

Face& Mesh::Circulator_on_faces::operator*()
{
    return m_mesh->m_faces.at(m_current_face_index);
}

//Prefix
Mesh::Circulator_on_faces& operator++(Mesh::Circulator_on_faces& operand)
{
    Face& current_face = *operand;

    Point p_a, p_b, p_c;
    p_a = operand.m_mesh->m_vertices[current_face.m_a].get_point();
    p_b = operand.m_mesh->m_vertices[current_face.m_b].get_point();
    p_c = operand.m_mesh->m_vertices[current_face.m_c].get_point();

    //The index of the vertex we're circulating around in the current face we're on
    int index_of_main_vertex_in_current_face;

    if (p_a == operand.m_vertex_circulating_around->get_point())
        index_of_main_vertex_in_current_face = 0;
    else if (p_b == operand.m_vertex_circulating_around->get_point())
        index_of_main_vertex_in_current_face = 1;
    else
        index_of_main_vertex_in_current_face = 2;


    int next_vertex_index;
    if (!operand.m_circulating_backwards)
        next_vertex_index = (index_of_main_vertex_in_current_face + 1) % 3;
    else
        next_vertex_index = (index_of_main_vertex_in_current_face + 2) % 3;

    operand.m_current_face_index = current_face.opposing_face(next_vertex_index);
    if(operand.m_current_face_index == -1 && operand.m_circulating_backwards)
        //We circulated around everything
        operand.m_current_face_index = -1;
    else if (operand.m_current_face_index == -1)
    {
        //We arrived outside of the mesh, we're going to circulate the other way around from the starting face
        operand.m_current_face_index = operand.m_start_face_index;
        operand.m_circulating_backwards = true;

        ++operand;
    }

    return operand;
}

//Postfix
Mesh::Circulator_on_faces operator++(Mesh::Circulator_on_faces& operand, int dummy)
{
    Mesh::Circulator_on_faces copy = operand;

    ++operand;

    return copy;
}

bool operator ==(const Mesh::Circulator_on_faces& a, const Mesh::Circulator_on_faces& b)
{
    if (a.m_current_face_index == -1 && b.m_current_face_index == -1)
        return true;
    else if (a.m_current_face_index != b.m_current_face_index)
        return false;
    else if (a.m_mesh != b.m_mesh)
        return false;
    else if (a.m_vertex_circulating_around != b.m_vertex_circulating_around)
        return false;
    else if (a.m_current_face_index == -1 && b.m_current_face_index == -1 && a.m_mesh == b.m_mesh)
        return true;

    return true;
}

bool operator !=(const Mesh::Circulator_on_faces& a, const Mesh::Circulator_on_faces& b)
{
    return !(a == b);
}


Mesh::Iterator_on_convex_hull_edges_in_order::Iterator_on_convex_hull_edges_in_order(Mesh& mesh) : m_mesh(&mesh)
{
    int edge_face_index = m_mesh->get_one_convex_hull_face();
    m_current_edge = Convex_hull_edge
    {
        edge_face_index,
        m_mesh->m_faces[edge_face_index].find_local_vertex_index_with_opposing_face(-1)
    };

    m_start_edge = m_current_edge;
}

const std::pair<int, int> Mesh::Iterator_on_convex_hull_edges_in_order::operator*()
{
    Face& face = m_mesh->m_faces[m_current_edge.face_index];

    int first, second;
    first = face.global_index_of_local_vertex_index((m_current_edge.local_index_of_vertex_opposite_to_edge + 1) % 3);
    second = face.global_index_of_local_vertex_index((m_current_edge.local_index_of_vertex_opposite_to_edge + 2) % 3);

    return std::make_pair(first, second);
}

//Prefix
Mesh::Iterator_on_convex_hull_edges_in_order& operator++(Mesh::Iterator_on_convex_hull_edges_in_order& operand)
{
    Mesh::Iterator_on_convex_hull_edges_in_order::Convex_hull_edge current_edge = operand.get_current_edge();

    Face& start_face = operand.m_mesh->m_faces[current_edge.face_index];
    int global_index_of_vertex_rotating_around = start_face.global_index_of_local_vertex_index((current_edge.local_index_of_vertex_opposite_to_edge + 2) % 3);
    int next_vertex_local_index = (start_face.local_index_of_global_vertex_index(global_index_of_vertex_rotating_around) + 2) % 3;
    int next_face_index = start_face.opposing_face(next_vertex_local_index);
    int current_face_index = next_face_index;
    if (next_face_index == -1)
    {
        //The next convex hull edge is on the same face
        Mesh::Iterator_on_convex_hull_edges_in_order::Convex_hull_edge new_current_edge
        {
            current_edge.face_index,
            (current_edge.local_index_of_vertex_opposite_to_edge + 1) % 3
        };

        operand.m_current_edge = new_current_edge;
        if (operand.m_current_edge == operand.m_start_edge)
            operand.m_past_the_end = true;

        return operand;
    }

    while (next_face_index != -1)
    {
        Face& current_face = operand.m_mesh->m_faces[current_face_index];
        current_face_index = next_face_index;

        next_vertex_local_index = (current_face.local_index_of_global_vertex_index(global_index_of_vertex_rotating_around) + 2) % 3;
        next_face_index = current_face.opposing_face(next_vertex_local_index);
    }

    // We found the next face that is on the convex hull, getting the edge
    Mesh::Iterator_on_convex_hull_edges_in_order::Convex_hull_edge new_current_edge
    {
        current_face_index,
        next_vertex_local_index
    };

    operand.m_current_edge = new_current_edge;
    if (operand.m_current_edge == operand.m_start_edge)
        operand.m_past_the_end = true;

    return operand;
}

// Postfix
Mesh::Iterator_on_convex_hull_edges_in_order operator++(Mesh::Iterator_on_convex_hull_edges_in_order& operand, int dummy)
{
    Mesh::Iterator_on_convex_hull_edges_in_order copy = operand;

    ++operand;

    return copy;
}

const Mesh::Iterator_on_convex_hull_edges_in_order::Convex_hull_edge &Mesh::Iterator_on_convex_hull_edges_in_order::get_current_edge()
{
    return m_current_edge;
}

bool operator ==(const Mesh::Iterator_on_convex_hull_edges_in_order& a, const Mesh::Iterator_on_convex_hull_edges_in_order& b)
{
    if (a.m_past_the_end == b.m_past_the_end)
        return true;
    else
        return (a.m_current_edge == b.m_current_edge);

}

bool operator !=(const Mesh::Iterator_on_convex_hull_edges_in_order& a, const Mesh::Iterator_on_convex_hull_edges_in_order& b)
{
    return !(a == b);
}

bool operator==(const Mesh::Iterator_on_convex_hull_edges_in_order::Convex_hull_edge &a, const Mesh::Iterator_on_convex_hull_edges_in_order::Convex_hull_edge &b)
{
    return a.face_index == b.face_index && a.local_index_of_vertex_opposite_to_edge == b.local_index_of_vertex_opposite_to_edge;
}



void Mesh::add_vertex(const Vertex& vertex)
{
    m_vertices.push_back(vertex);
}

void Mesh::add_face(const Face& face)
{
    m_faces.push_back(face);
}

int Mesh::get_one_convex_hull_face()
{
    int face_index = 0;
    for (Face& face : m_faces)
    {
        if (face.find_local_vertex_index_with_opposing_face(-1) != -1)
            return face_index;

        face_index++;
    }

}

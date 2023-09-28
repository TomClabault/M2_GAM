#include "mesh.h"

//Laplacien
/**
    TODO laplacien
    TODO faire attention au point quand on flip une edge que un point qui avait pour face incidente le point 0, point maintenant sur la face 1 ou l'inverse (pliutôt ql'inverse potentiellement)
 */

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

Mesh::Iterator_on_vertices Mesh::vertices_begin()
{
    return Mesh::Iterator_on_vertices(*this);
}

Mesh::Iterator_on_vertices Mesh::vertices_past_the_end()
{
    return Mesh::Iterator_on_vertices(*this, true);
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
    Mesh::Circulator_on_faces circulator = Mesh::Circulator_on_faces(*this, vertex);
    Mesh::Circulator_on_faces circulator_begin = circulator;

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
    } while (circulator != circulator_begin);

    return laplacien_sum / (2 * neighbor_faces_area_sum);
}

void Mesh::face_split(const int face_index, const Point& new_point)
{

    Face current_face_copy = m_faces[face_index];
    Face& new_face0 = m_faces[face_index];
    Face new_face1, new_face2;

    int new_face1_index = m_faces.size();
    int new_face2_index = m_faces.size() + 1;

    int face_touching_new_face_0 = current_face_copy.opposing_face(1);
    int face_touching_new_face_1 = current_face_copy.opposing_face(2);
    int face_touching_new_face_2 = current_face_copy.opposing_face(0);

    m_vertices.push_back(Vertex(face_index, new_point));

    //The current face we're inserting into is going to become 'face 0' of the
    //3 new faces created
    //The new point is always the vertex 0 of the new faces
    new_face0.m_a = m_vertices.size() - 1;
    new_face1.m_a = m_vertices.size() - 1;
    new_face2.m_a = m_vertices.size() - 1;

    new_face0.m_b = current_face_copy.m_c;
    new_face1.m_b = current_face_copy.m_a;
    new_face2.m_b = current_face_copy.m_b;

    new_face0.m_c = current_face_copy.m_a;
    new_face1.m_c = current_face_copy.m_b;
    new_face2.m_c = current_face_copy.m_c;



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

    m_faces.push_back(new_face1);
    m_faces.push_back(new_face2);
}

void Mesh::edge_flip(const int face_0_index, const int face_1_index)
{
    Face face_0_copy = m_faces[face_0_index];
    Face face_1_copy = m_faces[face_1_index];

    Face& face_0 = m_faces[face_0_index];
    Face& face_1 = m_faces[face_1_index];

    int local_vertex_on_face_0_opposing_to_face_1 = face_0_copy.find_local_vertex_index_with_opposing_face(face_1_index);
    int local_vertex_on_face_1_opposing_to_face_0 = face_1_copy.find_local_vertex_index_with_opposing_face(face_0_index);

    int local_vertex1_on_face_0 = (local_vertex_on_face_0_opposing_to_face_1 + 1) % 3;
    int local_vertex2_on_face_0 = (local_vertex1_on_face_0 + 1) % 3;
    int local_vertex2_on_face_1 = (local_vertex_on_face_1_opposing_to_face_0 + 1) % 3;
    int local_vertex1_on_face_1 = (local_vertex2_on_face_1 + 1) % 3;

    face_0.global_index_of_local_vertex_index(local_vertex1_on_face_0) = face_1_copy.global_index_of_local_vertex_index(local_vertex_on_face_1_opposing_to_face_0);
    face_1.global_index_of_local_vertex_index(local_vertex2_on_face_1) = face_0_copy.global_index_of_local_vertex_index(local_vertex_on_face_0_opposing_to_face_1);

    int local_vertex0_on_face_0 = local_vertex1_on_face_0 - 1;
    if (local_vertex0_on_face_0 == -1)
        local_vertex0_on_face_0 = 2;

    int face_2_index = face_0_copy.opposing_face(local_vertex1_on_face_0);
    int face_3_index = face_1_copy.opposing_face(local_vertex1_on_face_1);
    face_0.opposing_face(local_vertex0_on_face_0) = face_3_index;
    face_0.opposing_face((local_vertex0_on_face_0 + 1) % 3) = face_2_index;
    face_0.opposing_face((local_vertex0_on_face_0 + 2) % 3) = face_1_index;

    int face_4_index = face_0_copy.opposing_face(local_vertex2_on_face_1);
    int face_5_index = face_0_copy.opposing_face(local_vertex2_on_face_0);
    face_1.opposing_face(local_vertex1_on_face_1) = face_0_index;
    face_1.opposing_face((local_vertex1_on_face_1 + 1) % 3) = face_5_index;
    face_1.opposing_face((local_vertex1_on_face_1 + 2) % 3) = face_4_index;

    Face& face_3 = m_faces[face_3_index];
    int local_vertex4_on_face_3 = face_3.find_local_vertex_index_with_opposing_face(face_1_index);
    face_3.opposing_face(local_vertex4_on_face_3) = face_0_index;

    Face& face_5 = m_faces[face_5_index];
    int local_vertex6_on_face_5 = face_5.find_local_vertex_index_with_opposing_face(face_0_index);
    face_5.opposing_face(local_vertex6_on_face_5) = face_1_index;
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

    int next_vertex_index = (index_of_main_vertex_in_current_face + 1) % 3;

    operand.m_current_face_index = current_face.opposing_face(next_vertex_index);

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
    if (a.m_current_face_index != b.m_current_face_index)
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






void Mesh::add_vertex(const Vertex& vertex)
{
    m_vertices.push_back(vertex);
}

void Mesh::add_face(const Face& face)
{
    m_faces.push_back(face);
}

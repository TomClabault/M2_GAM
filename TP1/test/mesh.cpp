#include "mesh.h"

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

    operand.m_current_face_index = current_face.m_opposing_faces(next_vertex_index);

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

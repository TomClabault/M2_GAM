#include "mesh.h"
#include "segment.h"

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

std::pair<int, int> Mesh::Iterator_on_edges::operator*()
{
    return m_current_edge;
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

        operand.m_current_edge = std::make_pair(edge_p1_smallest ? edge_point_1_global_index : edge_point_2_global_index,
                                                edge_p1_smallest ? edge_point_2_global_index : edge_point_1_global_index);

        if (operand.m_current_edge_in_current_face == 2)
        {
            operand.m_current_edge_in_current_face = 0;
            operand.m_current_face_index++;
            if ((unsigned long long int)operand.m_current_face_index == mesh->m_faces.size())
                operand.m_past_the_end = true;
        }
        else
            operand.m_current_edge_in_current_face++;
    }
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

Point Mesh::barycenter_of_face(const Face &face) const
{
    return (m_vertices[face.m_a].get_point() + m_vertices[face.m_b].get_point() + m_vertices[face.m_c].get_point()) / 3;
}

Circle Mesh::get_circumscribed_circle_of_face(const Face &face) const
{
    Point circle_center;

    Point point_a, point_b, point_c;
    point_a = m_vertices[face.m_a].get_point();
    point_b = m_vertices[face.m_b].get_point();
    point_c = m_vertices[face.m_c].get_point();

    Vector AB = point_b - point_a;
    Vector BC = point_c - point_b;
    Vector CA = point_a - point_c;

    float tan_A, tan_B, tan_C;
    tan_A = length(cross(AB, -CA)) / dot(AB, -CA);
    tan_B = length(cross(BC, -AB)) / dot(BC, -AB);
    tan_C = length(cross(CA, -BC)) / dot(CA, -BC);

    circle_center = Point(tan_B + tan_C, tan_C + tan_A, tan_A + tan_B);

    return Circle(circle_center, length(circle_center - point_a));
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

    int global_index_of_vertex2_on_face_0 = face_0_copy.global_index_of_local_vertex_index(local_vertex2_on_face_0);
    int global_index_of_vertex0_on_face_1 = face_1_copy.global_index_of_local_vertex_index(local_vertex1_on_face_1);
    Vertex& vertex2_on_face_0 = m_vertices[global_index_of_vertex2_on_face_0];
    Vertex& vertex0_on_face_1 = m_vertices[global_index_of_vertex0_on_face_1];

    //If the first vertex on the edge between the face that have been flipped was pointing to one of the two
    //faces, we're going to have to check whether the face that it was pointing to is still adjacent to the vertex
    bool was_pointing_to_face_0 = false;
    if (vertex2_on_face_0.get_adjacent_face_index() == face_0_index)
    {
        was_pointing_to_face_0 = true;
        //We're going to update the adjacent face pointed to by the vertices that composED the edge before the flip
        //before now that the faces are flipped, the face pointed to by the vertex may not be adjacent anymore

        //Looking at the vertices of the face 0
        for (int i = 0; i < 3; i++)
        {
            if(face_0.global_index_of_local_vertex_index(i) == global_index_of_vertex2_on_face_0)
                //We found the vertex on the flipped face so there's nothing to update, the vertex is still
                //adjacent to the flipped face
                return;
        }

        //We didn't find the vertex on the flipped face. This means
    }
    else if (vertex2_on_face_0.get_adjacent_face_index() == face_1_index)
    {
        for (int i = 0; i < 3; i++)
            if(face_1.global_index_of_local_vertex_index(i) == global_index_of_vertex2_on_face_0)
                return;
    }

    if (was_pointing_to_face_0)
        //Because the face 0 is no longer adjacent the vertex, this means that the vertex has to be pointing
        //to the other face
        vertex2_on_face_0.set_adjacent_face_index(face_1_index);
    else //was_pointing_to_face_1
        vertex2_on_face_0.set_adjacent_face_index(face_0_index);



    //Same check for the other vertex of the edge that was flipped
    was_pointing_to_face_0 = false;
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
                return;
        }

        //We didn't find the vertex on the flipped face. This means
    }
    else if (vertex0_on_face_1.get_adjacent_face_index() == face_1_index)
    {
        for (int i = 0; i < 3; i++)
            if(face_1.global_index_of_local_vertex_index(i) == global_index_of_vertex0_on_face_1)
                return;
    }

    if (was_pointing_to_face_0)
        //Because the face 0 is no longer adjacent the vertex, this means that the vertex has to be pointing
        //to the other face
        vertex0_on_face_1.set_adjacent_face_index(face_1_index);
    else //was_pointing_to_face_1
        vertex0_on_face_1.set_adjacent_face_index(face_0_index);
}

void Mesh::insert_point_2D(const Point &point)
{
    bool is_inside_triangulation = false;

    //Determining whether the point that we want to insert is inside the triangulation or not
    Face& arbitrary_start_face = m_faces[0];
    Point start_face_barycenter = barycenter_of_face(arbitrary_start_face);
    Segment from_to_segment(start_face_barycenter, point);

    int current_face_index = 0;
    Face& current_face = arbitrary_start_face;
    for (int i = 0; i < 3; i++)
    {
        Point segment_point_a = m_vertices[current_face.global_index_of_local_vertex_index(i)].get_point();
        Point segment_point_b = m_vertices[current_face.global_index_of_local_vertex_index((i + 1) % 3)].get_point();

        Segment edge = Segment(segment_point_a, segment_point_b);
        if (edge.intersect(from_to_segment))
        {
            current_face_index = current_face.opposing_face((i + 2) % 3);
            current_face = m_faces[current_face_index];

            i = 0;
        }
    }

    Point a, b, c;
    a = m_vertices[current_face.m_a].get_point();
    b = m_vertices[current_face.m_b].get_point();
    c = m_vertices[current_face.m_c].get_point();

    if (Point::is_point_in_triangle(point, a, b, c))
        face_split(current_face_index, point);
    else
    {
        //We're going to have to add the point outside of the convex hull of the current
        //mesh

        //TODO parcourir les aretes au bord du mesh et voir si elles sont visibles par le point a inserer

        //DO WE NEED INFINITE FACES TO FIND THE EDGES THAT ARE AT THE BOUNDARY OF THE MESH ?
    }
}

bool Mesh::is_edge_locally_delaunay(int face1_index, int face2_index) const
{
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
    std::vector<std::pair<int, int>> to_flip;
    for (auto edge = edges_begin(); edge != edges_past_the_end(); edge++)
    {
        if (is_edge_locally_delaunay(*edge)) --------------
    }
}

void Mesh::compute_convex_hull_edges()
{
    for (int face_index = 0; face_index < m_faces.size(); face_index++)
    {
        Face& face = m_faces[face_index];

        if (face.m_fa == -1)
            m_convex_hull_edges.push_back(std::make_pair(face.m_b, face.m_c));
        if (face.m_fb == -1)
            m_convex_hull_edges.push_back(std::make_pair(face.m_c, face.m_a));
        if (face.m_fc == -1)
            m_convex_hull_edges.push_back(std::make_pair(face.m_a, face.m_b));
        m_convex_hull_edges_faces.push_back(face_index);
    }
}

void Mesh::insert_outside_convex_hull_2D(const Point &point)
{
    Vertex new_vertex = Vertex(m_faces.size(), point);
    int new_vertex_index = m_vertices.size();
    m_vertices.push_back(new_vertex);

    //[global vertex index of new edge 1, global vertex index of new edge 1] -> [global vertex index to update, face index of vertex to update]
    std::map<std::pair<int, int>, std::pair<int, int>> new_edges_to_vertex_to_update_and_face;
    //[global vertex index of new edge 1, global vertex index of new edge 1] -> is the edge part of the convex hull
    std::map<std::pair<int, int>, bool> part_of_new_convex_hull;
    //Face associated to the new edges that are now part of  the convex hull
    std::vector<int> part_of_new_convex_hull_associated_face;

    std::vector<std::list<std::pair<int, int>>::iterator> convex_hull_edges_to_remove;
    std::vector<std::list<int>::iterator> convex_hull_edges_faces_to_remove;

    int edge_index = 0;
    auto convex_hull_edge_ite = m_convex_hull_edges.begin();
    auto convex_hull_edge_face_ite = m_convex_hull_edges_faces.begin();
    for (; convex_hull_edge_ite != m_convex_hull_edges.end();
           convex_hull_edge_ite++, edge_index++)
    {
        std::pair<int, int> convex_hull_edge = *convex_hull_edge_ite;
        int convex_hull_edge_face_index = *convex_hull_edge_face_ite;

        Point edge_point1 = m_vertices[convex_hull_edge.first].get_point();
        Point edge_point2 = m_vertices[convex_hull_edge.second].get_point();

        //If the edge is visible from the point we want to insert
        if (Point::orientation_test(point, edge_point2, edge_point1) > 0)
        {
            Face new_face = Face(new_vertex_index, convex_hull_edge.first, convex_hull_edge.second, convex_hull_edge_face_index, -1, -1);
            int new_face_index = m_faces.size();
            m_faces.push_back(new_face);

            convex_hull_edges_to_remove.push_back(convex_hull_edge_ite);
            convex_hull_edges_faces_to_remove.push_back(convex_hull_edge_face_ite);

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
                Face& opposing_vertex_face = m_faces[edge_1_find->second.second];
                opposing_vertex_face.opposing_face(opposing_vertex_face.local_index_of_global_vertex_index(edge_1_find->second.first)) = new_face_index;

                part_of_new_convex_hull.find(new_edge_1)->second = false;
                part_of_new_convex_hull_associated_face.push_back(new_face_index);
            }
            else
            {
                new_edges_to_vertex_to_update_and_face.insert(std::make_pair(new_edge_1, std::make_pair(convex_hull_edge.second, new_face_index)));
                part_of_new_convex_hull.insert(std::make_pair(new_edge_1, true));
            }

            //Same for the other new edge
            auto edge_2_find = new_edges_to_vertex_to_update_and_face.find(new_edge_2);
            if (edge_2_find != new_edges_to_vertex_to_update_and_face.end())
            {
                Face& opposing_vertex_face = m_faces[edge_2_find->second.second];
                opposing_vertex_face.opposing_face(opposing_vertex_face.local_index_of_global_vertex_index(edge_2_find->second.first)) = new_face_index;

                part_of_new_convex_hull.find(new_edge_2)->second = false;
                part_of_new_convex_hull_associated_face.push_back(new_face_index);
            }
            else
            {
                new_edges_to_vertex_to_update_and_face.insert(std::make_pair(new_edge_2, std::make_pair(convex_hull_edge.second, new_face_index)));
                part_of_new_convex_hull.insert(std::make_pair(new_edge_2, true));
            }
        }
    }

    //Removing the edges that are not part of the convex hull anymore
    for (auto& convex_hull_edge_to_remove : convex_hull_edges_to_remove)
        m_convex_hull_edges.remove(*convex_hull_edge_to_remove);

    //And the associated faces
    for (auto& convex_hull_edge_face_to_remove : convex_hull_edges_faces_to_remove)
        m_convex_hull_edges_faces.remove(*convex_hull_edge_face_to_remove);

    //Adding the new edges that are now part of the convex hull
    int index = 0;
    for (auto new_edge : part_of_new_convex_hull)
    {
        if (new_edge.second == true)
        {
            m_convex_hull_edges.push_back(new_edge.first);

            //And the associated face
            m_convex_hull_edges_faces.push_back(part_of_new_convex_hull_associated_face[index]);
        }

        index++;
    }
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

void Mesh::push_convex_hull_edge(int index_vertex1, int index_vertex2)
{
    m_convex_hull_edges.push_back(std::make_pair(index_vertex1, index_vertex2));
}

void Mesh::push_convex_hull_edge_face(int face_index)
{
    m_convex_hull_edges_faces.push_back(face_index);
}

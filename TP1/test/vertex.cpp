#include "vertex.h"

#include <ostream>

std::ostream& operator <<(std::ostream& os, const Vertex& vertex)
{
    os << "Vertex[adj_face_index=" << vertex.m_adjacent_face_index << ", point=" << vertex.m_geometric_point << "]";

    return os;
}

Vector operator-(const Vertex &a, const Vertex &b)
{
    return a.get_point() - b.get_point();
}

#ifndef VERTEX_H
#define VERTEX_H

#include "point.h"

#include <ostream>

class Vertex
{
public:
    Vertex(int adjacent_face_index, Point geometric_coordinates) : m_adjacent_face_index(adjacent_face_index), m_geometric_point(geometric_coordinates) {}

    int get_adjacent_face_index() { return m_adjacent_face_index; }
    void set_adjacent_face_index(int adjacent_face_index) { m_adjacent_face_index = adjacent_face_index; }

    Point& get_point() { return m_geometric_point; }

    friend std::ostream& operator <<(std::ostream& os, const Vertex& vertex);

private:
    //Index of a face that contains the current vertex
    int m_adjacent_face_index;

    Point m_geometric_point;
};

#endif // VERTEX_H

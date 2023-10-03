#ifndef FACE_H
#define FACE_H

#include <ostream>

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
    Face() : m_a(-1), m_b(-1), m_c(-1), m_fa(-1), m_fb(-1), m_fc(-1) {}

    int local_index_of_global_vertex_index(int global_vertex_index)
    {
        return (m_a == global_vertex_index) ? 0 : (m_b == global_vertex_index) ? 1 : (m_c == global_vertex_index) ? 2 : -1;
    }

    /**
     * @return The local index of the vertex that has for opposing face 'face_index'
     * -1 if no vertex has for opposing 'face_index'
     */
    int find_local_vertex_index_with_opposing_face(int face_index) const
    {
        return (face_index == m_fa) ? 0 : (m_fb == face_index) ? 1 : (face_index == m_fc ? 2 : -1);
    }

    int opposing_face(int index) const
    {
        return *(&m_fa + index);
    }

    int& opposing_face(int index)
    {
        return *(&m_fa + index);
    }

    int global_index_of_local_vertex_index(int index) const
    {
        return *(&m_a + index);
    }


    int& global_index_of_local_vertex_index(int index)
    {
        return *(&m_a + index);
    }

    //Vertices of the face
    int m_a, m_b, m_c;

    //Adjacent faces. The faces are opposing the corresponding vertex:
    //fa is the face opposing the vertex 'a' (so it's the face that shares the 'bc' line)
    //...
    int m_fa, m_fb, m_fc;

    friend std::ostream& operator << (std::ostream& os, const Face& face);
};

#endif // FACE_H

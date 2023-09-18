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

    int& m_opposing_faces(int index) {return *(&m_fa + index);}
    int& m_vertices(int index) {return *(&m_a + index);}

    //Vertices of the face
    int m_a, m_b, m_c;

    //Adjacent faces. The faces are opposing the corresponding vertex:
    //fa is the face opposing the vertex 'a' (so it's the face that shares the 'bc' line)
    //...
    int m_fa, m_fb, m_fc;

    friend std::ostream& operator << (std::ostream& os, const Face& face);
};

#endif // FACE_H

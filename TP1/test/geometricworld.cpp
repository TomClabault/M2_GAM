#include "geometricworld.h"
#include "segment.h"
#include "TP1/offreader.h"

#include <iostream>

// Draw a Point
void glPointDraw(const Point & p) {
    glVertex3f(p.x, p.y, p.z);
}

void GeometricWorld::load_off(const char* filepath)
{
    //Exemple simple de boucle infinie lors de l'application de l'algorithme de Lawnson
    //à cause d'une face dégénérée
//    std::srand(5975);
//    for (int i = 0; i < 5; i++)
//    {
//        double x = std::rand() / (double)RAND_MAX;
//        double y = std::rand() / (double)RAND_MAX;

//        _mesh.insert_point_2D(Point(x, y, 0), true);
//        _mesh.save_as_off("a_random_mesh" + std::to_string(i) + ".off");
//    }

    //_mesh = OffReader::read_off(filepath);
    //std::cout << std::filesystem::current_path() << std::endl;
    //_mesh = OffReader::read_off("cube_maillage_triangles.off");
    //_mesh = OffReader::read_off("queen.off");
    //_mesh.insert_point_cloud("../test/data/alpes_random_2.txt");
    //_mesh.scale_to_min_max_points(Point(-3, -3, -3), Point(3, 3, 3));

    std::vector<Point> points;
    points.push_back(Point(0, 0, 0));
    points.push_back(Point(1, 0, 0));
    points.push_back(Point(2, 0.5, 0));
    points.push_back(Point(1.5, 1.5, 0));
    points.push_back(Point(0.5, 2, 0));
    points.push_back(Point(-0.5, 2, 0));
    points.push_back(Point(-0.5, 1, 0));

    for (const Point& point : points)
        _mesh.insert_point_2D(point, true);

    std::vector<Segment> constraint_segments;
    for (int i = 0; i < points.size() - 1; i++)
        constraint_segments.push_back(Segment(points[i + 0], points[i + 1]));
    constraint_segments.push_back(Segment(Point(0.5, 0.5, 0), Point(0.5, 1.5, 0)));

    _mesh.ruppert(constraint_segments);

    //_mesh.add_vertex(Vertex(0, Point(0, 0, 0)));
    //_mesh.add_vertex(Vertex(0, Point(1, 0, 0)));
    //_mesh.add_vertex(Vertex(0, Point(0.5, 0.5, 0)));
    //_mesh.add_vertex(Vertex(1, Point(1.5, 0.25, 0)));
    //_mesh.add_vertex(Vertex(2, Point(-0.5, 0.25, 0)));
    //_mesh.add_vertex(Vertex(2, Point(0.5, -0.5, 0)));

    //_mesh.add_face(Face(0, 1, 2, 1, 2, -1));
    //_mesh.add_face(Face(1, 3, 2, -1, 0, -1));
    //_mesh.add_face(Face(0, 2, 4, -1, -1, 0));
//    _mesh.add_face(Face(0, 1, 2, 1, 2, 3));
//    _mesh.add_face(Face(1, 3, 2, -1, 0, 4));
//    _mesh.add_face(Face(0, 2, 4, -1, -1, 0));
//    _mesh.add_face(Face(5, 1, 0, 0, -1, 4));
//    _mesh.add_face(Face(5, 3, 1, 1, 3, -1));


    //_mesh.insert_outside_convex_hull_2D(Point(0.5, -0.5, 0)); //5
    //_mesh.insert_outside_convex_hull_2D(Point(-0.5, -0.6, 0));//6
    //_mesh.insert_outside_convex_hull_2D(Point(1.0, -0.6, 0)); //7
    //_mesh.insert_outside_convex_hull_2D(Point(2.0, 0, 0));    //8
    //_mesh.insert_outside_convex_hull_2D(Point(-2.0, 0, 0));   //9
    //_mesh.insert_outside_convex_hull_2D(Point(0.0, -2.0, 0)); //10
    //_mesh.insert_outside_convex_hull_2D(Point(0.0, 2.0, 0));  //11
    //_mesh.insert_outside_convex_hull_2D(Point(3.0, 3.0, 0));  //11

    //_mesh.delaunayize_lawson();

    //_mesh.face_split(0, _mesh.barycenter_of_face(0), true);
}

void GeometricWorld::precompute_mesh_curvature()
{
    m_precomputed_curvatures.resize(_mesh.m_vertices.size());

    Vector max_curvature = Vector(std::numeric_limits<double>::min());
    for (int i = 0; i < _mesh.m_vertices.size(); i++)
    {
        Vector curvature = length(_mesh.laplacian_mean_curvature(i));
        max_curvature = max(max_curvature, curvature);

        m_precomputed_curvatures[i] = curvature;
    }

    for (Vector& curvature : m_precomputed_curvatures)
        curvature /= max_curvature;
}

#define COLOR_CURVATURE 0

//Example with a bBox
void GeometricWorld::draw()
{
    for (int i = 0; i < _mesh.m_faces.size(); i++)
    {
        const Face& face = _mesh.m_faces.at(i);
        Point& vertex_a = _mesh.m_vertices.at(face.m_a).get_point();
        Point& vertex_b = _mesh.m_vertices.at(face.m_b).get_point();
        Point& vertex_c = _mesh.m_vertices.at(face.m_c).get_point();

        glBegin(GL_TRIANGLES);
#if COLOR_CURVATURE
        Vector& curvature_a = m_precomputed_curvatures[face.m_a];
        glColor3d(curvature_a.x, curvature_a.y, curvature_a.z);
        glPointDraw(vertex_a);

        Vector& curvature_b = m_precomputed_curvatures[face.m_b];
        glColor3d(curvature_b.x, curvature_b.y, curvature_b.z);
        glPointDraw(vertex_b);

        Vector& curvature_c = m_precomputed_curvatures[face.m_c];
        glColor3d(curvature_c.x, curvature_c.y, curvature_c.z);
        glPointDraw(vertex_c);
#else
        glColor3d(1, 0, 0);
        glPointDraw(vertex_a);

        glColor3d(0, 1, 0);
        glPointDraw(vertex_b);

        glColor3d(0, 0, 1);
        glPointDraw(vertex_c);
#endif

        glEnd();
    }
}

//Example with a wireframe bBox
void GeometricWorld::drawWireFrame() {
    glColor3d(0,1,0);
    glBegin(GL_LINE_STRIP);
    glPointDraw(_bBox[0]);
    glPointDraw(_bBox[1]);
    glEnd();
    glColor3d(0,0,1);
    glBegin(GL_LINE_STRIP);
    glPointDraw(_bBox[0]);
    glPointDraw(_bBox[2]);
    glEnd();
    glColor3d(1,0,0);
    glBegin(GL_LINE_STRIP);
    glPointDraw(_bBox[0]);
    glPointDraw(_bBox[3]);
    glEnd();
}

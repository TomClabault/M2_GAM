#include "offreader.h"

#include <fstream>
#include <iostream>

Mesh OffReader::read_off(const char* filepath)
{
    Mesh mesh;

    std::ifstream file(filepath);

    //Skipping the first "OFF" line
    file.ignore(10000000, '\n');

    //The second line of the file is
    //X Y Z
    //X = nb vertices
    //Y = nb faces
    //Z = nb edges
    int trash;
    int nb_vertices, nb_faces;
    file >> nb_vertices;
    file >> nb_faces;
    file >> trash;//We don't need the number of edges

    mesh.m_vertices.reserve(nb_vertices);
    mesh.m_faces.reserve(nb_faces);

    //Reading all the vertices of the file
    for (int i = 0; i < nb_vertices; i++)
    {
        double x, y, z;

        file >> x;
        file >> y;
        file >> z;

        mesh.m_vertices.push_back(Vertex(-1, Point(x, y, z)));
    }

    //Reading the faces
    for (int i = 0; i < nb_faces; i++)
    {
        int index_first_vertex, index_second_vertex, index_third_vertex;
        int nb_vertices_in_face;

        file >> nb_vertices_in_face;
        if (nb_vertices_in_face != 3)
        {
            std::cout << "The OFF reader can only read triangular faces" << std::endl;

            break;
        }

        file >> index_first_vertex;
        file >> index_second_vertex;
        file >> index_third_vertex;

        //We're setting fa, fb and fc to -1 for now, we'll initialize them later
        mesh.m_faces.push_back(Face(index_first_vertex, index_second_vertex, index_third_vertex, -1, -1, -1));
    }

    //Now that we have read the vertices and the faces, we're going to:
    //
    //For each vertex, set the index of a face it belongs to
    //For each vertex of a face, set the index of the face opposing that vertex

    //This map maps a pair of vertex (an edge) in a face to the index of the current face
    //(the index of the face the edge belongs to) and
    //the index of the vertex opposing the edge within the vertex (local index of the vertex, not global)
    std::map<std::pair<int, int>, std::pair<int, int>> edges_to_face_and_opposing;

    for (int face_index = 0; face_index < mesh.m_faces.size(); face_index++)
    {
        Face& face = mesh.m_faces.at(face_index);

        //We use the current face to set the 'adjacent_face_index" attribute of all the
        //vertices of the face
        mesh.m_vertices.at(face.m_a).set_adjacent_face_index(face_index);
        mesh.m_vertices.at(face.m_b).set_adjacent_face_index(face_index);
        mesh.m_vertices.at(face.m_c).set_adjacent_face_index(face_index);





        //We want the keys of the map (the global indices of the two vertices of the current edge)
        //to be ordered in increasing order
        //i.e. we want that each key of the map is std::pair<int, int>(a, b) so that a < b

        int first_edge_smallest_index = face.m_a < face.m_b ? face.m_a : face.m_b;
        int first_edge_biggest_index = face.m_a < face.m_b ? face.m_b : face.m_a;

        auto key_pair = std::pair<int, int>(first_edge_smallest_index, first_edge_biggest_index);
        auto value_pair = std::pair<int, int>(face_index, 2);

        //If the key already exists, i.e. we already have found this edge in another face
        auto value_found = edges_to_face_and_opposing.find(key_pair);
        if (value_found != edges_to_face_and_opposing.end())
        {
            //If the key exists, we're setting the opposing face of the third vertex of the current face
            //to be the face that was stored in the map
            int other_face_index = value_found->second.first;
            int local_vertex_index_opposing_current_edge_in_other_face = value_found->second.second;

            face.m_fc = other_face_index;
            mesh.m_faces.at(other_face_index).opposing_face(local_vertex_index_opposing_current_edge_in_other_face) = face_index;
        }
        else
            edges_to_face_and_opposing.insert({key_pair, value_pair});





        int second_edge_smallest_index = face.m_b < face.m_c ? face.m_b : face.m_c;
        int second_edge_biggest_index = face.m_b < face.m_c ? face.m_c : face.m_b;

        key_pair = std::pair<int, int>(second_edge_smallest_index, second_edge_biggest_index);
        value_pair = std::pair<int, int>(face_index, 0);
        value_found = edges_to_face_and_opposing.find(key_pair);
        if (value_found != edges_to_face_and_opposing.end())
        {
            int other_face_index = value_found->second.first;
            int local_vertex_index_opposing_current_edge_in_other_face = value_found->second.second;

            face.m_fa = other_face_index;
            mesh.m_faces.at(other_face_index).opposing_face(local_vertex_index_opposing_current_edge_in_other_face) = face_index;
        }
        else
            edges_to_face_and_opposing.insert({key_pair, value_pair});





        int third_edge_smallest_index = face.m_c < face.m_a ? face.m_c : face.m_a;
        int third_edge_biggest_index = face.m_c < face.m_a ? face.m_a : face.m_c;
        key_pair = std::pair<int, int>(third_edge_smallest_index, third_edge_biggest_index);
        value_pair = std::pair<int, int>(face_index, 1);
        value_found = edges_to_face_and_opposing.find(key_pair);
        if (value_found != edges_to_face_and_opposing.end())
        {
            int other_face_index = value_found->second.first;
            int local_vertex_index_opposing_current_edge_in_other_face = value_found->second.second;

            face.m_fb = other_face_index;
            mesh.m_faces.at(other_face_index).opposing_face(local_vertex_index_opposing_current_edge_in_other_face) = face_index;
        }
        else
            edges_to_face_and_opposing.insert({key_pair, value_pair});
    }

    return mesh;
}

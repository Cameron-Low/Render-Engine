#include "Loader.h"

#include "Utils.h"
#include "Material.h"

#include <fstream>
#include <map>
#include <glm/glm.hpp>
#include <glm/gtx/string_cast.hpp>

const char* ws = " \t\n\r\f\v";

// trim from end of string (right)
inline std::string& rtrim(std::string& s, const char* t = ws) {
    s.erase(s.find_last_not_of(t) + 1);
    return s;
}

// trim from beginning of string (left)
inline std::string& ltrim(std::string& s, const char* t = ws) {
    s.erase(0, s.find_first_not_of(t));
    return s;
}

// trim from both ends of string (right then left)
inline std::string& trim(std::string& s, const char* t = ws) {
    return ltrim(rtrim(s, t), t);
}

std::map<std::string, Material> load_obj_mats(const std::string file) {
    std::map<std::string, Material> materials;
    std::string material_name;

    // Parse the obj file line by line
    std::ifstream obj_mat_file("models/" + file + ".mtl");
    std::string line;
    while (std::getline(obj_mat_file, line)) {
        // Split the line
        std::vector<std::string> data = split(trim(line), ' '); 
        if (!data[0].compare("newmtl")) {
            material_name = data[1];
            materials[material_name] = Material(material_name, "");
        } else if (!data[0].compare("Kd")) {
            glm::vec3 colour(255 * std::stof(data[1]), 255 * std::stof(data[2]), 255 * std::stof(data[3]));
            materials[material_name].diffuse_colour = colour;
        } else if (!data[0].compare("Ka")) {
            glm::vec3 colour(255 * std::stof(data[1]), 255 * std::stof(data[2]), 255 * std::stof(data[3]));
            materials[material_name].ambient_colour = colour;
        } else if (!data[0].compare("Ks")) {
            glm::vec3 colour(255 * std::stof(data[1]), 255 * std::stof(data[2]), 255 * std::stof(data[3]));
            materials[material_name].specular_colour = colour;
        } else if (!data[0].compare("Ke")) {
            glm::vec3 colour(255 * std::stof(data[1]), 255 * std::stof(data[2]), 255 * std::stof(data[3]));
            materials[material_name].emission_colour = colour;
        } else if (!data[0].compare("Ns")) {
            materials[material_name].shininess = std::stof(data[1]);
        } else if (!data[0].compare("Ni")) {
            materials[material_name].refractive_index = std::stof(data[1]);
        } else if (!data[0].compare("d")) {
            materials[material_name].transparency = 1 - std::stof(data[1]);
        } else if (!data[0].compare("Tr")) {
            materials[material_name].transparency = std::stof(data[1]);
        } else if (!data[0].compare("Tf")) {
            glm::vec3 colour(255 * std::stof(data[1]), 255 * std::stof(data[2]), 255 * std::stof(data[3]));
            materials[material_name].emission_colour = colour;
        } else if (!data[0].compare("illum")) {
            materials[material_name].model = std::stoi(data[1]);
        } else if (!data[0].compare("map_Kd")) {
            materials[material_name].texture = data[1];
        }
    }
    return materials;
}


std::vector<ModelTriangle> load_obj(const std::string file, float scale){
    std::map<std::string, Material> materials = load_obj_mats(file);
    std::vector<glm::vec3> vertices;
    std::vector<glm::vec3> normals;
    std::vector<glm::vec2> texture_coords;
    std::vector<ModelTriangle> triangles;
    Material material;
    material.diffuse_colour = glm::vec3(255,255,255);

    // Parse the obj file line by line
    std::ifstream obj_file("models/" + file + ".obj");
    std::string line;
    while (std::getline(obj_file, line)) {
        // Split the line
        std::vector<std::string> data = split(trim(line), ' ');
        if (!data[0].compare("usemtl")) {
            // Set the current colour
            material = materials[data[1]];
        } else if (!data[0].compare("v")) { 
            // Create a vertex from the coords and apply any scaling
            glm::vec3 vertex(stof(data[1]), stof(data[2]), stof(data[3]));
            vertex *= scale;

            // Add the vertex to the vertex list
            vertices.push_back(vertex);
        } else if (!data[0].compare("vt")) { 
            // Create a vertex from the coords and apply any scaling
            glm::vec2 coord(stof(data[1]), stof(data[2]));

            // Add the vertex to the vertex list
            texture_coords.push_back(coord);
        } else if (!data[0].compare("vn")) { 
            // Create a vertex from the coords and apply any scaling
            glm::vec3 vertex_normal(stof(data[1]), stof(data[2]), stof(data[3]));

            // Add the vertex to the vertex list
            normals.push_back(vertex_normal);
        } else if (!data[0].compare("f")) { // Split the line
            // Split each set of indices into the image/texture indices
            std::vector<int> indices;
            std::vector<int> indices_tex;
            std::vector<int> indices_normals;
            data.erase(data.begin());
            for (std::string i : data) {
                std::vector<std::string> temp = split(i, '/');
                indices.push_back(std::stoi(temp[0]) - 1);
                if (temp.size() > 1 && temp[1].compare("")) {
                    indices_tex.push_back(std::stoi(temp[1]) - 1);
                }
                if (temp.size() == 3 && temp[2].compare("")) {
                    indices_normals.push_back(std::stoi(temp[2]) - 1);
                }
            }

            // Find the vertices from their respective index and store them in a ModelTriangle
            ModelTriangle triangle;
            for (int i = 0; i < 3; i++) {
                triangle.vertices[i] = vertices[indices[i]];
                if (indices_tex.size() != 0) {
                    triangle.texturePoints[i] = texture_coords[indices_tex[i]];
                }
                if (indices_normals.size() != 0) {
                    triangle.has_vns = true;
                    triangle.vertex_normals[i] = glm::normalize(normals[indices_normals[i]]);
                }
            }
            triangle.material = material;
            // Compute normal
            glm::vec3 edge1 = triangle.vertices[1] - triangle.vertices[0];
            glm::vec3 edge2 = triangle.vertices[2] - triangle.vertices[0];
            triangle.surface_normal = glm::normalize(glm::cross(edge1, edge2)); 

            // Add the triangle to the rest
            triangles.push_back(triangle);
        }
    }
    return triangles;
}



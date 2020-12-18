#pragma once

#include <iostream>
#include <glm/glm.hpp>

struct Material {
        std::string name{};
        std::string texture{};
        glm::vec3 colour{};
        glm::vec3 ambient_colour{};
        glm::vec3 diffuse_colour{};
        glm::vec3 specular_colour{};
        glm::vec3 emission_colour{};
        float transparency{};
        float shininess{};
        float refractive_index{};
        int model{};
	Material();
	Material(std::string name, std::string texture);
};

std::ostream &operator<<(std::ostream &os, const Material &material);

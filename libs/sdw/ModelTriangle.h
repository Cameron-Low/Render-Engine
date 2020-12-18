#pragma once

#include <glm/glm.hpp>
#include <string>
#include <array>
#include "Material.h"

struct ModelTriangle {
	std::array<glm::vec3, 3> vertices{};
	std::array<glm::vec3, 3> vertex_normals{};
	std::array<glm::vec2, 3> texturePoints{};
	Material material{};
	glm::vec3 surface_normal{};
        bool has_vns;

	ModelTriangle();
	ModelTriangle(const glm::vec3 &v0, const glm::vec3 &v1, const glm::vec3 &v2, Material trigMaterial);
	friend std::ostream &operator<<(std::ostream &os, const ModelTriangle &triangle);
};

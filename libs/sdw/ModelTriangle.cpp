#include "ModelTriangle.h"
#include <utility>

ModelTriangle::ModelTriangle() = default;

ModelTriangle::ModelTriangle(const glm::vec3 &v0, const glm::vec3 &v1, const glm::vec3 &v2, Material trigMaterial) :
		vertices({{v0, v1, v2}}), vertex_normals({}), texturePoints(), material(std::move(trigMaterial)), surface_normal(), has_vns(false) {}

std::ostream &operator<<(std::ostream &os, const ModelTriangle &triangle) {
	os << "(" << triangle.vertices[0].x << ", " << triangle.vertices[0].y << ", " << triangle.vertices[0].z << ")\n";
	os << "(" << triangle.vertices[1].x << ", " << triangle.vertices[1].y << ", " << triangle.vertices[1].z << ")\n";
	os << "(" << triangle.vertices[2].x << ", " << triangle.vertices[2].y << ", " << triangle.vertices[2].z << ")\n";
	return os;
}

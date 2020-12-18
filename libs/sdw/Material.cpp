#include "Material.h"
#include <utility>
#include <glm/gtx/string_cast.hpp>

Material::Material() = default;
Material::Material(std::string n, std::string t) :
		name(std::move(n)),
		texture(std::move(t)) {}

std::ostream &operator<<(std::ostream &os, const Material &material) {
	os << "Kd " << glm::to_string(material.diffuse_colour) << std::endl
        << "Ka " << glm::to_string(material.ambient_colour) << std::endl
        << "Ks " << glm::to_string(material.specular_colour) << std::endl
        << "Ke " << glm::to_string(material.emission_colour) << std::endl
        << "Ns " << material.shininess << std::endl
        << "Ni " << material.refractive_index << std::endl
        << "Tr " << material.transparency << std::endl
        << "Illum " << material.model << std::endl;
	return os;
}

#include "RayTriangleIntersection.h"

RayTriangleIntersection::RayTriangleIntersection() = default;
RayTriangleIntersection::RayTriangleIntersection(const glm::vec3 &p, float d, const ModelTriangle &t, size_t index) :
		point(p),
		distance(d),
		triangle(t),
		triangle_index(index) {}

std::ostream &operator<<(std::ostream &os, const RayTriangleIntersection &intersection) {
	os << "Intersection is at [" << intersection.point[0] << "," << intersection.point[1] << "," <<
	   intersection.point[2] << "] on triangle " << intersection.triangle <<
	   " at a distance of " << intersection.distance;
	return os;
}

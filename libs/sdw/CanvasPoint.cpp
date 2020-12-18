#include "CanvasPoint.h"

CanvasPoint::CanvasPoint() :
                colour(255, 255, 255),
		texturePoint(-1, -1) {}

CanvasPoint::CanvasPoint(float xPos, float yPos) :
		x(xPos),
		y(yPos),
                colour(255, 255, 255),
		depth(0.0),
		brightness(1.0),
		texturePoint(-1, -1) {}

CanvasPoint::CanvasPoint(float xPos, float yPos, const glm::vec3 &c) :
		x(xPos),
		y(yPos),
                colour(c),
		depth(0.0),
		brightness(1.0),
		texturePoint(-1, -1) {}

CanvasPoint::CanvasPoint(float xPos, float yPos, const glm::vec3 &c, float pointDepth) :
		x(xPos),
		y(yPos),
                colour(c),
		depth(pointDepth),
		brightness(1.0),
		texturePoint(-1, -1) {}

CanvasPoint::CanvasPoint(float xPos, float yPos, const glm::vec3 &c, float pointDepth, float pointBrightness) :
		x(xPos),
		y(yPos),
                colour(c),
		depth(pointDepth),
		brightness(pointBrightness),
		texturePoint(-1, -1) {}

std::ostream &operator<<(std::ostream &os, const CanvasPoint &point) {
	os << "(" << point.x << ", " << point.y << ", " << point.depth << ") " << point.brightness;
	return os;
}

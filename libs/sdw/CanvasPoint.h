#pragma once

#include <glm/glm.hpp>
#include <iostream>

struct CanvasPoint {
	float x{};
	float y{};
        glm::vec3 colour{};
	float depth{};
	float brightness{};
        glm::vec2 texturePoint{};

	CanvasPoint();
	CanvasPoint(float xPos, float yPos);
	CanvasPoint(float xPos, float yPos, const glm::vec3 &colour);
	CanvasPoint(float xPos, float yPos, const glm::vec3 &colour, float pointDepth);
	CanvasPoint(float xPos, float yPos, const glm::vec3 &colour, float pointDepth, float pointBrightness);

	friend std::ostream &operator<<(std::ostream &os, const CanvasPoint &point);
};

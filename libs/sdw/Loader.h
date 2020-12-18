#pragma once

#include "ModelTriangle.h"

#include <string>
#include <vector>

std::vector<ModelTriangle> load_obj(const std::string file, float scale);

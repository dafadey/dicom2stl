#pragma once

#include <iostream>
#include "geo.h"


void writeVector(std::ostream& s, const vec3& p);

vec3 normal(const triangle& t);

void saveSTL(const char* filename, const geo& g);

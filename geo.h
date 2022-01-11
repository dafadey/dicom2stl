#pragma once

#include <array>
#include <vector>
#include <iostream>

typedef std::array<double,3> vec3;

typedef std::array<vec3*,3> triangle;

struct geo : public std::vector<triangle>
{
  std::vector<vec3> points;
  int memory();
};

vec3 operator+(const vec3& a, const vec3& b);

vec3 operator*(const vec3& a, const double b);

vec3 operator*(const double b, const vec3& a);

vec3 operator-(const vec3& a, const vec3& b);

vec3 operator/(const vec3& a, const double b);

vec3 cross_prod(const vec3& a, const vec3& b);

double operator*(const vec3& a, const vec3& b);

void normalize(vec3& a);

//vec3 max(const vec3& a, const vec3& b);

//vec3 min(const vec3& a, const vec3& b);

std::ostream& operator<<(std::ostream& s, const vec3& v);

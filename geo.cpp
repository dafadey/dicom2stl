#include <cmath>
#include "geo.h"

vec3 operator+(const vec3& a, const vec3& b) {
  return vec3{a[0]+b[0], a[1] + b[1] ,a[2] + b[2]};
}

vec3 operator*(const vec3& a, const double b) {
  return vec3{ a[0] * b, a[1] * b ,a[2] * b };
}

vec3 operator*(const double b, const vec3& a) {
  return vec3{ a[0] * b, a[1] * b ,a[2] * b };
}

vec3 operator-(const vec3& a, const vec3& b) {
  return vec3{ a[0] - b[0], a[1] - b[1] ,a[2] - b[2] };
}

vec3 operator/(const vec3& a, const double b) {
  return vec3{ a[0] / b, a[1] / b ,a[2] / b };
}

vec3 cross_prod(const vec3& a, const vec3& b) {
  return vec3{ a[1]*b[2]-a[2]*b[1] , -a[0] * b[2] + a[2] * b[0], a[0] * b[1] - a[1] * b[0]};
}

double operator*(const vec3& a, const vec3& b) {
  return a[0]*b[0]+ a[1] * b[1]+ a[2] * b[2];
}

void normalize(vec3& a) {
  a = a / std::sqrt(a*a);
}

/*
vec3 max(const vec3& a, const vec3& b) {
  return vec3{std::max(a[0],b[0]),std::max(a[1],b[1]),std::max(a[2],b[2]) };
}

vec3 min(const vec3& a, const vec3& b) {
  return vec3{ std::min(a[0],b[0]),std::min(a[1],b[1]),std::min(a[2],b[2]) };
}
*/

std::ostream& operator<<(std::ostream& s, const vec3& v) {
  return s << '(' << v[0] << ", " << v[1] << ", " << v[2] << ')';
}

int geo::memory() {
  return this->size() * sizeof(vec3*) * 3 + points.size() * sizeof(vec3);
}

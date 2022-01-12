#pragma once
#include "geo.h"
#include <string>
#include <iostream>

typedef std::array<int, 3> ivec;

typedef std::array<float,3> vec3f;

std::ostream& operator<<(std::ostream& o, const vec3f& v);

std::ostream& operator<<(std::ostream& o, const ivec& v);

struct vgeo : public std::vector<float> {
  vgeo() = delete;
  vgeo(int i, int j, int k);
  float min;
  float max;
  const ivec dim;
  const int dim10;
  vec3f d{ 1.f,1.f,1.f };
  vec3f o{ .0f,.0f,.0f };
  float& get(int _i, int _j, int _k);
  int memory();
  std::string toBase64() const;
};

vec3f get_origin_from_box(const vec3f& box_origin, const vgeo& g);

vgeo crop(vgeo& vg, const vec3f& center, const vec3f& scale);

void embed(const vgeo& from, vgeo& to);

vgeo smooth(vgeo& vg, float r0);

void sortPoints(vec3* pts, int n, const vec3& normal);

geo vgeo2surface(const vgeo& vg, float level, float level_side, const vgeo* = nullptr, const vgeo* = nullptr);

vgeo vgeo_from_base64(const std::string& base64str);

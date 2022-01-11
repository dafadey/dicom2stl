#pragma once

#include <array>
#include <vector>

typedef std::array<float, 2> fpt;

struct func : public std::vector<fpt> {
  // returns {xmin, xmax, ymin, ymax}
  std::array<float, 4> get_bounds();
  int segments_count();
  float get_value(float x);
};

fpt operator*(const fpt& a, const fpt& b);
fpt operator/(const fpt& a, const fpt& b);
fpt operator+(const fpt& a, const fpt& b);
fpt operator-(const fpt& a, const fpt& b);
fpt operator*(const fpt& a, const float& b);
fpt operator*(const float b, const fpt& a);
fpt operator/(const fpt& a, const float& b);

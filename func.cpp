#include "func.h"
#include <limits>

// returns {xmin, xmax, ymin, ymax}
std::array<float, 4> func::get_bounds() {
  std::array<float, 4> res{std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), -std::numeric_limits<float>::max()};
  for(auto& pt : *this) {
    for(int i=0; i<2 ; i++) {
      res[0+i*2] = res[0+i*2] < pt[i] ? res[0+i*2] : pt[i];
      res[1+i*2] = res[1+i*2] > pt[i] ? res[1+i*2] : pt[i];
    }
  }
  return res;
}

int func::segments_count() {
  return this->size() - 1;
}

float func::get_value(float x) {
  for (int i = 0; i < this->size() - 1; i++) {
    float x0 = (*this)[i][0];
    float x1 = (*this)[i+1][0];
    if((x0 - x) * (x1 - x) <= .0f) {
      float y0 = (*this)[i][1];
      float y1 = (*this)[i + 1][1];
      return  (y0 * (x1 - x) + y1 * (x - x0)) / (x1-x0);
    }
  }
  return 0;
}


fpt operator*(const fpt& a, const fpt& b) {return fpt{a[0]*b[0], a[1]*b[1]};}
fpt operator/(const fpt& a, const fpt& b) {return fpt{a[0]/b[0], a[1]/b[1]};}
fpt operator+(const fpt& a, const fpt& b) {return fpt{a[0]+b[0], a[1]+b[1]};}
fpt operator-(const fpt& a, const fpt& b) {return fpt{a[0]-b[0], a[1]-b[1]};}
fpt operator*(const fpt& a, const float& b) {return fpt{a[0]*b, a[1]*b};}
fpt operator*(const float b, const fpt& a) {return fpt{a[0]*b, a[1]*b};}
fpt operator/(const fpt& a, const float& b) {return fpt{a[0]/b, a[1]/b};}

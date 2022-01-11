#include <fstream>
#include "saveSTL.h"

void writeVector(std::ostream& s, const vec3& p) {
  for (int i = 0; i < 3; i++) {
    float val = static_cast<float>(p[i]);
    s.write(reinterpret_cast<char*>(&val), sizeof(float));
  }
}

vec3 normal(const triangle& t) {
  return cross_prod(*t[1]-*t[0], *t[2] - *t[0]);
}

void saveSTL(const char* filename, const geo& g) {
  std::ofstream f(filename, std::ios_base::binary);
  for(int i=0;i<80;i++)
    f << 'f';
  int tris_count = static_cast<int>(g.size());
  f.write(reinterpret_cast<char*>(&tris_count), sizeof(int));
  for (const auto& t : g) {
    auto n=normal(t);
    normalize(n);
    writeVector(f, n);
    writeVector(f, *t[0]);
    writeVector(f, *t[1]);
    writeVector(f, *t[2]);
    f << char(0) << char(0);
  }
  f.close();
}

#include <omp.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <cassert>

#include "vgeo.h"
#include "timer.h"

#include "base64/base64.h"

vgeo::vgeo(int i, int j, int k) : dim(ivec{i,j,k}), dim10(i*j), min(std::numeric_limits<float>::max()), max(-std::numeric_limits<float>::max()) {
  resize(i*j*k);
}
  
float& vgeo::get(int _i, int _j, int _k)
{
  int i=_i<0 ? _i+dim[0] : _i;
  i = i>=dim[0] ? i-dim[0] : i;
  int j = _j < 0 ? _j + dim[1] : _j;
  j = j >= dim[1] ? j - dim[1] : j;
  int k = _k < 0 ? _k + dim[2] : _k;
  k = k >= dim[2] ? k - dim[2] : k;
  return (*this)[i+j*dim[0]+k*dim10];
}

vec3f get_origin_from_box(const vec3f& box_center, const vgeo& g) {
  vec3f res;
  for (int i = 0; i != 3; i++)
    res[i] = std::floor((static_cast<float>(g.dim[i]) * .5 - box_center[i] / g.d[i])) * g.d[i];
  return res;
}

vgeo crop(vgeo& vg, const vec3f& center, const vec3f& scale) {
  timer tim("vgeo::crop time is ");
  vgeo res(std::floor(vg.dim[0] * scale[0]), std::floor(vg.dim[1] * scale[1]), std::floor(vg.dim[2] * scale[2]));
  res.d = vg.d;
  res.o = get_origin_from_box(center, res);
  //std::cout << "cropping strcuture of size " << vg.dim[0] << " x " << vg.dim[1] << " x " << vg.dim[2] << ", center: (" << vg.o[0] << ", " << vg.o[1] << ", " << vg.o[2] << ") dim is(" << vg.d[0] << ", " << vg.d[1] << ", " << vg.d[2] << ") -> to structue of size " << res.dim[0] << " x " << res.dim[1] << " x " << res.dim[2] << " center: (" << res.o[0] << ", " << res.o[1] << ", " << res.o[2] << ") dim is (" << res.d[0] << ", " << res.d[1] << ", " << res.d[2] << ")\n";
  #pragma omp parallel for
  for (int kk = 0; kk < res.dim[2]; kk++) {
    int k = std::floor(static_cast<double>(vg.dim[2]) * .5) - std::round(res.o[2] / vg.d[2]) + kk;
    for (int jj = 0; jj < res.dim[1]; jj++) {
      int j = std::floor(static_cast<double>(vg.dim[1]) * .5) - std::round(res.o[1] / vg.d[1]) + jj;
      for (int ii = 0; ii < res.dim[0]; ii++) {
        int i = std::floor(static_cast<double>(vg.dim[0]) * .5) - std::round(res.o[0] / vg.d[0]) + ii;
        res[ii+jj*res.dim[0]+kk*res.dim10] = vg.get(i, j, k);
      }
    }
  }
  return res;
}

void embed(const vgeo& from, vgeo& to) {
  timer tim("vgeo::embed time is ");
  std::cout << "embedding mask " << from.dim[0] << 'x' << from.dim[1] << 'x' << from.dim[2] << " to mask " << to.dim[0] << 'x' << to.dim[1] << 'x' << to.dim[2] << '\n';
  for (int kk = 0; kk < from.dim[2]; kk++) {
    int k = std::floor(static_cast<double>(to.dim[2]) * .5) - std::round(from.o[2] / to.d[2]) + kk;
    for (int jj = 0; jj < from.dim[1]; jj++) {
      int j = std::floor(static_cast<double>(to.dim[1]) * .5) - std::round(from.o[1] / to.d[1]) + jj;
      for (int ii = 0; ii < from.dim[0]; ii++) {
        int i = std::floor(static_cast<double>(to.dim[0]) * .5) - std::round(from.o[0] / to.d[0]) + ii;
        //std::cout << i << ' ' << j << ' ' << k << "<-" << ii << ' ' << jj << ' ' << kk << '\n';
        to.get(i,j,k) = from[ii+jj*from.dim[0]+kk*from.dim10];
      }
    }
  }
}

vgeo smooth(vgeo& vg, float r0) {
  timer tim("vgeo::smooth time is ");
  vgeo res(vg.dim[0], vg.dim[1], vg.dim[2]);
  for (int i = 0; i < 3; i++) {
    res.d[i] = vg.d[i];
    res.o[i] = vg.o[i];
  }
  //exp(-pow(r/r0,2)) = e-2;
  // r/r0 = sqrt(2)
  int tn[3];
  for(int i=0;i<3;i++)
    tn[i] = 2 * std::ceil(r0 * std::sqrt(2.f)/vg.d[i]);
  vgeo tool(tn[0], tn[1], tn[2]);
  for (int i = 0; i < 3; i++)
    tool.d[i] = vg.d[i];

  float r02 = r0*r0;

  std::cout << "tool size is " << tool.dim[0] << 'x' << tool.dim[1] << 'x' << tool.dim[2] << '\n';

  for (int kk = 0; kk < tn[2]; kk++) {
    for (int jj = 0; jj < tn[1]; jj++) {
      for (int ii = 0; ii < tn[0]; ii++) {
        float r2=std::pow(static_cast<float>(ii-tn[0]/2)/tool.d[0],2) + std::pow(static_cast<float>(jj - tn[1] / 2) / tool.d[1], 2)+ std::pow(static_cast<float>(kk - tn[2] / 2) / tool.d[2], 2);
        tool[ii+jj*tool.dim[0]+ kk * tool.dim10] = std::exp(-r2/r02);
      }
    }
  }

  std::cout << "smoothing...\n";
  #pragma omp parallel for
  for (int k = 0; k < vg.dim[2]; k++) {
    //std::cout << "smoothing " << static_cast<double>(k * 1000 / vg.dim[2]) / 10. << "%\n";
    for (int j = 0; j < vg.dim[1]; j++) {
      for (int i = 0; i < vg.dim[0]; i++) {
        float val=.0f;
        float norm=.0f;
        for (int kk = 0; kk < tn[2]; kk++) {
          for (int jj = 0; jj < tn[1]; jj++) {
            for (int ii = 0; ii < tn[0]; ii++) {
              val+= vg.get(i+ii-tn[0]/2,j+jj-tn[1]/2,k+kk-tn[2]/2) * tool[ii+ jj * tool.dim[0] + kk * tool.dim10];
              norm += tool[ii + jj * tool.dim[0] + kk * tool.dim10];
            }
          }
        }
        val/=norm;
        res[i+j*res.dim[0]+k*res.dim10] = val;
      }
    }
  }
  
  return res;
}

bool criteria(float ay, float ax, float by, float bx) {
  // instead of atan2(ax,ay) > atan2(bx,by])
  if(ay * by < 0)
    return ay > 0;
  if(ax * bx < 0)
    return ay > 0 ? ax < 0 : ax > 0;
  if(ax * bx > 0)
    return ay > 0 ? ax < bx : ax > bx;
  
  return std::atan2(ay,ax) > std::atan2(by,bx);
}

void sortPoints(vec3* pts, int n, const vec3& normal) {
  vec3 center{.0f,.0f,.0f};
  for(int i=0; i<n;i++)
    center = center + pts[i];
  center = center / static_cast<float>(n);
  vec3 u{.0f,.0f,.0f};
  for(int i=0;i<3;i++) {
    vec3 cand = cross_prod(normal, vec3{i == 0 ? 1.f : .0f, i == 1 ? 1.f : .0f, i == 2 ? 1.f : .0f });
    u = u*u < cand*cand ? cand : u;
  }
  vec3 v = cross_prod(u, normal);
//  std::sort(pts, pts+n, [&v, &u, &center](const vec3& a, const vec3& b){ return std::atan2((a-center) * u, (a - center) * v) > std::atan2((b - center) * u, (b - center) * v);});
  std::sort(pts, pts+n, [&v, &u, &center](const vec3& a, const vec3& b){ return criteria((a-center) * u, (a - center) * v, (b - center) * u, (b - center) * v);});

}

geo vgeo2surface(const vgeo& vg, float level, float level_side, const vgeo* level_vgeo, const vgeo* mask_vgeo) {
  timer tim("vgeo::vgeo2sufrace time is ");

  std::vector<float> xvalues((vg.dim[0] - 1) * vg.dim[1] * vg.dim[2], std::numeric_limits<float>::max());
  std::vector<float> yvalues(vg.dim[0] * (vg.dim[1] - 1) * vg.dim[2], std::numeric_limits<float>::max());
  std::vector<float> zvalues(vg.dim[0] * vg.dim[1] * (vg.dim[2] - 1), std::numeric_limits<float>::max());

  std::vector<int> boxes(vg.dim[0] * vg.dim[1] * vg.dim[2], 0);

  if(!level_vgeo) {
    //std::cout << "doing xvalues\n";
    #pragma omp parallel for
    for (int k = 0; k < vg.dim[2]; k++) {
      for (int j = 0; j < vg.dim[1]; j++) {
          float y = static_cast<float>(j) * vg.d[1] - vg.o[1];
          float ys = std::abs(2. * y / static_cast<float>(vg.dim[1]));
          for (int i = 0; i < vg.dim[0] - 1; i++) {
            float x0 = static_cast<float>(i) * vg.d[0] - vg.o[0];
            float x1 = static_cast<float>(i + 1) * vg.d[0] - vg.o[0];
            float xs0 = std::abs(2. * x0 / static_cast<float>(vg.dim[0]));
            float xs1 = std::abs(2. * x1 / static_cast<float>(vg.dim[0]));
            const auto v0 = vg[i + j * vg.dim[0] + k * vg.dim10] - (level * (1. - xs0) + level_side * xs0);
            const auto v1 = vg[i+1 + j * vg.dim[0] + k * vg.dim10] - (level * (1. - xs1) + level_side * xs1);
            if(v1*v0 < 0.f) {
              xvalues[i+j*(vg.dim[0]-1)+ k * (vg.dim[0] - 1)* vg.dim[1]] = (x1 * v0 - x0 * v1) / (v0 - v1);
              boxes[i+j*vg.dim[0]+k*vg.dim10] = 1;
              boxes[i+(j?j-1:0)*vg.dim[0]+k*vg.dim10] = 1;
              boxes[i+j*vg.dim[0]+(k?k-1:k)*vg.dim10] = 1;
              boxes[i+(j?j-1:0)*vg.dim[0]+(k?k-1:k)*vg.dim10] = 1;
            }
        }
      }
    }
    //std::cout << "doing yvalues\n";
    #pragma omp parallel for
    for (int k = 0; k < vg.dim[2]; k++) {
      for (int j = 0; j < vg.dim[1] - 1; j++) {
        float y = static_cast<float>(j) * vg.d[1] - vg.o[1];
        float ys = std::abs(2. * y / static_cast<float>(vg.dim[1]));
        for (int i = 0; i < vg.dim[0]; i++) {
          float x = static_cast<float>(i) * vg.d[0] - vg.o[0];
          float xs = std::abs(2. * x / static_cast<float>(vg.dim[0]));
          const auto v0 = vg[i + j * vg.dim[0] + k * vg.dim10] - (level * (1 - xs) + level_side * xs);
          const auto v1 = vg[i + (j+1) * vg.dim[0] + k * vg.dim10] - (level * (1 - xs) + level_side * xs);
          if (v1 * v0 < 0.f) {
            float y0 = static_cast<float>(j) * vg.d[1] - vg.o[1];
            float y1 = static_cast<float>(j + 1) * vg.d[1] - vg.o[1];
            yvalues[i + j * vg.dim[0] + k * vg.dim[0] * (vg.dim[1]-1)] = (y1 * v0 - y0 * v1) / (v0 - v1);
            boxes[i+j*vg.dim[0]+k*vg.dim10] = 1;
            boxes[(i?i-1:i)+j*vg.dim[0]+k*vg.dim10] = 1;
            boxes[i+j*vg.dim[0]+(k?k-1:k)*vg.dim10] = 1;
            boxes[(i?i-1:i)+j*vg.dim[0]+(k?k-1:k)*vg.dim10] = 1;
          }
        }
      }
    }
    //std::cout << "doing zvalues\n";
    #pragma omp parallel for
    for (int k = 0; k < vg.dim[2] - 1; k++) {
      for (int j = 0; j < vg.dim[1]; j++) {
        float y = static_cast<float>(j) * vg.d[1] - vg.o[1];
        float ys = std::abs(2. * y / static_cast<float>(vg.dim[1]));
        for (int i = 0; i < vg.dim[0]; i++) {
          float x = static_cast<float>(i) * vg.d[0] - vg.o[0];
          float xs = std::abs(2. * x / static_cast<float>(vg.dim[0]));
          const auto v0 = vg[i + j * vg.dim[0] + k * vg.dim10] - (level * (1 - xs) + level_side * xs);
          const auto v1 = vg[i + j * vg.dim[0] + (k + 1) * vg.dim10] - (level * (1 - xs) + level_side * xs);
          if (v1 * v0 < 0.f) {
            float z0 = static_cast<float>(k) * vg.d[2] - vg.o[2];
            float z1 = static_cast<float>(k + 1) * vg.d[2] - vg.o[2];
            zvalues[i + j * vg.dim[0] + k * vg.dim10] = (z1 * v0 - z0 * v1) / (v0 - v1);
            boxes[i+j*vg.dim[0]+k*vg.dim10] = 1;
            boxes[(i?i-1:i)+j*vg.dim[0]+k*vg.dim10] = 1;
            boxes[i+(j?j-1:j)*vg.dim[0]+k*vg.dim10] = 1;
            boxes[(i?i-1:i)+(j?j-1:j)*vg.dim[0]+k*vg.dim10] = 1;
          }
        }
      }
    }
  } else { //3d level
    //std::cout << "iso surface with 3d level\n";
    //std::cout << "doing xvalues\n";
    #pragma omp parallel for
    for (int k = 0; k < vg.dim[2]; k++) {
      for (int j = 0; j < vg.dim[1]; j++) {
        for (int i = 0; i < vg.dim[0] - 1; i++) {
          float x0 = static_cast<float>(i) * vg.d[0] - vg.o[0];
          float x1 = static_cast<float>(i + 1) * vg.d[0] - vg.o[0];
          const auto v0 = vg[i + j * vg.dim[0] + k * vg.dim10] - (*level_vgeo)[i + j * vg.dim[0] + k * vg.dim10];
          const auto v1 = vg[i + 1 + j * vg.dim[0] + k * vg.dim10] - (*level_vgeo)[i + 1 + j * vg.dim[0] + k * vg.dim10];
          if (v1 * v0 < 0.f) {
            xvalues[i + j * (vg.dim[0] - 1) + k * (vg.dim[0] - 1) * vg.dim[1]] = (x1 * v0 - x0 * v1) / (v0 - v1);
            boxes[i+j*vg.dim[0]+k*vg.dim10] = 1;
            boxes[i+(j?j-1:0)*vg.dim[0]+k*vg.dim10] = 1;
            boxes[i+j*vg.dim[0]+(k?k-1:k)*vg.dim10] = 1;
            boxes[i+(j?j-1:0)*vg.dim[0]+(k?k-1:k)*vg.dim10] = 1;
          }
        }
      }
    }
    //std::cout << "doing yvalues\n";
    #pragma omp parallel for
    for (int k = 0; k < vg.dim[2]; k++) {
      for (int j = 0; j < vg.dim[1] - 1; j++) {
        for (int i = 0; i < vg.dim[0]; i++) {
          const auto v0 = vg[i + j * vg.dim[0] + k * vg.dim10] - (*level_vgeo)[i + j * vg.dim[0] + k * vg.dim10];
          const auto v1 = vg[i + (j + 1) * vg.dim[0] + k * vg.dim10] - (*level_vgeo)[i + (j + 1) * vg.dim[0] + k * vg.dim10];
          if (v1 * v0 < 0.f) {
            float y0 = static_cast<float>(j) * vg.d[1] - vg.o[1];
            float y1 = static_cast<float>(j + 1) * vg.d[1] - vg.o[1];
            yvalues[i + j * vg.dim[0] + k * vg.dim[0] * (vg.dim[1] - 1)] = (y1 * v0 - y0 * v1) / (v0 - v1);
            boxes[i+j*vg.dim[0]+k*vg.dim10] = 1;
            boxes[(i?i-1:i)+j*vg.dim[0]+k*vg.dim10] = 1;
            boxes[i+j*vg.dim[0]+(k?k-1:k)*vg.dim10] = 1;
            boxes[(i?i-1:i)+j*vg.dim[0]+(k?k-1:k)*vg.dim10] = 1;
          }
        }
      }
    }
    //std::cout << "doing zvalues\n";
    #pragma omp parallel for
    for (int k = 0; k < vg.dim[2] - 1; k++) {
      for (int j = 0; j < vg.dim[1]; j++) {
        for (int i = 0; i < vg.dim[0]; i++) {
          const auto v0 = vg[i + j * vg.dim[0] + k * vg.dim10] - (*level_vgeo)[i + j * vg.dim[0] + k * vg.dim10];
          const auto v1 = vg[i + j * vg.dim[0] + (k + 1) * vg.dim10] - (*level_vgeo)[i + j * vg.dim[0] + (k + 1) * vg.dim10];
          if (v1 * v0 < 0.f) {
            float z0 = static_cast<float>(k) * vg.d[2] - vg.o[2];
            float z1 = static_cast<float>(k + 1) * vg.d[2] - vg.o[2];
            zvalues[i + j * vg.dim[0] + k * vg.dim10] = (z1 * v0 - z0 * v1) / (v0 - v1);
            boxes[i+j*vg.dim[0]+k*vg.dim10] = 1;
            boxes[(i?i-1:i)+j*vg.dim[0]+k*vg.dim10] = 1;
            boxes[i+(j?j-1:j)*vg.dim[0]+k*vg.dim10] = 1;
            boxes[(i?i-1:i)+(j?j-1:j)*vg.dim[0]+k*vg.dim10] = 1;
          }
        }
      }
    }
  }

  std::cout << "vgeo::vgeo2sufrace elapsed" << tim.elapsed() << " s\n";

  int boxes_in=0;
  for(auto& b : boxes)
    boxes_in = b ? boxes_in + 1 : boxes_in;
  std::cout << static_cast<float>(boxes_in)/static_cast<float>(boxes.size()) << " % of boxes are filled\n";

  const int emap[4][2]{{0,0}, {0,1}, {1,0}, {1,1}};

  int ngeos{};
  #pragma omp parallel
  {
    ngeos = omp_get_num_threads();
  }

  //std::cout << "ngeos is " << ngeos << '\n';
  std::vector<geo> geos(ngeos);
  //geo g;
  //std::cout << "loop over boxes\n";
  #pragma omp parallel for schedule(dynamic, 1)
  for (int k = 0; k < vg.dim[2] - 1; k++) {
    vec3 normal;
    std::array<vec3, 12> pts;
    const int threadid = omp_get_thread_num();
    //std::cout << static_cast<double>(k*1000/vg.dim[2])/10. << "%\n";
    for (int j = 0; j < vg.dim[1] - 1; j++) {
      for (int i = 0; i < vg.dim[0] - 1; i++) {
        if(!boxes[i+j*vg.dim[0]+k*vg.dim10])
          continue;
        if(mask_vgeo && (*mask_vgeo)[i+j*vg.dim[0]+k*vg.dim10] == 0.f)
          continue;
        
        float y[2]{ static_cast<float>(j) * vg.d[1] - vg.o[1], static_cast<float>(j + 1) * vg.d[1] - vg.o[1] };
        float z[2]{ static_cast<float>(k) * vg.d[2] - vg.o[2], static_cast<float>(k + 1) * vg.d[2] - vg.o[2] };
        float x[2]{ static_cast<float>(i) * vg.d[0] - vg.o[0], static_cast<float>(i + 1) * vg.d[0] - vg.o[0] };
        
        int eig=0;
        for(int ei=0; ei<4; ei++) {
          float val = xvalues[i + (j + emap[ei][0]) * (vg.dim[0] - 1) + (k + emap[ei][1]) * (vg.dim[0] - 1) * vg.dim[1]];
          if(val != std::numeric_limits<float>::max()) {
            pts[eig] = vec3{val, y[emap[ei][0]], z[emap[ei][1]]};
            eig++;
          }
        }
        for (int ei = 0; ei < 4; ei++) {
          float val = yvalues[(i + emap[ei][0]) + j * vg.dim[0] + (k + emap[ei][1]) * vg.dim[0] * (vg.dim[1] - 1)];
          if (val != std::numeric_limits<float>::max()) {
            pts[eig] = vec3{ x[emap[ei][0]], val, z[emap[ei][1]] };
            eig++;
          }
        }
        for (int ei = 0; ei < 4; ei++) {
          float val = zvalues[(i + emap[ei][0]) + (j + emap[ei][1]) * vg.dim[0] + k * vg.dim10];
          if (val != std::numeric_limits<float>::max()) {
            pts[eig] = vec3{ x[emap[ei][0]], y[emap[ei][1]], val };
            eig++;
          }
        }

        if (eig < 3)
          continue;
        
        normal = vec3{(vg[i+1+j*vg.dim[0]+k*vg.dim10]-vg[i + j * vg.dim[0] + k * vg.dim10])/vg.d[0], (vg[i + (j+1) * vg.dim[0] + k * vg.dim10] - vg[i + j * vg.dim[0] + k * vg.dim10]) / vg.d[1], (vg[i + j * vg.dim[0] + (k+1) * vg.dim10] - vg[i + j * vg.dim[0] + k * vg.dim10]) / vg.d[2] };

        normalize(normal);
        sortPoints(pts.data(), eig, normal);
          
        vec3* adr = static_cast<vec3*>(nullptr) + geos[threadid].points.size();
        
        for (int i = 0; i < eig; ++i)
          geos[threadid].points.push_back(pts[i]);

        for (int i = 2; i < eig; ++i) {
          int i1=i-1;
          geos[threadid].push_back(triangle{adr+i,adr+i1,adr});
        }
      }
    }
  }
  
  //std::cout << "collecting to one surface\n";
  
  size_t points_size{};
  for(const auto& pg : geos) {
    //std::cout << "pg.points.size()=" << pg.points.size() << '\n';
    points_size += pg.points.size();
  }

  geo g;
  
  g.points.resize(points_size);
  const size_t adr0 = reinterpret_cast<size_t>(g.points.data());

  size_t pt_id = 0;
  for(const auto& pg : geos) {
    size_t pt_id0 = pt_id;
    for(const auto& pt : pg.points) {
      g.points[pt_id] = pt;
      pt_id++;
    }
    for(const triangle& t : pg) {
      triangle new_t;
      for(int i=0; i < 3; i++) {
        size_t adr = reinterpret_cast<size_t>(t[i]) + pt_id0 * sizeof(vec3);
        new_t[i] = reinterpret_cast<vec3*>(adr0 + adr);
      }
      g.push_back(new_t);
    }
  }
  
  return g;
}

int vgeo::memory() {
  return this->size() * sizeof(float);
}

static std::vector<float> compress_vgeo(const vgeo* vg) {
  std::vector<float> out;
  for(int i=0;i<3;i++) {
    out.push_back(vg->dim[i]);
    out.push_back(1.f);
  }
  float v = (*vg)[0];
  int count = 1;
  for (int i=1;i<vg->size();i++) {
    if((*vg)[i] == v)
      count++;
    else {
      out.push_back(v);
      out.push_back(static_cast<float>(count));
      v = (*vg)[i];
      count = 1;
    }
  }
  out.push_back(v);
  out.push_back(static_cast<float>(count));
  return out;
}

static vgeo decompress_vgeo(const std::vector<float>& compressed) {
  int dim[3];
  for(int i=0; i < 3; i++)
    dim[i] = static_cast<int>(compressed[i*2]);
  vgeo out(dim[0], dim[1], dim[2]);
  int id = 0;
  for(int i=6;i<compressed.size();i+=2) {
    for(int j=0;j<static_cast<int>(compressed[i+1]);j++, id++)
      out[id] = compressed[i];
  }
  return out;
}

std::string vgeo::toBase64() const {
  std::vector<float> compressed = compress_vgeo(this);
  return base64::encode(reinterpret_cast<BYTE const*>(compressed.data()), compressed.size() * sizeof(float));
}

vgeo vgeo_from_base64(const std::string& base64str) {
  std::vector<BYTE> compressed_byte = base64::decode(base64str);
  std::vector<float> compressed;
  const int n = compressed_byte.size() / sizeof(float);
  for(int i=0;i<n;i++)
    compressed.push_back(reinterpret_cast<float*>(compressed_byte.data())[i]);
  return decompress_vgeo(compressed);
}

std::ostream& operator<<(std::ostream& o, const vec3f& v) {
  return o << '(' << v[0] << ", " << v[1] << ", " << v[2] << ')';
}

std::ostream& operator<<(std::ostream& o, const ivec& v) {
  return o << '[' << v[0] << "x" << v[1] << "x" << v[2] << ']';
}

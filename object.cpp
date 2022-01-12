#include <fstream>
#include <omp.h>

#include "object.h"
#include "readDICOM.h"
#include "saveSTL.h"
#include "draw.h"
#include "timer.h"

object::object(vgeo* source_) : source(source_), smoothed(nullptr), cropped(nullptr), mask(nullptr), mask_cropped(nullptr), surface_nomask(nullptr), surface(nullptr), box_center{.0f, .0f, .0f}, box_scale{1.f, 1.f, 1.f}, smooth_factor(1.7), level{633} {
  for(int i=0; i< levels.size(); i++) {
    auto& lvl = levels[i];
    lvl.push_back(fpt{-source->o[i], 0.f});
    lvl.push_back(fpt{ source->d[i] * static_cast<float>(source->dim[i]) - source->o[i], 0.f});
  }
  mask = new vgeo(source->dim[0], source->dim[1], source->dim[2]);
  mask->o = source->o;
  mask->d = source->d;
  for(auto& v : *mask)
    v=1.f;
}

object::~object()
{
  if (cropped)
    delete cropped;
  if (smoothed)
    delete smoothed;
  if(surface_nomask)
    delete surface_nomask;
  if(surface)
    delete surface;
  if(mask)
    delete mask;
  if (mask_cropped)
    delete mask_cropped;
}

void object::reset_mask() {
  std::cout << "resetting mask\n";
  for (auto& v : *mask_cropped)
    v = 1.f;
}

void object::invert_mask() {
  for (auto& v : *mask_cropped)
    v = v == 0. ? 1.f : .0f;
}

void object::crop_surface() {
  timer tim("object::crop_surface time is ");
  if(cropped)
    delete cropped;
  
  cropped = new vgeo(crop(*smoothed, box_center, box_scale));
  
  if (mask_cropped) {
    embed(*mask_cropped, *mask);
    delete mask_cropped;
  }
  mask_cropped = new vgeo(crop(*mask, box_center, box_scale));
}

void object::smooth_surface() {
  timer tim("object::smooth_surface time is ");
  if(smoothed)
    delete smoothed;
  smoothed = new vgeo(smooth(*source, smooth_factor));
}

void object::make_surface(renderer* ren) {
  timer tim("object::make_surface time is ");
  
  if(!cropped)
    return;
  
  bool first_time = false;
  if(surface_nomask)
    delete surface_nomask;
  else
    first_time = true;

  vgeo level_vg(cropped->dim[0], cropped->dim[1], cropped->dim[2]);
  level_vg.o = cropped->o;
  level_vg.d = cropped->d;
  std::array<std::vector<float>,3> levels_discretized;
  for(int lid=0; lid <levels_discretized.size(); lid++) {
    levels_discretized[lid].resize(cropped->dim[lid], .0f);
    for (int i = 0; i < levels_discretized[lid].size(); i++)
      levels_discretized[lid][i] = levels[lid].get_value(static_cast<float>(i) * level_vg.d[lid] - level_vg.o[lid]);
  }
  
  #pragma omp parallel for
  for (int k = 0; k < level_vg.dim[2]; k++) {
    float z = static_cast<float>(k) * level_vg.d[2] - level_vg.o[2];
    for (int j = 0; j < level_vg.dim[1]; j++) {
      float y = static_cast<float>(j) * level_vg.d[1] - level_vg.o[1];
      for (int i = 0; i < level_vg.dim[0]; i++) {
        float x = static_cast<float>(i) * level_vg.d[0] - level_vg.o[0];
        level_vg[i + j * level_vg.dim[0] + k * level_vg.dim10] = level + levels_discretized[0][i] + levels_discretized[1][j] + levels_discretized[2][k];
      }
    }
  }
  
  surface_nomask = new geo(vgeo2surface(*cropped, level, level, &level_vg));
  //surface = new geo(vgeo2surface(*cropped, level, level, &level_vg, mask_cropped));
  if(first_time) {
    ogl_surface.init(ren);
    ren->items.push_back(&ogl_surface);
  }
  apply_mask();
  
  ogl_surface.update_model(surface);
}

void object::apply_mask() {
  timer tim("object::apply_mask time is ");
  //std::cout << "apply mask : " << mask_cropped->dim[0] << 'x' << mask_cropped->dim[1] << 'x' << mask_cropped->dim[2] << " size=" << mask_cropped->size() << '\n';
  if(surface)
    delete surface;
  surface = new geo;
  for(auto& pt : surface_nomask->points)
    surface->points.push_back(pt);

  int ngeos{};
  #pragma omp parallel
  {
    ngeos = omp_get_num_threads();
  }
  std::vector<geo> geos(ngeos);
  #pragma omp parallel for
  for(int t_id = 0; t_id < surface_nomask->size(); t_id++) {
    const int threadid = omp_get_thread_num();
    auto& t = (*surface_nomask)[t_id];
    //coord_i = mask_icoord_i * mask.d[i] - mask.o[i]
    bool skip_triangle = false;
    for(int t_id=0;t_id<3;t_id++) {
      int icoords[3];
      for(int i=0;i<3;i++)
        icoords[i] = std::floor(static_cast<float>((*(t[t_id]))[i] + mask_cropped->o[i]) / mask_cropped->d[i]);
      if(mask_cropped->get(icoords[0], icoords[1], icoords[2]) == 0.f) {
        skip_triangle = true;
        break;
      }
    }
    
    if(skip_triangle)
      continue;
      
    triangle new_t;
    for(int i=0;i<3;i++)
      new_t[i] = t[i] - surface_nomask->points.data() + surface->points.data();
    geos[threadid].push_back(new_t);
  }
  
  for(auto& g : geos) {
    for(auto& t : g)
      surface->push_back(t);
  }
}

void object::save_surface() {
  saveSTL((output_name+".stl").c_str(), *surface);
  std::ofstream ini((output_name + ".ini").c_str());
  ini << "name = \"" << name << "\"\n";
  ini << "source_file = \"" << source_name << "\"\n";
  ini << "level = " << level << '\n';
  ini << "box_center = (" << box_center[0] << ", " << box_center[1] << ", " << box_center[2] << ")\n";
  ini << "box_scale = (" << box_scale[0] << ", " << box_scale[1] << ", " << box_scale[2] << ")\n";
  ini << "smoothness = " << smooth_factor << '\n';
  ini << "x_levels = ";
  for(const auto& pt : levels[0])
    ini << " (" << pt[0] << ", " << pt[1] << ')';
  ini << '\n';
  ini << "y_levels = ";
  for (const auto& pt : levels[1])
    ini << " (" << pt[0] << ", " << pt[1] << ')';
  ini << '\n';
  ini << "z_levels = ";
  for (const auto& pt : levels[2])
    ini << " (" << pt[0] << ", " << pt[1] << ')';
  ini << '\n'; 
  if(mask) {
    if(mask_cropped)
      embed(*mask_cropped, *mask);
    ini << "mask = " << mask->toBase64() << "\n";
  }
  ini << '\n';
}

static std::string remove_chars(const std::string& input, char c) {
  std::string out;
  for (size_t i = 0; i < input.size(); i++) {
    if (input.c_str()[i] == c)
      continue;
    out += input.c_str()[i];
  }
  return out;
}

static std::vector<std::string> split(const std::string& input, const std::string& delimiter) {
  size_t pos = 0;
  std::vector<std::string> out;
  while (true) {
    size_t end = input.find(delimiter, pos);
    out.push_back(input.substr(pos, end - pos));
    pos = end + delimiter.size();
    if (end == std::string::npos)
      break;
  }
  return out;
}


static vec3f read_vector(const std::string& input) {
  vec3f output;
  auto items = split(remove_chars(remove_chars(input, '('), ')'), ",");
  for (int i = 0; i < 3; i++)
    output[i] = std::stof(items[i]);
  return output;
}

static func read_func(const std::string& input) {
  func f;
  auto points = split(remove_chars(input, '('), ")");
  for (auto p : points) {
    if(p.empty())
      continue;
    auto coords = split(p,",");
    f.push_back(fpt{std::stof(coords[0]), std::stof(coords[1])});
  }
  return f;
}

static std::string leave_text_inside_quotes(const std::string& input) {
  std::string out;
  bool start = false;
  for (size_t i = 0; i < input.size(); i++) {
    if (input.c_str()[i] == '\"') {
      start = true;
      continue;
    }
    if (input.c_str()[i] == '\"' && start)
      break;
    if(start)
      out += input.c_str()[i];
  }
  return out;
}

object* load_surface(const std::string& ini_name, std::vector<vgeo>& vgs) {
  std::ifstream ini(ini_name.c_str());
  std::string line;
  std::string source_file;
  vgeo* vg{nullptr};
  if(vgs.size())
    vg = &vgs[0];
  
  float level;
  std::array<func, 3> levels;
  vec3f box_center;
  vec3f box_scale;
  float smooth_factor;
  std::string name;
  std::string source_name;
  vgeo* mask{};
  while (std::getline(ini, line)) {
    auto tokens = split(line, "=");
    if(tokens.size() < 2)
      continue;
    auto first_token = remove_chars(tokens[0], ' ');
    auto second_token = remove_chars(tokens[1], ' ');
    if (first_token == "source_file" && !vg) {
      source_name = leave_text_inside_quotes(tokens[1]);
      std::cout << "READ INI: source_name = " << source_name << '\n';
      DICOMreader r(source_name.c_str(), "");
      r.read(vgs);
      vg = &vgs[0];
    }
    if (first_token == "name")
      name = leave_text_inside_quotes(tokens[1]);
    if (first_token == "level")
      level = std::stof(second_token);
    if (first_token == "smoothness")
      smooth_factor = std::stof(second_token);
    if (first_token == "box_center")
        box_center = read_vector(second_token);
    if (first_token == "box_scale")
        box_scale = read_vector(second_token);
    if (first_token == "x_levels")
      levels[0] = read_func(second_token);
    if (first_token == "y_levels")
      levels[1] = read_func(second_token);
    if (first_token == "z_levels")
      levels[2] = read_func(second_token);
    if (first_token == "mask")
      mask=new vgeo(vgeo_from_base64(second_token));
  }

  ini.close();
  object* obj = new object(vg);
  if(obj->mask)
    delete obj->mask;
  obj->mask = mask;
  obj->mask->d = vg->d;
  obj->mask->o = vg->o;
  //obj->mask->o = get_origin_from_box(box_center, *mask);
  obj->level = level;
  obj->box_center = box_center;
  obj->box_scale = box_scale;
  obj->levels = levels;
  obj->name = name;
  obj->source_name = source_name;
  obj->smooth_factor = smooth_factor;
  std::cout << "read mask with o at " << obj->mask->o[0] << ',' << obj->mask->o[1] << ',' << obj->mask->o[2] << ", dim " << obj->mask->dim[0] << ',' << obj->mask->dim[1] << ',' << obj->mask->dim[2] << ", d " << obj->mask->d[0] << ','<< obj->mask->d[1] << ',' << obj->mask->d[2] << ", size: " << obj->mask->size() << '\n';
  return obj;
}

void object::draw_mask(float x0, float y0, float radius_x, float radius_y, const float* view_mat, const float* proj_mat, bool erase) {
  timer tim("object::draw_mask time is ");
  /*
  v_k = vm_i,k  * v_i

  v_j = pm_k,j * v_k
  
  v_j = pm_k,j * vm_i,k v_i
  
  m_i,j = pm_k,j * vm_i, k
  */
  float matrix[16];
  for (int j = 0; j < 4; j++) {
    for(int i = 0; i < 4; i++) {
      float sum = .0f;
      for(int k = 0; k < 4; k++)
        sum += proj_mat[k*4 + j] * view_mat[i*4 + k];
      matrix[i * 4 + j] = sum;
    }
  }

  const float radius_x_1 = 1.f / radius_x;
  const float radius_y_1 = 1.f / radius_y;

  int removed = 0;

  #pragma omp parallel for schedule(dynamic, 1)
  for (int k = 0; k < mask_cropped->dim[2]; k++) {
    const float z = mask_cropped->d[2] * static_cast<float>(k) - mask_cropped->o[2];
    for (int j = 0; j < mask_cropped->dim[1]; j++) {
      const float y = mask_cropped->d[1] * static_cast<float>(j) - mask_cropped->o[1];
      const int base_addr = j * mask_cropped->dim[0] + k * mask_cropped->dim10;
      float* it = mask_cropped->data() + base_addr;
      for (int i = 0; i < mask_cropped->dim[0]; i++, it++) {
        if(*it == erase ? .0f : 1.f)
          continue;
        const float x = mask_cropped->d[0] * static_cast<float>(i) - mask_cropped->o[0];
        const float xs = x * matrix[0] + y * matrix[4] + z * matrix[8] + matrix[12];
        const float ys = x * matrix[1] + y * matrix[5] + z * matrix[9] + matrix[13];
        #define SQR(X) (X*X)
        const float dist2 = SQR((xs - x0) * radius_x_1) + SQR((ys - y0) * radius_y_1);
        #undef SQR
        if(dist2 < 1.f) {
          *it = erase ? 0.f : 1.f;
          removed++;
        }
      }
    }
  }
/*
  if(removed) {
    apply_mask();
    update_model();
  }
*/
}

void object::update_model() {
  ogl_surface.update_model(surface);
}

int object::memory() {
  int res{};
  if(source)
    res += source->memory();
  if(smoothed)
    res += smoothed->memory();
  if(cropped)
    res += cropped->memory();
  if(mask)
    res += mask->memory();
  if (mask_cropped)
    res += mask_cropped->memory();
  if(surface_nomask)
    res += surface_nomask->memory();
  if(surface)
    res += surface->memory();
  res += ogl_surface.memory();
  return res;
}

#pragma once

#include "vgeo.h"
#include "OGLitem.h"
#include "func.h"

struct renderer;
struct vgeo;
struct geo;

struct object {
  ~object();
  object() = delete;
  object(vgeo* source_);
  
  void crop_surface();
  void smooth_surface();
  void make_surface(renderer* ren);
  void save_surface();
  void reset_mask();
  void invert_mask();
  void apply_mask();
  void update_model();

  void draw_mask(float x, float y, float radius_x, float radius_y, const float* view_mat, const float* proj_mat);
  
  vgeo* source;
  vgeo* smoothed;
  vgeo* cropped;
  vgeo* mask_cropped;
  vgeo* mask;
  geo* surface_nomask;
  geo* surface;
  OGLsurface ogl_surface;
  std::array<func,3> levels;
  float level;
  vec3f box_center;
  vec3f box_scale;
  float smooth_factor;
  std::string name;
  std::string output_name;
  std::string source_name;
  int memory();
};

object* load_surface(const std::string&, std::vector<vgeo>& vgs);

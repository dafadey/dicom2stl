#pragma once

#ifdef _WIN32
#include <windows.h>
#endif

#include <GL/glew.h>
#include <array>
#include <vector>
#include <cfloat>
#include <string>
#include <map>
#include "geo.h"
#include "OGLitem.h"

class GLFWwindow;

struct object;

struct renderer {
  enum e_mouse_state{HOVER, ROTATE, PAN, ZOOM, SELECT, DRAWMASK};
  
  ~renderer();
  GLuint vao_cursor2d{};
  GLuint VBO_cursor2d{};
  GLuint shader_cursor2d{};
  GLint cursor2d_pos_location{};
  GLint cursor2d_color_location{};
  GLint cursor2d_resolution_location{};
  GLint cursor2d_radii_location{};

  renderer() = default;
  
  std::vector<OGLitem*> items;
  
  std::vector<object*>* objects;
  
  std::map<std::array<std::string,3>, GLuint> shaders;
  
  vec3 cam_pos{};
  vec3 fp_pos{};
  vec3 cam_up{};
  vec3 light_dir{};

  std::array<GLfloat,16> view_matrix;
  std::array<GLfloat,16> proj_matrix;
 
  float mask_r_x{.1f};
  float mask_r_y{.1f};

  GLFWwindow* win;
  
  e_mouse_state mouse_state{HOVER};
  double mouse_pos[2];
  
  bool nocallbacks{};

  bool init(GLFWwindow*, std::vector<object*>*);
  
  void render();

  void set_callbacks(GLFWwindow* window);
  
  GLuint getShader(const std::string& vs_name, const std::string& fs_name, const std::string& gs_name = "");

  void reset_camera();

  void remove_item(OGLitem*);
  
  void draw_mask(double x, double y);
  
  void apply_mask();

  bool parallel=true;

};


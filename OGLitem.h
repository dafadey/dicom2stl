#pragma once

#include <GL/glew.h>
#include <vector>
#include <cfloat>

#include "geo.h"

struct renderer;

struct geo;

struct OGLitem {
  OGLitem();
  virtual void init(renderer* ren) = 0;
  virtual void draw(GLfloat* view_matrix, GLfloat* proj_matrix, GLfloat* light_dir) = 0;
  virtual std::array<vec3,2> get_bounds() = 0;

  virtual ~OGLitem();

  GLuint VBO{};
  GLuint iBO{};
  GLuint vao{};
  GLuint shader{};

  std::vector<GLfloat> VBOdata;
  bool visible = true;
  
  virtual int memory() = 0;
};

struct OGLsurface : public OGLitem {
  virtual void init(renderer* ren) override;
  virtual void draw(GLfloat* view_matrix, GLfloat* proj_matrix, GLfloat* light_dir) override;
  virtual std::array<vec3,2> get_bounds() override;
  void update_model(geo* surface);

  GLint verts_location{};
  GLint norms_location{};
  GLint color_location{};
  GLint light_location{};
  GLint shiny_location{};
  GLint proj_matrix_location{};
  GLint view_matrix_location{};
  vec3 geo_min{ DBL_MAX, DBL_MAX, DBL_MAX };
  vec3 geo_max{ -DBL_MAX, -DBL_MAX, -DBL_MAX };
  
  GLfloat surface_color[3]{.5f, .5f, .5f};
  GLfloat shiny{0.3};
  
  GLint tris_count{};
  
  virtual int memory() override;

};

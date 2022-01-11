#include <iostream>
#include <omp.h>

#include "OGLitem.h"
#include "draw.h"
#include "timer.h"

OGLitem::OGLitem() : vao(0), VBO(0), iBO(0), shader(0) {}

OGLitem::~OGLitem() {
  if(VBO)
    glDeleteBuffers(1, &VBO);
  if(iBO)
    glDeleteBuffers(1, &iBO);
  if(vao)
    glDeleteVertexArrays(1, &vao);
}

void OGLsurface::init(renderer* ren) {
  if(!vao) {
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);
  } else
    std::cout << "WARNING: VAO is already initialized to " << vao << '\n';
  
  if(!VBO)
    glGenBuffers(1, &VBO);
  else
    std::cout << "WARNING: VBO is already initialized to " << VBO << '\n';
  
  if(!iBO)
    glGenBuffers(1, &iBO);
  else
    std::cout << "WARNING: VBO is already initialized to " << iBO << '\n';

  glBindBuffer(GL_ARRAY_BUFFER, VBO);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, iBO);

  shader = ren->getShader("sha.vs", "sha.fs");
  
  verts_location = glGetAttribLocation(shader, "vertex_pos");
  norms_location = glGetAttribLocation(shader, "normal_pos");

  light_location = glGetUniformLocation(shader, "light_dir");
  color_location = glGetUniformLocation(shader, "surface_color");
  proj_matrix_location = glGetUniformLocation(shader, "proj_matrix");
  view_matrix_location = glGetUniformLocation(shader, "view_matrix");
  shiny_location = glGetUniformLocation(shader, "shiny");


  std::cout << "\tattr: verts_location=" << verts_location << '\n';
  std::cout << "\tattr: norms_location=" << norms_location << '\n';
  
  std::cout << "\tuniform: color_location=" << color_location << '\n';
  std::cout << "\tuniform: light_location=" << light_location << '\n';
  std::cout << "\tuniform: shiny_location=" << shiny_location << '\n';
  std::cout << "\tuniform: proj_matrix_location=" << proj_matrix_location << '\n';
  std::cout << "\tuniform: view_matrix_location=" << view_matrix_location << '\n';
  glBindVertexArray(0);
}



void OGLsurface::update_model(geo* g)
{
  timer tim("OGLsurface::update_model time is");
  if(VBO == 0 || iBO == 0) {
    std::cerr << "cannot update model VBO or iBO are not generated, call init() for renderer\n";
    return;
  }
  if(!g) {
    std::cerr << "the model is not set\n";
    return;
  }
  
  glBindVertexArray(vao);

  glBindBuffer(GL_ARRAY_BUFFER, VBO);
  VBOdata.clear();
  VBOdata.resize(g->points.size() * 3 * 2);
  geo_min = vec3{ DBL_MAX, DBL_MAX, DBL_MAX };
  geo_max = vec3{ -DBL_MAX, -DBL_MAX, -DBL_MAX };
  for(int i=0; i<g->points.size(); i++) {
    const vec3& pt = g->points[i];
    for(int c=0; c<3; c++) {
      VBOdata[i*3+c] = pt[c];
      geo_min[c] = geo_min[c] < pt[c] ? geo_min[c] : pt[c];
      geo_max[c] = geo_max[c] > pt[c] ? geo_max[c] : pt[c];
    }
  }

  //std::cout << "geo_min: " << geo_min << ", geo_max: " << geo_max << '\n';

  tris_count = g->size();
  
  std::vector<GLuint> indices(g->size()*3);
  size_t normals_offset = g->points.size() * 3;
  #pragma omp parallel for
  for(int i=0; i<g->size(); i++) {
    vec3 a = *(*g)[i][1] - *(*g)[i][0];
    vec3 b = *(*g)[i][2] - *(*g)[i][0];
    vec3 n = cross_prod(a, b);
    normalize(n);
    for(int vid=0; vid<3; vid++) {
      int id = (reinterpret_cast<size_t>((*g)[i][vid]) - reinterpret_cast<size_t>(g->points.data()))/sizeof(vec3);
      indices[i * 3 + vid] = id;
      for(int c=0; c < 3; c++)
        VBOdata[normals_offset + id * 3 + c] = n[c];
    }
  }

  glBufferData(GL_ARRAY_BUFFER, VBOdata.size() * sizeof(GLfloat), VBOdata.data(), GL_STATIC_DRAW);

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, iBO);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size()*sizeof(GLuint), indices.data(), GL_STATIC_DRAW);

  glVertexAttribPointer(verts_location, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (void*) 0);
  glVertexAttribPointer(norms_location, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (void*) (normals_offset * sizeof(GLfloat)));
  glEnableVertexAttribArray(verts_location);
  glEnableVertexAttribArray(norms_location);

  glBindVertexArray(0);
}

void OGLsurface::draw(GLfloat* view_matrix, GLfloat* proj_matrix, GLfloat* light_dir) {
  if(!visible)
    return;

  glUseProgram(shader);

  glEnable(GL_DEPTH_TEST);
  
  glUniformMatrix4fv(view_matrix_location, 1, false, view_matrix);
  glUniformMatrix4fv(proj_matrix_location, 1, false, proj_matrix);

  glUniform3fv(color_location, 1, surface_color);
  glUniform3fv(light_location, 1, light_dir);
  glUniform1f(shiny_location, shiny);

  glBindVertexArray(vao);

  glDrawElements(GL_TRIANGLES, tris_count * 3, GL_UNSIGNED_INT, 0);
  glBindVertexArray(0);
  glUseProgram(0);
}

std::array<vec3,2> OGLsurface::get_bounds() {
  return std::array<vec3,2>{geo_min, geo_max};
}

int OGLsurface::memory() {
  return VBOdata.size() * sizeof(GLfloat);
}

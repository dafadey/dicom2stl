#include "draw.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

#include <GLFW/glfw3.h>

#include "object.h"

static GLuint loadShaders(const std::string& VertexShaderCode, const std::string& FragmentShaderCode, const std::string& GeometryShaderCode = ""){

    GLuint VertexShaderID = glCreateShader(GL_VERTEX_SHADER);
    GLuint FragmentShaderID = glCreateShader(GL_FRAGMENT_SHADER);
    GLuint GeometryShaderID = 0;
    if(!GeometryShaderCode.empty())
      GeometryShaderID = glCreateShader(GL_GEOMETRY_SHADER);

    GLint Result = GL_FALSE;
    int InfoLogLength;

    std::cout << "building vertex shader...\n";
    char const * VertexSourcePointer = VertexShaderCode.c_str();
    glShaderSource(VertexShaderID, 1, &VertexSourcePointer , NULL);
    glCompileShader(VertexShaderID);

    glGetShaderiv(VertexShaderID, GL_COMPILE_STATUS, &Result);
    glGetShaderiv(VertexShaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);
    if ( InfoLogLength > 0 ){
      std::vector<char> VertexShaderErrorMessage(InfoLogLength+1);
      glGetShaderInfoLog(VertexShaderID, InfoLogLength, NULL, &VertexShaderErrorMessage[0]);
      std::cerr << "building vertex shader failed with " << VertexShaderErrorMessage.data() << '\n';
    }

    if(GeometryShaderID) {
      std::cout << "building geometry shader...\n";
      char const* GeometrySourcePointer = GeometryShaderCode.c_str();
      glShaderSource(GeometryShaderID, 1, &GeometrySourcePointer, NULL);
      glCompileShader(GeometryShaderID);

      glGetShaderiv(GeometryShaderID, GL_COMPILE_STATUS, &Result);
      glGetShaderiv(GeometryShaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);
      if (InfoLogLength > 0) {
        std::vector<char> GeometryShaderErrorMessage(InfoLogLength + 1);
        glGetShaderInfoLog(GeometryShaderID, InfoLogLength, NULL, &GeometryShaderErrorMessage[0]);
        std::cerr << "building geometry shader failed with " << GeometryShaderErrorMessage.data() << '\n';
      }
    }

    std::cout << "building fragment shader...\n";
    char const * FragmentSourcePointer = FragmentShaderCode.c_str();
    glShaderSource(FragmentShaderID, 1, &FragmentSourcePointer , NULL);
    glCompileShader(FragmentShaderID);

    glGetShaderiv(FragmentShaderID, GL_COMPILE_STATUS, &Result);
    glGetShaderiv(FragmentShaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);
    if ( InfoLogLength > 0 ){
      std::vector<char> FragmentShaderErrorMessage(InfoLogLength+1);
      glGetShaderInfoLog(FragmentShaderID, InfoLogLength, NULL, &FragmentShaderErrorMessage[0]);
      std::cerr << "building fragment shader failed with " << FragmentShaderErrorMessage.data() << '\n';
    }

    std::cout << "building shader program...\n";
    GLuint ProgramID = glCreateProgram();
    glAttachShader(ProgramID, VertexShaderID);
    if (GeometryShaderID)
      glAttachShader(ProgramID, GeometryShaderID);
    glAttachShader(ProgramID, FragmentShaderID);
    glLinkProgram(ProgramID);

    glGetProgramiv(ProgramID, GL_LINK_STATUS, &Result);
    glGetProgramiv(ProgramID, GL_INFO_LOG_LENGTH, &InfoLogLength);
    if ( InfoLogLength > 0 ){
      std::vector<char> ProgramErrorMessage(InfoLogLength+1);
      glGetProgramInfoLog(ProgramID, InfoLogLength, NULL, &ProgramErrorMessage[0]);
      std::cerr << "building shader program failed with " << ProgramErrorMessage.data() << '\n';
    }

    glDeleteShader(VertexShaderID);
    glDeleteShader(GeometryShaderID);
    glDeleteShader(FragmentShaderID);

    return ProgramID;
}

static void cursor_position_callback(GLFWwindow* window, double xpos, double ypos) {

  renderer* ren = (renderer*) glfwGetWindowUserPointer(window);
  if (ren->nocallbacks)
    return;

  if(ren->mouse_state == renderer::e_mouse_state::PAN) {
    vec3 cam_y = ren->cam_up;
    normalize(cam_y);
    vec3 cam_z = ren->cam_pos - ren->fp_pos;
    normalize(cam_z);
    vec3 cam_x = cross_prod(cam_y, cam_z);
    ren->cam_pos = ren->cam_pos - cam_x * (xpos-ren->mouse_pos[0]) + cam_y * (ypos-ren->mouse_pos[1]);
    ren->fp_pos = ren->fp_pos - cam_x * (xpos-ren->mouse_pos[0]) + cam_y * (ypos-ren->mouse_pos[1]);
  }

  if(ren->mouse_state == renderer::e_mouse_state::ROTATE) {
    vec3 cam_y = ren->cam_up;
    normalize(cam_y);
    vec3 cam_z = ren->cam_pos - ren->fp_pos;
    normalize(cam_z);
    vec3 cam_x = cross_prod(cam_y, cam_z);
    vec3 cam_pos_rel = ren->cam_pos - ren->fp_pos;

    cam_pos_rel = vec3{cam_pos_rel * cam_x, cam_pos_rel * cam_y, cam_pos_rel * cam_z};
    
    double phi = -(xpos-ren->mouse_pos[0])/113.;
    double theta = (ypos-ren->mouse_pos[1])/113.;

    vec3 cam_up_rel{-std::sin(phi) * std::sin(theta), std::cos(theta), -std::cos(phi) * std::sin(theta)};
    ren->cam_up = cam_up_rel[0] * cam_x + cam_up_rel[1] * cam_y + cam_up_rel[2] * cam_z;
    normalize(ren->cam_up);

    cam_pos_rel = vec3{cam_pos_rel[2] * std::sin(phi) * std::cos(theta), cam_pos_rel[2] * std::sin(theta), cam_pos_rel[2] * std::cos(phi) * std::cos(theta)};
    
    ren->cam_pos = cam_pos_rel[0] * cam_x + cam_pos_rel[1] * cam_y + cam_pos_rel[2] * cam_z;
    ren->cam_pos = cross_prod(ren->cam_up, cross_prod(ren->cam_pos, ren->cam_up));
    ren->cam_pos = ren->cam_pos + ren->fp_pos;

  }
  
  if(ren->mouse_state == renderer::e_mouse_state::DRAWMASK) {
    int win_geo[2];
    glfwGetWindowSize(window, win_geo, &win_geo[1]);
    ren->draw_mask((xpos - win_geo[0]/2)/static_cast<double>(win_geo[0]) * 2.f, (win_geo[1]/2 - ypos)/static_cast<double>(win_geo[1]) * 2.f);
  }

  ren->mouse_pos[0] = xpos;
  ren->mouse_pos[1] = ypos;
}

static void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)  {
  renderer* ren = (renderer*) glfwGetWindowUserPointer(window);
  if (ren->nocallbacks)
    return;

  vec3 distance = ren->cam_pos - ren->fp_pos;
  
  GLfloat factor = std::pow(1.1, -yoffset);

  distance = distance * factor;
  
  ren->cam_pos = distance + ren->fp_pos;

}

static void mouse_button_callback(GLFWwindow* window, int button, int action, int mods) {
  renderer* r = (renderer*)glfwGetWindowUserPointer(window);
  if (r->nocallbacks)
    return;
  if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS)
    r->mouse_state = renderer::e_mouse_state::ROTATE;

  if (button == GLFW_MOUSE_BUTTON_MIDDLE && action == GLFW_PRESS)
    r->mouse_state = renderer::e_mouse_state::PAN;

  if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS) {
    r->mouse_state = renderer::e_mouse_state::DRAWMASK;
    double xpos, ypos;
    glfwGetCursorPos(window, &xpos, &ypos);
    cursor_position_callback(window, xpos, ypos);
  }

  if (action == GLFW_RELEASE) {
    if (r->mouse_state == renderer::e_mouse_state::DRAWMASK)
      r->apply_mask();
    r->mouse_state = renderer::e_mouse_state::HOVER;
  }

}

static void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
  renderer* r = (renderer*)glfwGetWindowUserPointer(window);
  glViewport(0, 0, width, height);
}

bool renderer::init(GLFWwindow* win_, std::vector<object*>* _objects) {
  GLint GLM{};
  GLint GLm{};
  glGetIntegerv(GL_MAJOR_VERSION, &GLM);
  glGetIntegerv(GL_MINOR_VERSION, &GLm);
  std::cout << "GL:" << GLM << '.' << GLm << '\n';

  objects = _objects;

  win = win_;

  auto glew_err = glewInit();
  if (glew_err != GLEW_OK) {
    std::cerr << "renderer::init: ERROR: failed init glew with error " << glewGetErrorString(glew_err) << '\n';
    return false;
  }

  //make identity matrices
  for(int i=0; i<4; i++) {
    for(int j=0; j<4; j++) {
      view_matrix[i * 4 + j] = i == j ? 1.f : 0.f;
      proj_matrix[i * 4 + j] = i == j ? 1.f : 0.f;
    }
  }

  if (!vao_cursor2d) {
    glGenVertexArrays(1, &vao_cursor2d);
    glBindVertexArray(vao_cursor2d);
  }
  else
    std::cout << "WARNING: VAO for cursord2d is already initialized to " << vao_cursor2d << '\n';

  if (!VBO_cursor2d)
    glGenBuffers(1, &VBO_cursor2d);
  else
    std::cout << "WARNING: VBO for cursor2d is already initialized to " << VBO_cursor2d << '\n';

  glBindBuffer(GL_ARRAY_BUFFER, VBO_cursor2d);

  shader_cursor2d = getShader("sha_circle.vs", "sha_circle.fs", "sha_circle.gs");

  cursor2d_pos_location = glGetAttribLocation(shader_cursor2d, "circle_pos");
  cursor2d_color_location = glGetAttribLocation(shader_cursor2d, "circle_color");
  cursor2d_resolution_location = glGetAttribLocation(shader_cursor2d, "circle_resolution");
  cursor2d_radii_location = glGetAttribLocation(shader_cursor2d, "circle_radii");

  std::cout << "\tattr: cursor2d pos_location=" << cursor2d_pos_location << '\n';
  std::cout << "\tattr: cursor2d color_location=" << cursor2d_color_location << '\n';
  std::cout << "\tattr: cursor2d resolution_location=" << cursor2d_resolution_location << '\n';
  std::cout << "\tattr: cursor2d radii_location=" << cursor2d_radii_location << '\n';

  glBindVertexArray(0);

  return true;
}

renderer::~renderer() {
  if (VBO_cursor2d)
    glDeleteBuffers(1, &VBO_cursor2d);
  if (vao_cursor2d)
    glDeleteVertexArrays(1, &vao_cursor2d);
}

void renderer::set_callbacks(GLFWwindow* window) {
  glfwSetCursorPosCallback(window, cursor_position_callback);
  glfwSetMouseButtonCallback(window, mouse_button_callback);
  glfwSetScrollCallback(window, scroll_callback);
  glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
}

void renderer::render() {
  int win_geo[2];
  glfwGetWindowSize(win, win_geo, &win_geo[1]);

  vec3 cam_y = cam_up;
  normalize(cam_y);
  vec3 cam_z = cam_pos - fp_pos;
  normalize(cam_z);
  vec3 cam_x = cross_prod(cam_y, cam_z);

  //x_new = (p - cam_pos) * cam_x;
  //y_new = (p - cam_pos) * cam_y;
  //z_new = (p - cam_pos) * cam_z;
  
  vec3 shift = vec3{cam_pos * cam_x, cam_pos * cam_y, cam_pos * cam_z};
  
  for(int i=0;i<3;i++) {
    view_matrix[i*4] = cam_x[i];
    view_matrix[i*4+1] = cam_y[i];
    view_matrix[i*4+2] = cam_z[i];
  }
  for(int i=0;i<3;i++)
    view_matrix[3*4+i] = - shift[i];

  light_dir = cam_z + cam_x * 0.23 + cam_y * 0.33;
  normalize(light_dir);
 
  GLfloat f = FLT_MAX;
  GLfloat n = -FLT_MAX;
  
  for(auto& item : items) {
    auto item_bounds = item->get_bounds();
    for(int i=0; i<8; i++) {
      vec3 outline_vertex;
      for(int c=0; c<3; c++)
        outline_vertex[c] = item_bounds[(i>>c) & 1][c];
      GLfloat proj = static_cast<GLfloat>((outline_vertex - cam_pos) * cam_z);
      f = f < proj ? f : proj;
      n = n > proj ? n : proj;
    }
  }

  GLfloat scale;
  
  if(parallel)
    scale = 4.f * std::tan(56.f / 2.f * 3.14f / 180.f) / std::sqrt((cam_pos - fp_pos)*(cam_pos - fp_pos));
  else
    scale = 1.f/std::tan(56.f / 2.f * 3.14f / 180.f); //https://www.scratchapixel.com/lessons/3d-basic-rendering/perspective-and-orthographic-projection-matrix/building-basic-perspective-projection-matrix

  
  /* 
   *  0  1  2  3
   *  4  5  6  7
   *  8  9 10 11
   * 12 13 14 15
   * they multiply vector by COLUMN of matrix! yes, they do!
   */
  
  if(parallel) {
    proj_matrix[0] = scale;
    proj_matrix[5] = scale/static_cast<GLfloat>(win_geo[1]) * static_cast<GLfloat>(win_geo[0]);
    proj_matrix[10] = 2/(f-n);
    proj_matrix[14] = -(f+n)/(f-n);
    proj_matrix[11] = .0f;
    proj_matrix[15] = 1.f;
  } else { // perspective
    proj_matrix[0] = scale;
    proj_matrix[5] = scale/static_cast<GLfloat>(win_geo[1]) * static_cast<GLfloat>(win_geo[0]);
    proj_matrix[10] = -f/(f-n);
    proj_matrix[11] = -1.f;
    proj_matrix[14] = f*n/(f-n);
    proj_matrix[15] = 0.f;
  }

  GLfloat light[3]{static_cast<GLfloat>(light_dir[0]), static_cast<GLfloat>(light_dir[1]), static_cast<GLfloat>(light_dir[2])};

  for(auto& item : items)
    item->draw(view_matrix.data(), proj_matrix.data(), light);

  // draw mask brush
  // https://open.gl/geometry
  glBindVertexArray(vao_cursor2d);

  glBindBuffer(GL_ARRAY_BUFFER, VBO_cursor2d);

  float cursor2d_pos_x = (mouse_pos[0] - win_geo[0] / 2) / static_cast<double>(win_geo[0]) * 2.f;
  float cursor2d_pos_y = (win_geo[1] / 2 - mouse_pos[1]) / static_cast<double>(win_geo[1]) * 2.f;
  float VBOdata[8] { cursor2d_pos_x, cursor2d_pos_y, mask_erases ? 1.f : 0.f, mask_erases ? .0f : .7f, .0f, 64.0, mask_r_x, mask_r_y / static_cast<GLfloat>(win_geo[1]) * static_cast<GLfloat>(win_geo[0])};
  glBufferData(GL_ARRAY_BUFFER, 8 * sizeof(GLfloat), VBOdata, GL_STATIC_DRAW);
  glVertexAttribPointer(cursor2d_pos_location, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(GLfloat), (void*)0);
  glVertexAttribPointer(cursor2d_color_location, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (void*)(2 * sizeof(GLfloat)));
  glVertexAttribPointer(cursor2d_resolution_location, 1, GL_FLOAT, GL_FALSE, 1 * sizeof(GLfloat), (void*)(5 * sizeof(GLfloat)));
  glVertexAttribPointer(cursor2d_radii_location, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(GLfloat), (void*)(6 * sizeof(GLfloat)));

  glEnableVertexAttribArray(cursor2d_pos_location);
  glEnableVertexAttribArray(cursor2d_color_location);
  glEnableVertexAttribArray(cursor2d_resolution_location);
  glEnableVertexAttribArray(cursor2d_radii_location);

  glBindVertexArray(0);
  glUseProgram(shader_cursor2d);

  glDisable(GL_DEPTH_TEST);

  glBindVertexArray(vao_cursor2d);

  glDrawArrays(GL_POINTS, 0, 1);
  glBindVertexArray(0);
  glUseProgram(0);
}

GLuint renderer::getShader(const std::string& vs_name, const std::string& fs_name, const std::string& gs_name) {
  GLuint shader;
  auto it = shaders.find(std::array<std::string,3>{vs_name, fs_name, gs_name});
  if(it == shaders.end()) {
    std::stringstream vs;
    std::stringstream gs;
    std::stringstream fs;
    std::ifstream ivs(vs_name.c_str());
    vs << ivs.rdbuf();
    ivs.close();
    if(!gs_name.empty()) {
      std::ifstream igs(gs_name.c_str());
      gs << igs.rdbuf();
      igs.close();
    }
    std::ifstream ifs(fs_name);
    fs << ifs.rdbuf();
    ifs.close();

    shader = loadShaders(vs.str(), fs.str(), gs.str());
    if(!shader) {
      std::cerr << "load shaders failed\n";
      return false;
    }
    shaders[std::array<std::string,3>{vs_name, fs_name, gs_name}] = shader;
    return shader;
  } else
    return it->second;
}

void renderer::reset_camera() {
  vec3 geo_d{.0, .0, .0};
  fp_pos = vec3{.0, .0, .0};
  for(auto& item : items) {
    auto item_bounds = item->get_bounds();
    
    for(int c=0; c<3; c++) {
      double d = item_bounds[1][c] - item_bounds[0][c];
      geo_d[c] = geo_d[c] > d ? geo_d[c] : d;
    }
    fp_pos = fp_pos + (item_bounds[1] + item_bounds[0])*.5;
  }
  fp_pos = fp_pos / static_cast<double>(items.size());

  cam_pos = fp_pos + vec3{0.f, 0.f, 3*geo_d[2]};

  cam_up = vec3{0.f,1.f,0.f};

}

void renderer::remove_item(OGLitem* item_to_remove) {
  for(int i=0;i<items.size(); i++) {
    if(items[i] == item_to_remove) {
      items.erase(items.begin()+i);
      break;
    }
  }
}

void renderer::draw_mask(double x, double y) {
  int win_geo[2];
  glfwGetWindowSize(win, win_geo, &win_geo[1]);
  for(auto& obj : *objects)
    obj->draw_mask(static_cast<float>(x), static_cast<float>(y), mask_r_x, mask_r_y / static_cast<GLfloat>(win_geo[1]) * static_cast<GLfloat>(win_geo[0]), view_matrix.data(), proj_matrix.data(), mask_erases);
}

void renderer::apply_mask() {
  for (auto& obj : *objects) {
    obj->apply_mask();
    obj->update_model();
  }
}

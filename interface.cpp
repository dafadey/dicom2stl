#include "interface.h"
#include <iostream>

static void glfw_error_callback(int error, const char* description)
{
    std::cerr << "Glfw Error " << error << ": " << description << '\n' << std::flush;
}

bool imgui_interface::init() {

  glfwSetErrorCallback(glfw_error_callback);
  
 
  if (!glfwInit()) {
    std::cerr << "interface::init: ERROR: failed init glfw\n";
    return false;
  }
      
  // GL 3.2 + GLSL 150
  const char* glsl_version = "#version 150";
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);  // 3.2+ only
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);            // Required on Mac
  
  // Create window with graphics context
  window = glfwCreateWindow(1024, 768, "dicom2stl", NULL, NULL);
  if (window == NULL) {
    std::cerr << "interface::init: ERROR: failed to create glfw window\n";
    return false;
  }
  glfwMakeContextCurrent(window);
  glfwSwapInterval(1); // Enable vsync

  IMGUI_CHECKVERSION();
  ImGui::CreateContext();
  #ifndef NOIMPLOT
  ImPlot::CreateContext();
  #endif

  // Setup Dear ImGui style
  //ImGui::StyleColorsDark();
  //ImGui::StyleColorsClassic();
  ImGui::StyleColorsLight();
  ImGui::GetStyle().FrameBorderSize=1;
  ImGui::GetStyle().FrameRounding=3;

  // Setup Platform/Renderer backends
  ImGui_ImplGlfw_InitForOpenGL(window, true);
  ImGui_ImplOpenGL3_Init(glsl_version);

  return true;
}

void imgui_interface::close() {
  // Cleanup
  ImGui_ImplOpenGL3_Shutdown();
  ImGui_ImplGlfw_Shutdown();
  #ifndef NOIMPLOT
  ImPlot::DestroyContext();
  #endif
  ImGui::DestroyContext();

  glfwDestroyWindow(window);
  glfwTerminate();
}

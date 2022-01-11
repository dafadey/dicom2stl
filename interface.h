#pragma once

#include <GLFW/glfw3.h>
#include <imgui.h>
#include <imgui_internal.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>
#ifndef NOIMPLOT
#include <implot.h>
#endif

struct imgui_interface {

  GLFWwindow* window{};
  
  bool init();
  void close();

};

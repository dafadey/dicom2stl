#pragma once
#include <string>

#include <imgui_internal.h>
#include <imgui.h>
#include "func.h"


namespace ImGui {
  ImVec2 operator*(const ImVec2& a, const ImVec2& b);
  ImVec2 operator/(const ImVec2& a, const ImVec2& b);
  ImVec2 operator+(const ImVec2& a, const ImVec2& b);
  ImVec2 operator-(const ImVec2& a, const ImVec2& b);
  ImVec2 operator*(const ImVec2& a, const float& b);
  ImVec2 operator*(const float& b, const ImVec2& a);
  ImVec2 operator/(const ImVec2& a, const float& b);
 
  bool DrawPlot(func& f, const std::string& name, const ImColor& color = ImColor(0,0,0), const std::array<float, 2>* active_bounds=nullptr);
}

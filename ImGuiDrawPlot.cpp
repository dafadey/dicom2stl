#include "ImGuiDrawPlot.h"
#include <map>
#include <iostream>

namespace ImGui {
  ImVec2 operator*(const ImVec2& a, const ImVec2& b) {return ImVec2(a.x*b.x, a.y*b.y);}
  ImVec2 operator/(const ImVec2& a, const ImVec2& b) {return ImVec2(a.x/b.x, a.y/b.y);}
  ImVec2 operator+(const ImVec2& a, const ImVec2& b) {return ImVec2(a.x+b.x, a.y+b.y);}
  ImVec2 operator-(const ImVec2& a, const ImVec2& b) {return ImVec2(a.x-b.x, a.y-b.y);}
  ImVec2 operator*(const ImVec2& a, const float& b) {return ImVec2(a.x*b, a.y*b);}
  ImVec2 operator*(const float& b, const ImVec2& a) {return ImVec2(a.x*b, a.y*b);}
  ImVec2 operator/(const ImVec2& a, const float& b) {return ImVec2(a.x/b, a.y/b);}

  static ImVec2 mapPoint(const fpt& pt, const ImRect& fbb, const ImRect& bb) {
    return (ImVec2(pt[0], pt[1]) - fbb.Min) / (fbb.Max - fbb.Min) * (bb.Max - bb.Min) + bb.Min;
  }
  
  static bool in3pix(ImVec2 v1, ImVec2 v2) {
    return std::abs(v1.x-v2.x)<13 && std::abs(v1.y-v2.y)<13;
  }
  
  struct DrawPlotState {
    DrawPlotState() = default;
    enum eState{HOVER, DRAGFPT, DRAGBOUNDS};
    int call_count{0};
    eState state = eState::HOVER;
    ImRect fbb;
    int pt_id;
  };
  
  bool DrawPlot(func& f, const std::string& name, const ImColor& color, const std::array<float, 2>* active_bounds) {
    static std::map<std::string, DrawPlotState> drawPlotStatesMap;
    if(drawPlotStatesMap.find(name) == drawPlotStatesMap.end())
      drawPlotStatesMap[name] = DrawPlotState();
    
    auto& drawPlotState = drawPlotStatesMap[name];
    
    //std::cout << name << ':' << drawPlotState.call_count << ' ' << drawPlotState.pt_id << ' ' << (drawPlotState.state == DrawPlotState::eState::DRAGFPT ? "drag" : "hover") << '\n';
    
    if(f.size() <= 1)
      return false;
    const ImGuiStyle& Style = GetStyle();
    const ImGuiIO& IO = GetIO();
    ImDrawList* drawList = GetWindowDrawList();
    ImGuiWindow* window = GetCurrentWindow();
    if (window->SkipItems)
      return false;
    
    // prepare canvas
    const float avail = GetContentRegionAvail().x;
    const float dim = ImMin(avail, 128.f);
    ImVec2 canvas(avail, dim);

    ImRect bb(window->DC.CursorPos, window->DC.CursorPos + canvas);
    
    //std::cout << name << " : " << bb.Min.x << ' ' << bb.Min.y << " - " << bb.Max.x << ' ' << bb.Max.y << '\n';
    
    ItemSize(bb);
    if (!ItemAdd(bb, NULL))
      return false;

    const ImGuiID id = window->GetID("surfaces");
    const int hovered = ItemHoverable(bb, id);

    RenderFrame(bb.Min, bb.Max, active_bounds ? ImU32(ImColor(192, 192, 192)) : GetColorU32(ImGuiCol_FrameBg, 1), true, Style.FrameRounding);
    RenderFrameBorder(bb.Min, bb.Max, Style.FrameRounding);
    //if(hovered)
    //  drawList->PushClipRectFullScreen();
    
    auto func_bounds = f.get_bounds();
    float fh = func_bounds[3] - func_bounds[2];
    fh = fh < 133.f ? 133.f : fh;
    func_bounds[2] -= fh * 0.25;
    func_bounds[3] += fh * 0.25;
    ImRect fbb(ImVec2(func_bounds[0],func_bounds[2]), ImVec2(func_bounds[1],func_bounds[3]));
    
    if(drawPlotState.state == DrawPlotState::eState::DRAGFPT) // restore from state - do not update during drag
      fbb = drawPlotState.fbb;
    
    drawPlotState.fbb = fbb;

    if (active_bounds)
      drawList->AddRectFilled(ImVec2(mapPoint(fpt{ (*active_bounds)[0], 0.f }, fbb, bb).x,bb.Min.y+1), ImVec2(mapPoint(fpt{ (*active_bounds)[1], 0.f }, fbb, bb).x, bb.Max.y-1), GetColorU32(ImGuiCol_FrameBg, 1));
    
    drawList->AddLine(ImVec2(bb.Min.x+1, mapPoint(fpt{ 0.f, 0.f }, fbb, bb).y), ImVec2(bb.Max.x - 1, mapPoint(fpt{ 0.f, 0.f }, fbb, bb).y), ImColor(128, 128, 128));

    for(int i=0;i<f.segments_count();i++)
      drawList->AddLine(mapPoint(f[i], fbb, bb), mapPoint(f[i+1], fbb, bb), color, 1);
    
    for(auto& pt : f)
      drawList->AddCircleFilled(mapPoint(pt, fbb, bb), 3, color, 8);
      
    bool modified = false;

    if(hovered) {
      if(!(IsMouseClicked(0) || IsMouseDragging(0)))
        drawPlotState.state = DrawPlotState::eState::HOVER;

      drawPlotState.call_count++;
  
      ImVec2 mouse_pos = GetIO().MousePos;
      ImVec2 pti = mapPoint(std::array<float, 2>{mouse_pos.x, mouse_pos.y}, bb, fbb);
      fpt pos = fpt{pti.x, pti.y};
      int segment = 0;
      int pt_id=0;
      for(; pt_id < f.size(); pt_id++) {
        if(in3pix(mapPoint(f[pt_id], fbb, bb), mouse_pos))
          break;
      }
      if(pt_id == f.size()) {
        for(; segment<f.size()-1;segment++) {
          ImVec2 a = mapPoint(f[segment], fbb, bb);
          ImVec2 b = mapPoint(f[segment+1], fbb, bb);
          if(mouse_pos.x < b.x && mouse_pos.x > a.x) {
            if(in3pix(mouse_pos, a + (mouse_pos.x - a.x) / (b.x - a.x) * (b - a)))
              break;
          }
        }
      }

      if(drawPlotState.state == DrawPlotState::eState::DRAGFPT && drawPlotState.pt_id >= 0 && drawPlotState.pt_id < f.size())
        pt_id = drawPlotState.pt_id;

      if(pt_id != f.size()) {
        drawList->AddCircleFilled(mapPoint(f[pt_id], fbb, bb), 5, ImColor(0,0,0), 13);
        if (IsMouseClicked(0) || IsMouseDragging(0)) { // modify existing point
          drawPlotState.state = DrawPlotState::eState::DRAGFPT;
          drawPlotState.pt_id = pt_id;
          if(pt_id!=0 && pt_id!=f.segments_count())
            f[pt_id][0] = pos[0];
          f[pt_id][1] = pos[1];
          modified = true;
        }
        SetTooltip("(%4.3f, %4.3f)", f[pt_id][0], f[pt_id][1]);
        if(IsMouseClicked(1)) { // remove existing point
          f.erase(f.begin()+pt_id);
          modified = true;
        }
      } else {
        if(segment != f.segments_count()) {
          drawList->AddCircle(mapPoint(pos, fbb, bb), 5, ImColor(0,0,0), 13, 1);
          if(IsMouseClicked(0)) {
            f.insert(f.begin() + segment + 1, pos);
            modified = true;
          }
        }
      }
    }
    
    //if(hovered)
    //  drawList->PopClipRect();
    return modified;
  }
}

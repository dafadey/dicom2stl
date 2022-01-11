#include <string>
#include <iostream>
#include <cmath>

#include "ImGuiDrawPlot.h"
#include "imgui_controls.h"
#include "object.h"
#include "draw.h"

static void _update_camera(const vec3& dir, const vec3& up, renderer* ren) {
  if(dir*dir) {
    ren->cam_pos = ren->fp_pos + std::sqrt((ren->cam_pos - ren->fp_pos) * (ren->cam_pos - ren->fp_pos)) * dir;
    ren->cam_up = up;
  }
}

namespace ImGui {
  
  bool DicomSurfaceControl(object* obj, renderer* ren) {
    static std::map<std::string, DicomSurfaceControlState> dicomSurfaceControlStateMap;

    if(dicomSurfaceControlStateMap.find(obj->name) == dicomSurfaceControlStateMap.end())
      dicomSurfaceControlStateMap[obj->name] = DicomSurfaceControlState();

    auto& dicomSurfaceControlState = dicomSurfaceControlStateMap[obj->name];

    char output_name[113];
    if(obj->output_name.size() < 113) {
      int i=0;
      for(;i < obj->output_name.size(); i++)
        output_name[i] = obj->output_name.c_str()[i];
      output_name[i] = char(0);
    }

    if(InputText("output", output_name, 113))
      obj->output_name = std::string(output_name);
    SameLine();
    if(Button("save"))
      obj->save_surface();
    Checkbox("visible", &(obj->ogl_surface.visible));
    SameLine();
    ColorEdit3("color", obj->ogl_surface.surface_color);
    SliderFloat("shiny", &(obj->ogl_surface.shiny), 0.f, 1.f);
    
    bool recrop_surface = false;
    recrop_surface |= SliderFloat("xc", &(obj->box_center[0]), -obj->source->o[0], static_cast<float>(obj->source->dim[0]) * obj->source->d[0] - obj->source->o[0]);
    recrop_surface |= SliderFloat("yc", &(obj->box_center[1]), -obj->source->o[1], static_cast<float>(obj->source->dim[1]) * obj->source->d[1] - obj->source->o[1]);
    recrop_surface |= SliderFloat("zc", &(obj->box_center[2]), -obj->source->o[2], static_cast<float>(obj->source->dim[2]) * obj->source->d[2] - obj->source->o[2]);
    recrop_surface |= SliderFloat3("size %", obj->box_scale.data(), 0.f , 1.f);
    dicomSurfaceControlState.smoothness_changed |= SliderFloat("smoothness", &(obj->smooth_factor), 0.f , 3.3f);
    //SameLine();

    if (dicomSurfaceControlState.smoothness_changed) {
      SameLine();
      if(Button("smooth")) {
        obj->smooth_surface();
        dicomSurfaceControlState.smoothness_changed = false;
        recrop_surface = true;
      }
    }

    bool remake_surface = false;

    if(recrop_surface) {
      obj->crop_surface();
      remake_surface = true;
    }

    remake_surface |= SliderFloat("level", &(obj->level),0,1700);

    if(obj->cropped) {
      auto& cs = *(obj->cropped);
      std::array<float, 2> xbounds{-cs.o[0], cs.d[0] * static_cast<float>(cs.dim[0]) - cs.o[0]};
      remake_surface |= DrawPlot(obj->levels[0], "x", ImColor(223,0,0), &xbounds);
      std::array<float, 2> ybounds{ -cs.o[1], cs.d[1] * static_cast<float>(cs.dim[1]) - cs.o[1] };
      remake_surface |= DrawPlot(obj->levels[1], "y", ImColor(0,159,0), &ybounds);
      std::array<float, 2> zbounds{ -cs.o[2], cs.d[2] * static_cast<float>(cs.dim[2]) - cs.o[2] };
      remake_surface |= DrawPlot(obj->levels[2], "z", ImColor(0,0,255), &zbounds);
    }

    if(Button("reset mask")) {
      remake_surface |= true;
      obj->reset_mask();
    }
    SameLine();
    if (Button("invert mask")) {
      remake_surface |= true;
      obj->invert_mask();
    }

    if(remake_surface)
      obj->make_surface(ren);

    return true;
  }

  bool CameraControl(renderer* ren) {
    bool modified = false;
    vec3 dir{0,0,0};
    vec3 up;
    if(Button("x+")) {
      modified = true;
      dir = vec3{1,0,0};
      up = vec3{0,1,0};
    }
    SameLine();
    if(Button("x-")) {
      modified = true;
      dir = vec3{-1,0,0};
      up = vec3{0,1,0};
    }
    SameLine();
    if(Button("y+")) {
      modified = true;
      dir = vec3{0,1,0};
      up = vec3{1,0,0};
    }
    SameLine();
    if(Button("y-")) {
      modified = true;
      dir = vec3{0,-1,0};
      up = vec3{1,0,0};
    }
    SameLine();
    if(Button("z+")) {
      modified = true;
      dir = vec3{0,0,1};
      up = vec3{1,0,0};
    }
    SameLine();
    if(Button("z-")) {
      modified = true;
      dir = vec3{0,0,-1};
      up = vec3{1,0,0};
    }

    ::_update_camera(dir, up, ren);

    modified |= Checkbox("parallel", &(ren->parallel));
    static bool mask_symmetric = true;
    Checkbox("mask painter symmetric", &mask_symmetric);
    if (mask_symmetric) {
      if(InputFloat("mask brush radius", &(ren->mask_r_x)))
        ren->mask_r_y = ren->mask_r_x;
    } else {
      InputFloat("mask brush radius x", &(ren->mask_r_x));
      InputFloat("mask brush radius y", &(ren->mask_r_y));
    }
    
    Checkbox("mask brush erases", &(ren->mask_erases));
    
    return modified;
  }
  
}

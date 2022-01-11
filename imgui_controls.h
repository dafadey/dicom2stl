#pragma once

struct renderer;

struct object;

namespace ImGui {
  
  struct DicomSurfaceControlState {
    DicomSurfaceControlState() = default;
    bool smoothness_changed = false;
  };

  bool DicomSurfaceControl(object* obj, renderer* ren);
  
  bool CameraControl(renderer* ren);
}

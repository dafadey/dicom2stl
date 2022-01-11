#include <iostream>
#include <limits>
#include "readDICOM.h"
#include "draw.h"
#include "interface.h"
#include "ImGuiDrawPlot.h"
#include "object.h"
#include "imgui_controls.h"

int main(int argc, char* argv[]) {
  std::cout << "Hallo!\n";


  vec3f box_center{ 0, 0, 0 };
  vec3f box_scale{ 1., 1., 1. };
  float level{std::numeric_limits<float>::max()};
  float smoothness{1.7};

  bool err{false};

  std::string folder;
  std::string series;
  std::string output;
  std::string input_ini;

  if(argc > 2) {
    for (int i = 1; i < argc; i++)
    {
      if (std::string(argv[i]) == "--box") {
        if(i + 6 > argc) {
          err = true;
          break;
        }
        for(int ii=0; ii<3; ii++)
          box_center[ii] = atof(argv[++i]);
        for (int ii = 0; ii < 3; ii++)
          box_scale[ii] = atof(argv[++i]);
      }
      if (std::string(argv[i]) == "--dicom") {
        if (i + 1 > argc) {
          err = true;
          break;
        }
        folder = argv[++i];
      }
      if (std::string(argv[i]) == "--series") {
        if (i + 1 > argc) {
          err = true;
          break;
        }
        series = argv[++i];
      }
      if (std::string(argv[i]) == "--level") {
        if (i + 1 > argc) {
          err = true;
          break;
        }
        level = atof(argv[++i]);
      }
      if (std::string(argv[i]) == "--smooth") {
        if (i + 1 > argc) {
          err = true;
          break;
        }
        smoothness = atof(argv[++i]);
      }
      if (std::string(argv[i]) == "--output") {
        if (i + 1 > argc) {
          err = true;
          break;
        }
        output = argv[++i];
      }
      if (std::string(argv[i]) == "--ini") {
        if (i + 1 > argc) {
          err = true;
          break;
        }
        input_ini = argv[++i];
      }
    }
  } else
    err=true;
  if(level == std::numeric_limits<float>::max() && input_ini.empty())
    err = true;

  if(folder.empty() && input_ini.empty())
    err = true;

  if (err) {
    std::cout << "--dicom <path folder to folder containing DICOM data>\n--series <series name>\n--box <center_x center_y center_z scale_x scale_y scale_z>\n\texample: --box -11 13 17 0.5 0.6 0.7\n\tfull box: --box 0 0 0 1 1 1\n--level <level in center>\n--smooth <smoothness factor>\n\talternativelly set --ini <ini file> to load configuration";
    return 0;
  }

  std::cout << "DICOM folder: " << folder << "\nSeries name: " << series << "\nbox: (" << box_center[0] << ", " << box_center[1] << ", " << box_center[2] << ") scale: (" << box_scale[0] << ", " << box_scale[1] << ", " << box_scale[2] << ")\n";

  std::vector<vgeo> vgs;

  if(!folder.empty()) {
    DICOMreader r(folder.c_str(), series.c_str());
    r.read(vgs);
  }

/*
* 371
* 831
  vec3 box_center{2, 7, -7};
  vec3 box_scale{ 0.371, 0.411, 0.454};
*/

  std::vector<object*> objects;

  if(!input_ini.empty())
    objects.push_back(load_surface(input_ini, vgs));

  for(int vgid=0; vgid< vgs.size(); vgid++) {
    auto& vg = vgs[vgid];
    std::cout << "doing vg #" << vgid << " [" << vg.dim[0] << ", " << vg.dim[1] << ", " << vg.dim[2] << "]\n";
    //vgeo vg_cropped = crop(vg, box_center, box_scale);
    //vgeo svg = smooth(vg_cropped, 1.7f);
    //geo g = vgeo2surface(svg, level, level_side);
    //std::cout << "writing " << g.size() << " triangles, with " << g.points.size() << " points\n";
    //std::string outFile = output+"_"+std::to_string(vgid) + ".stl";
    //saveSTL(outFile.c_str(), g);
    
    imgui_interface iface;
    iface.init();
    
    renderer ren;
    
    ren.init(iface.window, &objects);

    if(objects.size() == 0) {
      object* obj = new object(&vg);
      obj->name = "input";
      obj->output_name = obj->name;
      obj->source_name = folder;
      obj->box_center = box_center;
      obj->box_scale = box_scale;
      obj->level=level;
      obj->smooth_factor = smoothness;
      objects.push_back(obj);
    }
    objects[0]->smooth_surface();
    objects[0]->crop_surface();
    objects[0]->make_surface(&ren);

    ren.reset_camera();

    glfwSetWindowUserPointer(iface.window, (void*) &ren);
    
    ren.set_callbacks(iface.window);
    
    char new_surface_name[113]="surf";

    while (!glfwWindowShouldClose(iface.window))  {
      glfwWaitEvents();
      //glfwPollEvents();
      ren.nocallbacks = false;

      static int counter{0};

      ImGui_ImplOpenGL3_NewFrame();
      ImGui_ImplGlfw_NewFrame();
      ImGui::NewFrame();

      ImGui::Begin("surfaces");
      ImGui::Text(("source: " + folder).c_str());
      ImGui::InputText("name", new_surface_name, 113);
      ImGui::SameLine();
      if(ImGui::Button("+")) {
        std::string name = std::string(new_surface_name);
        bool available = true;
        for(auto obj : objects) {
          if(obj->name == name)
            available = false;
        }
        if(available) {
          object* obj = new object(&vg);
          obj->name = name;
          obj->output_name = obj->name;
          obj->source_name = folder;
          obj->box_center = box_center;
          obj->box_scale = box_scale;
          obj->level = level;
          objects.push_back(obj);
          obj->smooth_surface();
          obj->crop_surface();
          obj->make_surface(&ren);

        }
      }

      object* obj_to_remove = nullptr;
      for(auto& obj : objects) {
        if(ImGui::TreeNode(obj->name.c_str())) {
          ImGui::DicomSurfaceControl(obj, &ren);
          if(ImGui::Button("x"))
            obj_to_remove = obj;
          ImGui::TreePop();
        }
      }
      
      for(int i=0;i<objects.size();i++) {
        if(objects[i] == obj_to_remove) {
          objects.erase(objects.begin() + i);
          break;
        }
      }
      ren.remove_item(&(obj_to_remove->ogl_surface));
      delete obj_to_remove;
      
      //ImGui::Text("%.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
      ImGui::End();

      ImGui::Begin("camera");
      ImGui::CameraControl(&ren);
      ImGui::Text("%.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
      int memory_total=0;
      for(auto& obj : objects)
        memory_total += obj->memory();
      ImGui::Text("Memory: %d Mb", memory_total/1024/1024);
      ImGui::End();

      ImGui::Render();

      if (ImGui::GetIO().WantCaptureMouse)
        ren.nocallbacks = true;
      glClearColor(.3, .3, .3, 1.);
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
      ren.render();
      ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
      glfwSwapBuffers(iface.window);
      glFlush();
    }
    
    iface.close();
    break;
  }
}

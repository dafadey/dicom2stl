EXE = dicom2stl
#BASE_PATH = /home/dan
IMGUI_DIR = $(BASE_PATH)/imgui
IMPLOT_DIR = $(BASE_PATH)/implot
ITK_DIR = $(BASE_PATH)/ITK_5.2.1/install
SOURCES = main.cpp readDICOM.cpp vgeo.cpp geo.cpp saveSTL.cpp interface.cpp draw.cpp OGLitem.cpp func.cpp ImGuiDrawPlot.cpp object.cpp imgui_controls.cpp timer.cpp base64/base64.cpp
SOURCES += $(IMGUI_DIR)/imgui.cpp $(IMGUI_DIR)/imgui_draw.cpp $(IMGUI_DIR)/imgui_tables.cpp $(IMGUI_DIR)/imgui_widgets.cpp
SOURCES += $(IMGUI_DIR)/backends/imgui_impl_glfw.cpp $(IMGUI_DIR)/backends/imgui_impl_opengl3.cpp
OBJS = $(addsuffix .o, $(basename $(notdir $(SOURCES))))
UNAME_S := $(shell uname -s)
LINUX_GL_LIBS = -lGL -lGLEW -fopenmp -lpthread

CXXFLAGS = -I$(ITK_DIR)/include/ITK-5.2/ -I$(IMGUI_DIR) -I$(IMGUI_DIR)/backends
CXXFLAGS += -g -O3 -fopenmp -DNOIMPLOT --fast-math -march=native -mtune=native
LIBS = $(ITK_DIR)/lib/libITKIOBMP-5.2.a $(ITK_DIR)/lib/libITKIOGDCM-5.2.a $(ITK_DIR)/lib/libITKIOGIPL-5.2.a $(ITK_DIR)/lib/libITKIOJPEG-5.2.a $(ITK_DIR)/lib/libITKIOMeshBYU-5.2.a $(ITK_DIR)/lib/libITKIOMeshFreeSurfer-5.2.a $(ITK_DIR)/lib/libITKIOMeshGifti-5.2.a $(ITK_DIR)/lib/libITKIOMeshOBJ-5.2.a $(ITK_DIR)/lib/libITKIOMeshOFF-5.2.a $(ITK_DIR)/lib/libITKIOMeshVTK-5.2.a $(ITK_DIR)/lib/libITKIOMeta-5.2.a $(ITK_DIR)/lib/libITKIONIFTI-5.2.a $(ITK_DIR)/lib/libITKIONRRD-5.2.a $(ITK_DIR)/lib/libITKIOPNG-5.2.a $(ITK_DIR)/lib/libITKIOTIFF-5.2.a $(ITK_DIR)/lib/libITKIOVTK-5.2.a $(ITK_DIR)/lib/libITKCommon-5.2.a $(ITK_DIR)/lib/libITKIOImageBase-5.2.a $(ITK_DIR)/lib/libITKTestKernel-5.2.a $(ITK_DIR)/lib/libITKIOBMP-5.2.a $(ITK_DIR)/lib/libITKIOGDCM-5.2.a $(ITK_DIR)/lib/libitkgdcmMSFF-5.2.a $(ITK_DIR)/lib/libitkgdcmDICT-5.2.a $(ITK_DIR)/lib/libitkgdcmIOD-5.2.a $(ITK_DIR)/lib/libitkgdcmDSED-5.2.a $(ITK_DIR)/lib/libitkgdcmCommon-5.2.a $(ITK_DIR)/lib/libitkgdcmjpeg8-5.2.a $(ITK_DIR)/lib/libitkgdcmjpeg12-5.2.a $(ITK_DIR)/lib/libitkgdcmjpeg16-5.2.a $(ITK_DIR)/lib/libitkgdcmopenjp2-5.2.a $(ITK_DIR)/lib/libitkgdcmcharls-5.2.a $(ITK_DIR)/lib/libitkgdcmuuid-5.2.a $(ITK_DIR)/lib/libITKIOGIPL-5.2.a $(ITK_DIR)/lib/libITKIOJPEG-5.2.a $(ITK_DIR)/lib/libITKIOMeshBYU-5.2.a $(ITK_DIR)/lib/libITKIOMeshFreeSurfer-5.2.a $(ITK_DIR)/lib/libITKIOMeshGifti-5.2.a $(ITK_DIR)/lib/libITKgiftiio-5.2.a $(ITK_DIR)/lib/libITKEXPAT-5.2.a $(ITK_DIR)/lib/libITKIOMeshOBJ-5.2.a $(ITK_DIR)/lib/libITKIOMeshOFF-5.2.a $(ITK_DIR)/lib/libITKIOMeshVTK-5.2.a $(ITK_DIR)/lib/libITKIOMeshBase-5.2.a $(ITK_DIR)/lib/libITKQuadEdgeMesh-5.2.a $(ITK_DIR)/lib/libITKMesh-5.2.a $(ITK_DIR)/lib/libITKIOMeta-5.2.a $(ITK_DIR)/lib/libITKMetaIO-5.2.a $(ITK_DIR)/lib/libITKIONIFTI-5.2.a $(ITK_DIR)/lib/libITKniftiio-5.2.a $(ITK_DIR)/lib/libITKznz-5.2.a $(ITK_DIR)/lib/libITKTransform-5.2.a $(ITK_DIR)/lib/libITKIONRRD-5.2.a $(ITK_DIR)/lib/libITKNrrdIO-5.2.a $(ITK_DIR)/lib/libITKIOPNG-5.2.a $(ITK_DIR)/lib/libitkpng-5.2.a $(ITK_DIR)/lib/libITKIOTIFF-5.2.a $(ITK_DIR)/lib/libitktiff-5.2.a $(ITK_DIR)/lib/libitkzlib-5.2.a $(ITK_DIR)/lib/libitkjpeg-5.2.a $(ITK_DIR)/lib/libITKIOVTK-5.2.a $(ITK_DIR)/lib/libITKIOImageBase-5.2.a $(ITK_DIR)/lib/libITKCommon-5.2.a $(ITK_DIR)/lib/libitkdouble-conversion-5.2.a $(ITK_DIR)/lib/libitksys-5.2.a $(ITK_DIR)/lib/libITKVNLInstantiation-5.2.a $(ITK_DIR)/lib/libitkvnl_algo-5.2.a $(ITK_DIR)/lib/libitkvnl-5.2.a $(ITK_DIR)/lib/libitkv3p_netlib-5.2.a $(ITK_DIR)/lib/libitkvcl-5.2.a

ifeq ($(UNAME_S), Linux) #LINUX
	ECHO_MESSAGE = "Linux"
	LIBS += $(LINUX_GL_LIBS) `pkg-config --static --libs glfw3`

	CXXFLAGS += `pkg-config --cflags glfw3`
	CFLAGS = $(CXXFLAGS)
endif

ifeq ($(UNAME_S), Darwin) #APPLE
	ECHO_MESSAGE = "Mac OS X"
	LIBS += -framework OpenGL -framework Cocoa -framework IOKit -framework CoreVideo
	LIBS += -L/usr/local/lib -L/opt/local/lib -L/opt/homebrew/lib
	#LIBS += -lglfw3
	LIBS += -lglfw

	CXXFLAGS += -I/usr/local/include -I/opt/local/include -I/opt/homebrew/include
	CFLAGS = $(CXXFLAGS)
endif

ifeq ($(OS), Windows_NT)
	ECHO_MESSAGE = "MinGW"
	LIBS += -lglfw3 -lgdi32 -lopengl32 -limm32

	CXXFLAGS += `pkg-config --cflags glfw3`
	CFLAGS = $(CXXFLAGS)
endif

##---------------------------------------------------------------------
## BUILD RULES
##---------------------------------------------------------------------

%.o:%.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

%.o:$(IMGUI_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

%.o:$(IMGUI_DIR)/backends/%.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

%.o:base64/%.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<


all: $(EXE)
	@echo Build complete for $(ECHO_MESSAGE)

$(EXE): $(OBJS)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LIBS)

clean:
	rm -f $(EXE) $(OBJS)


cmake_minimum_required(VERSION 2.8)

PROJECT(MitoGraph)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(MitoGraph MACOSX_BUNDLE MitoGraph MitoThinning ssThinning)

if(VTK_LIBRARIES)
  target_link_libraries(MitoGraph ${VTK_LIBRARIES})
else()
  target_link_libraries(MitoGraph vtkHybrid vtkWidgets)
endif()

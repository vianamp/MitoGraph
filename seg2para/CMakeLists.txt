cmake_minimum_required(VERSION 2.8)

PROJECT(seg2para)

include_directories('.')

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(seg2para MACOSX_BUNDLE seg2para)

if(VTK_LIBRARIES)
  target_link_libraries(seg2para ${VTK_LIBRARIES})
else()
  target_link_libraries(seg2para vtkHybrid vtkWidgets)
endif()

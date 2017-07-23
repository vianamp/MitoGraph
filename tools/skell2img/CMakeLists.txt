
cmake_minimum_required(VERSION 2.8)

PROJECT(skell2img)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(skell2img MACOSX_BUNDLE skell2img)

if(VTK_LIBRARIES)
  target_link_libraries(skell2img ${VTK_LIBRARIES})
else()
  target_link_libraries(skell2img vtkHybrid vtkWidgets)
endif()

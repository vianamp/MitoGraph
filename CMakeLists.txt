
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(MitoGraph)

FIND_PACKAGE(VTK)
INCLUDE(${VTK_USE_FILE})
INCLUDE_DIRECTORIES(${VTK_INCLUDE_DIRS})

ADD_EXECUTABLE(MitoGraph MitoGraph MitoThinning ssThinning)

if(VTK_LIBRARIES)
  TARGET_LINK_LIBRARIES(MitoGraph ${VTK_LIBRARIES})
else()
  TARGET_LINK_LIBRARIES(MitoGraph vtkHybrid vtkWidgets)
endif()


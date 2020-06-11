CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(MitoGraph)

# Ensure VTK is present
FIND_PACKAGE(VTK)
if(VTK_FOUND)
  message (STATUS "VTK Version Found: ${VTK_VERSION}")
  INCLUDE(${VTK_USE_FILE})
else()
  message (STATUS "VTK was not found!")
  return()
endif()

# Specify the executables to produce
ADD_EXECUTABLE(MitoGraph MitoGraph MitoThinning ssThinning)

# Link VTK components
if(VTK_LIBRARIES)
  TARGET_LINK_LIBRARIES(MitoGraph ${VTK_LIBRARIES})
else()
  TARGET_LINK_LIBRARIES(MitoGraph vtkHybrid vtkWidgets)
endif()

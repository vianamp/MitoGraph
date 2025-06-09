CMAKE_MINIMUM_REQUIRED(VERSION 3.12)

PROJECT(MitoGraph)

# Find VTK installation
FIND_PACKAGE(VTK REQUIRED)

# VTK 9 uses modules, specify the components we need
IF(VTK_VERSION VERSION_GREATER_EQUAL "9.0")
    FIND_PACKAGE(VTK REQUIRED COMPONENTS
        CommonCore
        CommonDataModel
        CommonExecutionModel
        FiltersCore
        FiltersGeneral
        FiltersGeometry
        FiltersModeling
        FiltersSources
        IOCore
        IOImage
        IOLegacy
        IOXML
        ImagingCore
        ImagingGeneral
        ImagingHybrid
        ImagingMath
        InteractionStyle
        RenderingCore
        RenderingFreeType
        RenderingOpenGL2
    )
ENDIF()

# Include VTK use file for older versions
IF(VTK_VERSION VERSION_LESS "8.90.0")
    INCLUDE(${VTK_USE_FILE})
ENDIF()

# Include directories
INCLUDE_DIRECTORIES(${VTK_INCLUDE_DIRS})

# Add executable
ADD_EXECUTABLE(MitoGraph MitoGraph MitoThinning ssThinning)

# Link libraries
IF(VTK_VERSION VERSION_GREATER_EQUAL "8.90.0")
    # VTK 9.x
    VTK_MODULE_AUTOINIT(
        TARGETS MitoGraph
        MODULES ${VTK_LIBRARIES}
    )
    TARGET_LINK_LIBRARIES(MitoGraph ${VTK_LIBRARIES})
ELSE()
    # VTK 8.x and older
    IF(VTK_LIBRARIES)
        TARGET_LINK_LIBRARIES(MitoGraph ${VTK_LIBRARIES})
    ELSE()
        TARGET_LINK_LIBRARIES(MitoGraph vtkHybrid vtkWidgets)
    ENDIF()
ENDIF()

# Set C++ standard
SET_PROPERTY(TARGET MitoGraph PROPERTY CXX_STANDARD 11)
SET_PROPERTY(TARGET MitoGraph PROPERTY CXX_STANDARD_REQUIRED ON)


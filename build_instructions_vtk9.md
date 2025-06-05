# MitoGraph VTK 9 Compatibility Build Instructions

## Overview

This document describes how to build MitoGraph 3.1 with VTK 9 on macOS. The original version requires VTK 8.2.0, but through modifications to CMakeLists.txt and API calls in the source code, it is now compatible with modern VTK 9.

## System Requirements

- macOS 10.15+
- Homebrew
- CMake 3.12+
- VTK 9.x (installed via Homebrew)
- C++11 compatible compiler

## Installation Steps

### 1. Install Dependencies

```bash
# Install CMake and VTK
brew install cmake vtk

# Verify VTK version
brew list vtk --versions
```

### 2. Build MitoGraph

```bash
# Navigate to project directory
cd /path/to/MitoGraph-3.1

# Create and enter build directory
mkdir -p build && cd build

# Configure project
cmake ..

# Compile
make
```

### 3. Test Installation

```bash
# Test basic functionality
./MitoGraph -xy 0.0645 -z 0.2 -path . -export_image_binary
```

## Major Modifications

### CMakeLists.txt Updates

- Upgraded minimum CMake version requirement to 3.12
- Added VTK 9 modular support
- Maintained backward compatibility with VTK 8.x
- Set C++11 standard

### Source Code API Fixes

1. **FindPoint API Change**: VTK 9 changed `FindPoint(x,y,z)` to `FindPoint(point[3])`
2. **vtkLongArray Deprecation**: Replaced with `vtkTypeInt64Array`
3. **register Keyword**: Removed deprecated register storage class specifier

### Compatibility

- Supports VTK 8.2.0+
- Supports VTK 9.x
- Backward compatible with original functionality
- Maintains same command-line interface

## Troubleshooting

### Common Issues

1. **VTK Not Found**: Ensure VTK is properly installed via Homebrew
2. **CMake Version Too Old**: Upgrade to CMake 3.12+
3. **Compilation Warnings**: sprintf warnings do not affect functionality and can be ignored

### Verify Installation

```bash
# Check executable
ls -la build/MitoGraph

# Run help command
./build/MitoGraph
```

## Performance Notes

- Compiled version has identical functionality to original VTK 8.2.0 version
- Performance characteristics remain consistent
- Supports all original command-line arguments and features

## License

This modified version follows the original MitoGraph license terms. 
## Recipe for MacOS 15.5 Chipset M2 and VTK9

```
brew install cmake vtk
git clone https://github.com/vianamp/MitoGraph.git && cd MitoGraph
mkdir build && cd build
cmake ..
make
```

## Recipe for MacOS 12.5.1 Chipset M2 and VTK8

```
mkdir mitograph
cd mitograph
wget https://www.vtk.org/files/release/8.2/VTK-8.2.0.tar.gz
tar -xvf VTK-8.2.0.tar.gz && cd VTK-8.2.0
mkdir build && cd build
brew install cmake
brew install libpng
export CFLAGS="-DPNG_ARM_NEON_OPT=0"
export CXXFLAGS="-DPNG_ARM_NEON_OPT=0"
cmake -DPNG_PNG_INCLUDE_DIR=/opt/homebrew/include -DPNG_LIBRARY_RELEASE=/opt/homebrew/lib/libpng.dylib ..
make
sudo make install
cd ../..
git clone https://github.com/vianamp/MitoGraph.git && cd MitoGraph
mkdir build && cd build
cmake ..
make
```

## Recipe for Windows 10 and VTK8

1. Install CMake
2. Install Visual Studio 17 2022
3. Download VTK 8.2.0 Source Code
4. Open CMake (cmake-gui).
5. Set the source directory to the VTK source directory.
6. Set the build directory to where you want to build the binaries (e.g., C:/vtk-build).
7. Click Configure and select your Visual Studio version
8. Click Generate and then Open Project.
9. In Visual Studio, build the project (ALL_BUILD).


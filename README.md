<p align="center">
  <img src="doc/mitograph.png" width="auto" height="380" title="MoCo Logo">
</p>

MitoGraph is a fully automated image processing method and software dedicated to calculating the three-dimensional morphology of mitochondria in live cells. MitoGraph is currently optimized and validated only for quantifying the volume and topology of tubular mitochondrial networks in budding yeast [1,2]. However, MitoGraph can also be applied to mitochondria in other cell types and possibly other intracellular (or tissue) structures, with proper validation. MitoGraph is continuously being updated. Please contact us if you have questions that go beyond those discussed in [1,2] or any request that would make MitoGraph a better tool for your research.

Matheus Viana - vianamp@gmail.com - 01.19.2018

### <a href="https://github.com/vianamp/MitoGraph/releases/tag/v3.0">v3.0 Release</a>

✓ Adaptive threshold

✓ Support binary and VTK input

✓ Support 2D image data

✓ Export skeleton coordinates as text file

✓ Save skeleton with information of mitochondrial tubules width, pixel intensity, edge length and nodes position

### <a href="https://github.com/vianamp/MitoGraph/releases/tag/v2.1">v2.1 Release</a>

✓ Extra flag that lets you use __R__ to perform basic analysis of MitoGraph output

---

### How to Install

Please follow the steps bellow to install MitoGraph.

1. <a href="https://github.com/vianamp/MitoGraph/releases/download/v3.0/MitoGraph-OSX-10.12.zip">Download the latest version of MitoGraph available for Mac OSX 10.10 or later</a>

2. Create a folder in your Desktop and name it MitoGraph

3. Move the downloaded file to the folder MitoGraph in your and unzip it

---

### Testing Data, ImageJ Macros and Other Tools

Before running MitoGraph on your own dataset, you may want to check the tools we made available for preparing your data and our test datasets.

<p align="center">
  <img src="https://sites.google.com/site/vianamp/_/rsrc/1418664353567/mitograph/mitoexamples.png" width="auto" height="128" title="Example Dataset">
</p>

* <a href="https://github.com/vianamp/MitoGraphTools/blob/master/README.md">MitoGraph Tools</a>: in this repository you will find an example dataset to test MitoGraph on. There is also two ImageJ scripts that help you to prepare your data to make them similar to our test dataset.

---

### How to Run MitoGraph

Download the <a href="https://github.com/vianamp/MitoGraphTools/blob/master/README.md">test dataset</a> mentioned in the previous section and unzip the downloaded file in your Desktop. Type the following command in the terminal of your Mac OS (spotlight + terminal)

`cd ~/Desktop/MitoGraph`

`./MitoGraph -xy 0.056 -z 0.2 -path ~/Desktop/MitoGraphTools-1.0`

The flag `-xy` specifies the pixel size in microns and the flag `-z` specifies the z-spacing also in microns of your images. Finally, the flag `-path` indicates the folder that contains the images you want to analyze.

#### Optional Flags

`-scales a b c` where a, b and c are numeric values specifying the initial, final and total number of scales which should be used by MitoGraph [1]. For example: `-scales 1.5 2.0 6`. Default values are a=1.0, b=1.5 and c=6.

`-threshold a` where a is numeric value in the range [0,1] for the post-divergence threshold. For example: `-threshold 0.1`. Default value is a=0.1667.

`-adaptive a` where a is a numeric integer value specifying that the input image should be split into a x a blocks before the segmentation. This is useful for images with high variablity of pixel intensity.

`-binary`: indicates that the input is a binary image.

`-vtk`: indicates that the input data is of __VTK Imagedata__ type instead of TIFF.

`-labels_off`: turns off the labels of nodes in the output file __filename_nodes.vtk__.

`-analyze`: uses R to do a per connected component analysis of MitoGraph output. This flag requires __R__ and its package __igraph__ to be installed in your system.

---

### MitoGraph Outputs

The output of MitoGraph will be saved in the directory specified with `-path`.

**Single files:**

* __mitograph.config__ - Stores the parameters used to run MitoGraph and the date and time when the analysis was complete.

**Text files that can be open in any Text Editor (one file per sample):**

* __filename.gnet__ - Connection list of type `node_i` `node_j` of the graph that represents the mitochondria. First line of this file gives the total number of nodes.
* __filename.coo__ - Coordinates xyz of nodes in the graph.
* __filename.mitograph__ - Mitochondrial attributes: volume from voxels, average width (µm), std width (µm),  total length (µm) and volume from length (µm3). This file also shows per component statistics when the flag `-analyze` is used.
* __filename.txt__ - Coordinates of all the points along the mitochondrial skeleton as well as the local width (µm) and original pixel intensity.
* __filename.cc__ - Connected component which each node belong to and the component volume (µm3). Only when the flag `-analyze` is used.

**Image files (one file per sample):**

* __filename.png__ - Max projection of mitochondria after MitoGraph binarization. Should be used for fast assessment of MitoGraph segmentation result.

**VTK files that can be open using Paraview [3] (one file per sample):**

* __filename_Nodes.vtk__ - Nodes of the graph that represents the mitochondria and their labels.
* __filename_skeleton.vtk__ - Mitochondrial skeleton. This file contains information about mitochondrial tubule width, pixel intensity, edge length and nodes position that can be viewd in Paraview by setting the coloring mode.
* __filename_mitosurface.vtk__ - Mitochondrial surface.

### Building MitoGraph from Source
#### Windows
:construction: TODO :construction:

___Temporary Solution___: _Windows users can install a linux based OS in a virtual machine to build and use the program from there._

#### MacOS
:construction: TODO :construction:

#### Linux
1. Install [CMake](https://cmake.org/download/) and [VTK 7+](https://vtk.org/download/) if they are not already present on your system. For Debian/Ubuntu based systems, these can be obtained using the following command:</br>
<code> sudo apt install cmake vtklib7-dev </code>

2. Download the repository from GitHub and unzip the folder.

3. Create a build directory in the unzipped MitoGraph directory as follows (assuming you have navigated into the MitoGraph directory):</br>
<code>mkdir build && cd build</code>

5. Generate the build files with CMake with the following command:</br>
<code>cmake -DCMAKE_BUILD_TYPE=Release ../</code>

5. Build the executables with make with the following command.</br>
<code>make</code>

6. Now when you list the files (ls), you should see the MitoGraph executable. You can add this build folder to your PATH variable if you want to and then execute the program from anywhere.

---

### Contributions to MitoGraph

* With help of Hill Lab, we also provide an example dataset of mitochondria in mammalian cells and R scripts that will help you to analyze the data generated by MitoGraph. Please, check this out in <a href="https://github.com/Hill-Lab/MitoGraph-Contrib-RScripts">MitoGraph-Contrib-RScripts</a>.

---

### References

[1] - Matheus P Viana, Swee Lim, Susanne M Rafelski, Quantifying Mitochondrial Content in Living Cells (2015), Biophysical Methods in Cell Biology, (125) - 77-93 (http://www.sciencedirect.com/science/article/pii/S0091679X14000041)

[2] - Under revision.

[3] - https://www.paraview.org/

#ifndef MITOTHINNING_H
#define MITOTHINNING_H

#include <list>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <vtkMath.h>
#include <vtkPolyLine.h>
#include <vtkCellArray.h>
#include <vtkLongArray.h>
#include <vtkImageData.h>
#include <vtkDataArray.h>
#include <vtkTransform.h>
#include <vtkDataObject.h>
#include <vtkFloatArray.h>
#include <vtkVectorText.h>
#include <vtkInformation.h>
#include <vtkSphereSource.h>
#include <vtkSmartPointer.h>
#include <vtkImageResample.h>
#include <vtkStructuredPoints.h>
#include <vtkInformationVector.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnsignedShortArray.h>
#include <vtkStructuredPointsWriter.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkImageExtractComponents.h>
#include <vtkStructuredPointsReader.h>
#include <vtkImageGaussianSmooth.h>
#include <vtkAppendPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkContourFilter.h>
#include <vtkDoubleArray.h>
#include <vtkTIFFReader.h>
#include <vtkTIFFWriter.h>
#include <vtkPointData.h>
#include <vtkImageRFFT.h>
#include <vtkImageCast.h>
#include <vtkPoints.h>

#include "ssThinning.h"

    extern double _rad;
    extern double _dxy, _dz;
    extern bool _export_graph_files;
    extern bool _scale_polydata_before_save;
    extern bool _export_nodes_label;
    extern double _div_threshold;
    extern bool _improve_skeleton_quality;

    extern int ssdx_sort[26];
    extern int ssdy_sort[26];
    extern int ssdz_sort[26];

	// Thinning algorithm based on the paper: "A 3D 6-subinteration thinning
	// algorithm for extracting medial lines", by Kálman Palágyi and Attila
	// Kuba.
	vtkSmartPointer<vtkPolyData> Thinning3D(vtkSmartPointer<vtkImageData> ImageData, const char FileName[], double *attributes);

	// Routine used to save an ImageData
	void SaveImageData(vtkSmartPointer<vtkImageData> Image, const char FileName[], bool _resample = false);

	// Routine used to save a PolyData
	void SavePolyData(vtkSmartPointer<vtkPolyData> PolyData, const char FileName[], bool scale = _scale_polydata_before_save);

	// Calculate the length of a given edge.
	double GetEdgeLength(vtkIdType edge, vtkSmartPointer<vtkPolyData> PolyData);

	// Label connected components in Image. Results are stored
	// in Volume as negative labels. The routine returns the total
	// number of connected components.
	long int LabelConnectedComponents(vtkSmartPointer<vtkImageData> ImageData, vtkSmartPointer<vtkDataArray> Volume, std::vector<long int> &CSz, int ngbh, double threshold);

	// Routine to delete all voxels located at boundaries of the
	// image volume. Otherwise the algorithm will have problems
	// checking the 3D neighborhood of those voxels.
	void CleanImageBoundaries(vtkSmartPointer<vtkImageData> ImageData);

#endif

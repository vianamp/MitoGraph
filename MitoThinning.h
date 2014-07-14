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

	// Thinning algorithm based on the paper: "A 3D 6-subinteration thinning
	// algorithm for extracting medial lines", by Kálman Palágyi and Attila
	// Kuba.
	int Thinning3D(vtkImageData *ImageData, const char FileName[], double *attributes);

	// Routine used to save an ImageData
	void SaveImageData(vtkImageData *Image, const char FileName[]);

	// Routine used to save a PolyData
	void SavePolyData(vtkPolyData *PolyData, const char FileName[], bool scale = _scale_polydata_before_save);

#endif

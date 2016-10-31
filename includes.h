#include <list>
#include <cmath>
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <dirent.h>

#include <vtkMath.h>
#include <vtkImageFlip.h>
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
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkImageExtractComponents.h>
#include <vtkStructuredPointsReader.h>
#include <vtkImageGaussianSmooth.h>
#include <vtkKdTreePointLocator.h>
#include <vtkAppendPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkContourFilter.h>
#include <vtkCleanPolyData.h>
#include <vtkDoubleArray.h>
#include <vtkTIFFReader.h>
#include <vtkImageShiftScale.h>
#include <vtkTIFFWriter.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkImageRFFT.h>
#include <vtkImageCast.h>
#include <vtkPoints.h>

#ifndef _MITOGRAPH_ENV_VARS

	#define _MITOGRAPH_ENV_VARS

	struct attribute { std::string name; double value; };

	struct _mitoObject {
		bool _binary_input;
	    std::string Type, FileName;
	    double Ox, Oy, Oz;
	    double _sigmai, _sigmaf, _dsigma;
	    std::vector<attribute> attributes;
	};

#endif

#define DEBUG
#include "ssThinning.h"
#include "MitoThinning.h"

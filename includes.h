#include <list>
#include <cmath>
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#ifdef _WIN32
#include "includes/dirent.h"
#else
#include <dirent.h>
#endif

#include <vtkMath.h>
#include <vtkImageFlip.h>
#include <vtkPolyLine.h>
#include <vtkCellArray.h>
#include <vtkTypeInt64Array.h>
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
#include <vtkImageMedian3D.h>
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
	    std::string Type;
	    std::string Folder;
	    std::string FileName;
		bool _analyze;
		bool _binary_input;
	    double Ox;
	    double Oy;
	    double Oz;
	    double _sigmai;
	    double _sigmaf;
	    double _dsigma;
	    int _nsigma;
		bool _adaptive_threshold;
    	int _nblks;

	    std::vector<attribute> attributes;
	};

#endif

// #define DEBUG
#include "ssThinning.h"
#include "MitoThinning.h"
